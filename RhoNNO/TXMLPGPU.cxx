// TXMLPGPU
//
// Float-shadow-trained MLP whose Recall() runs on AMD GPU via HIP + hipBLAS.
// Layer forward: y = ReLU?(W x + b) executed as
//   1) copy bias vector into y device buffer
//   2) hipblasSgemm: y += W x   (alpha=1, beta=1)
//   3) optional elementwise ReLU kernel
//
// Part of the Neural Network Objects package (NNO)

#include "TXMLPGPU.h"
#include "VNeuralNetPlotter.h"
#include "TPerceptron.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

#ifdef RHONNO_WITH_GPU
#include <hip/hip_runtime.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {
inline double relu(double x) { return x > 0.0 ? x : 0.0; }
inline double relu_deriv(double y) { return y > 0.0 ? 1.0 : 0.0; }
}

// ----- GPU state (Pimpl) ----------------------------------------------------

struct TXMLPGPUGpuState {
#ifdef RHONNO_WITH_GPU
    // device buffers per layer
    std::vector<float*> d_W;      // [out*in] per layer
    std::vector<float*> d_B;      // [out] per layer
    std::vector<float*> d_Y;      // [out] per layer (gemm output)
    float* d_cur = nullptr;       // current activation buffer
    int    d_cur_cap = 0;         // allocated size in floats
    std::vector<float>  h_x;      // host staging

    void free_all() {
        for (float* p : d_W) if (p) hipFree(p);
        for (float* p : d_B) if (p) hipFree(p);
        for (float* p : d_Y) if (p) hipFree(p);
        d_W.clear(); d_B.clear(); d_Y.clear();
        if (d_cur) hipFree(d_cur);
        d_cur = nullptr; d_cur_cap = 0;
    }
    ~TXMLPGPUGpuState() { free_all(); }
#endif
    bool weights_dirty = true;    // re-upload at next Recall
    bool gpu_usable = false;
};

#ifdef RHONNO_WITH_GPU
namespace {
// y[o] = b[o] + sum_i  W[o*in+i] * x[i]; optional ReLU after.
// Uses one thread per output neuron. For these layer sizes (up to a few
// hundred outputs), simple parallelism beats any library GEMV that may
// or may not have kernels compiled for the local GPU arch.
__global__ void k_gemv(const float* W, const float* x, const float* b,
                       float* y, int in, int out, int do_relu) {
    int o = blockIdx.x * blockDim.x + threadIdx.x;
    if (o < out) {
        float s = b[o];
        const float* wr = W + (size_t)o * in;
        for (int i = 0; i < in; ++i) s += wr[i] * x[i];
        if (do_relu) s = fmaxf(s, 0.0f);
        y[o] = s;
    }
}

// Elementwise ReLU over an in-place float buffer
__global__ void k_relu(float* y, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) y[i] = fmaxf(y[i], 0.0f);
}

// Fill a buffer with a constant vector (bias)
__global__ void k_copy(const float* src, float* dst, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) dst[i] = src[i];
}
}
#endif

// ----- ctor/dtor ------------------------------------------------------------

void TXMLPGPU::BuildLayers(int layers, int innodes,
                           const std::vector<int>& nodes,
                           const std::vector<double>& steps,
                           bool reluOutput)
{
    if (layers < 1) Errorf((char*)"(TXMLPGPU) at least one layer necessary");
    if ((int)nodes.size() < layers) Errorf((char*)"(TXMLPGPU) nodes size < layers");
    if ((int)steps.size() < layers) Errorf((char*)"(TXMLPGPU) learnSteps size < layers");

    fParm.fLayers = layers;
    fParm.fInNodes = innodes;
    fParm.fOutNodes = nodes[layers - 1];
    fParm.fTransferId = TNeuralNetParameters::TR_RELU;
    fReLUOutput = reluOutput;

    fLayers.clear();
    fLayers.resize(layers);
    int prev = innodes;
    for (int L = 0; L < layers; ++L) {
        LayerF& ly = fLayers[L];
        ly.fInNodes = prev;
        ly.fOutNodes = nodes[L];
        ly.fLearnStep = steps[L];
        ly.fReLU = (L < layers - 1) ? true : reluOutput;
        ly.fWFloat.assign(static_cast<size_t>(ly.fOutNodes) * ly.fInNodes, 0.0);
        ly.fBFloat.assign(ly.fOutNodes, 0.0);
        prev = ly.fOutNodes;
    }

    fOut = new double[fParm.fOutNodes];
    TestPointer(fOut);
}

TXMLPGPU::TXMLPGPU()
{
    // out-of-line so unique_ptr<TXMLPGPUGpuState> sees a complete type
}

TXMLPGPU::TXMLPGPU(int layers, double inputRange, std::string netFile,
                   int innodes, const std::vector<int>& nodes,
                   const std::vector<double>& learnSteps, bool reluOutput)
    : VSupervisedNet("XMLPGPU", innodes, 0, netFile)
{
    fParm.fInScale = 1.0 / inputRange;
    BuildLayers(layers, innodes, nodes, learnSteps, reluOutput);
    InitNet();
}

TXMLPGPU::TXMLPGPU(int layers, double inputRange, std::string netFile,
                   int innodes, int n0, int n1, int n2,
                   double s0, double s1, double s2, bool reluOutput)
    : VSupervisedNet("XMLPGPU", innodes, 0, netFile)
{
    if (layers != 3) Errorf((char*)"(TXMLPGPU) 3-arg node ctor needs layers==3");
    fParm.fInScale = 1.0 / inputRange;
    BuildLayers(3, innodes, {n0, n1, n2}, {s0, s1, s2}, reluOutput);
    InitNet();
}

TXMLPGPU::TXMLPGPU(std::string netFile) : VSupervisedNet(netFile)
{
    ReadNet("XMLPGPU");
}

TXMLPGPU::TXMLPGPU(TXMLP& src, std::string netFile)
    : VSupervisedNet("XMLPGPU", src.GetParameters().fInNodes,
                     src.GetParameters().fOutNodes, netFile)
{
    const TNeuralNetParameters& sp = src.GetParameters();
    fParm.fInScale = sp.fInScale;
    fParm.fLayers = sp.fLayers;
    fParm.fInNodes = sp.fInNodes;
    fParm.fOutNodes = sp.fOutNodes;
    fParm.fTransferId = TNeuralNetParameters::TR_RELU;
    fReLUOutput = false;

    std::vector<int> nodes(sp.fLayers);
    std::vector<double> steps(sp.fLayers, 0.01);
    for (int L = 0; L < sp.fLayers; ++L) {
        TPerceptron* p = src.GetPerceptron(L);
        nodes[L] = p->GetParameters().fOutNodes;
        steps[L] = p->GetParameters().fLearnStep;
    }
    BuildLayers(sp.fLayers, sp.fInNodes, nodes, steps, false);

    for (int L = 0; L < sp.fLayers; ++L) {
        TPerceptron* p = src.GetPerceptron(L);
        LayerF& ly = fLayers[L];
        ly.fReLU = (L < sp.fLayers - 1) ||
                   (p->GetParameters().fTransferId == TNeuralNetParameters::TR_RELU);
        if (L == sp.fLayers - 1) fReLUOutput = ly.fReLU;
        for (int o = 0; o < ly.fOutNodes; ++o) {
            PerceptronUnit* up = &p->fU[o];
            ly.fBFloat[o] = -up->fThreshold;
            for (int i = 0; i < ly.fInNodes; ++i)
                ly.fWFloat[static_cast<size_t>(o) * ly.fInNodes + i] = up->fVector[i];
        }
    }
}

TXMLPGPU::~TXMLPGPU()
{
    if (fFilename != "" && fShouldSave) Save();
}

void TXMLPGPU::AllocNet()
{
    std::vector<int> nodes(fParm.fLayers, 1);
    std::vector<double> steps(fParm.fLayers, 0.01);
    BuildLayers(fParm.fLayers, fParm.fInNodes, nodes, steps, fReLUOutput);
}

void TXMLPGPU::InitNet()
{
    for (auto& ly : fLayers) {
        const double sigma = std::sqrt(2.0 / std::max(1, ly.fInNodes));
        for (size_t k = 0; k < ly.fWFloat.size(); ++k) {
            double u1 = 0.0, u2 = 0.0;
            while (u1 < 1e-12) u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
            u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
            double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            ly.fWFloat[k] = z * sigma;
        }
        std::fill(ly.fBFloat.begin(), ly.fBFloat.end(), 0.0);
    }
    InvalidateGpu();
}

void TXMLPGPU::SetMomentumTerm(double f)
{
    fMomentum = f;
    fParm.fMu = f;
}

// ----- HIP init / upload / forward ------------------------------------------

void TXMLPGPU::InvalidateGpu()
{
    if (fGpu) fGpu->weights_dirty = true;
    fOnGPU = false;
}

bool TXMLPGPU::HasGpuSupport()
{
#ifdef RHONNO_WITH_GPU
    return true;
#else
    return false;
#endif
}

void TXMLPGPU::SetDisableGPU(bool v)
{
    if (fDisableGPU == v) return;
    fDisableGPU = v;
    InvalidateGpu();
}

#ifdef RHONNO_WITH_GPU
static bool gpu_upload(TXMLPGPUGpuState& g, const std::vector<TXMLPGPU::LayerF>& layers,
                       int max_width)
{
    // Free any previous weights
    for (float* p : g.d_W) if (p) hipFree(p);
    for (float* p : g.d_B) if (p) hipFree(p);
    for (float* p : g.d_Y) if (p) hipFree(p);
    g.d_W.clear(); g.d_B.clear(); g.d_Y.clear();

    // Allocate per-layer tensors
    for (const auto& ly : layers) {
        float* w = nullptr;
        float* b = nullptr;
        float* y = nullptr;
        if (hipMalloc(&w, sizeof(float) * ly.fWFloat.size()) != hipSuccess) return false;
        if (hipMalloc(&b, sizeof(float) * ly.fBFloat.size()) != hipSuccess) { hipFree(w); return false; }
        if (hipMalloc(&y, sizeof(float) * ly.fOutNodes) != hipSuccess) { hipFree(w); hipFree(b); return false; }

        std::vector<float> wf(ly.fWFloat.begin(), ly.fWFloat.end());
        std::vector<float> bf(ly.fBFloat.begin(), ly.fBFloat.end());
        if (hipMemcpy(w, wf.data(), wf.size() * sizeof(float), hipMemcpyHostToDevice) != hipSuccess) { hipFree(w); hipFree(b); hipFree(y); return false; }
        if (hipMemcpy(b, bf.data(), bf.size() * sizeof(float), hipMemcpyHostToDevice) != hipSuccess) { hipFree(w); hipFree(b); hipFree(y); return false; }

        g.d_W.push_back(w);
        g.d_B.push_back(b);
        g.d_Y.push_back(y);
    }

    if (g.d_cur_cap < max_width) {
        if (g.d_cur) hipFree(g.d_cur);
        if (hipMalloc(&g.d_cur, sizeof(float) * max_width) != hipSuccess) { g.d_cur = nullptr; g.d_cur_cap = 0; return false; }
        g.d_cur_cap = max_width;
    }

    g.h_x.assign(max_width, 0.0f);
    return true;
}

static bool gpu_forward(TXMLPGPUGpuState& g,
                        const std::vector<TXMLPGPU::LayerF>& layers,
                        const float* x_host, int // x size = layers[0].fInNodes
                        , double* out, int out_n)
{
    // Upload activation to d_cur
    if (!g.d_cur) return false;
    const auto& first = layers[0];
    if (hipMemcpy(g.d_cur, x_host, sizeof(float) * first.fInNodes, hipMemcpyHostToDevice) != hipSuccess)
        return false;

    int block = 256;

    for (size_t L = 0; L < layers.size(); ++L) {
        const auto& ly = layers[L];
        // y = b + W*x  (row-major [out][in]); optional ReLU
        int grid = (ly.fOutNodes + block - 1) / block;
        hipLaunchKernelGGL(k_gemv, dim3(grid), dim3(block), 0, 0,
                           g.d_W[L], g.d_cur, g.d_B[L], g.d_Y[L],
                           ly.fInNodes, ly.fOutNodes, ly.fReLU ? 1 : 0);
        if (hipGetLastError() != hipSuccess) return false;

        // next layer input = this layer's output
        if (L + 1 < layers.size()) {
            if (hipMemcpy(g.d_cur, g.d_Y[L], sizeof(float) * ly.fOutNodes,
                          hipMemcpyDeviceToDevice) != hipSuccess) return false;
        }
    }

    // Download final activation from d_Y[last]
    std::vector<float> h(out_n);
    if (hipMemcpy(h.data(), g.d_Y.back(), sizeof(float) * out_n,
                  hipMemcpyDeviceToHost) != hipSuccess) return false;
    for (int o = 0; o < out_n; ++o) out[o] = static_cast<double>(h[o]);
    return true;
}
#endif // RHONNO_WITH_GPU

// ----- public forward paths -------------------------------------------------

void TXMLPGPU::ForwardFloat(const double* inScaled, double* out) const
{
    std::vector<double> cur(inScaled, inScaled + fParm.fInNodes);
    std::vector<double> nxt;
    for (size_t L = 0; L < fLayers.size(); ++L) {
        const LayerF& ly = fLayers[L];
        nxt.assign(ly.fOutNodes, 0.0);
        for (int o = 0; o < ly.fOutNodes; ++o) {
            double sum = ly.fBFloat[o];
            const double* w = &ly.fWFloat[static_cast<size_t>(o) * ly.fInNodes];
            for (int i = 0; i < ly.fInNodes; ++i) sum += w[i] * cur[i];
            if (ly.fReLU) sum = relu(sum);
            nxt[o] = sum;
        }
        cur.swap(nxt);
    }
    for (int o = 0; o < fParm.fOutNodes; ++o) out[o] = cur[o];
}

double* TXMLPGPU::Recall(NNO_INTYPE* in, NNO_OUTTYPE* out)
{
    // Scale input once on host
    std::vector<float> x(fParm.fInNodes);
    for (int i = 0; i < fParm.fInNodes; ++i)
        x[i] = static_cast<float>(static_cast<double>(in[i]) * fParm.fInScale);

#ifdef RHONNO_WITH_GPU
    if (!fDisableGPU) {
        if (!fGpu) fGpu = std::make_unique<TXMLPGPUGpuState>();
        if (fGpu->weights_dirty) {
            int maxw = fParm.fInNodes;
            for (auto& ly : fLayers) maxw = std::max(maxw, ly.fOutNodes);
            fGpu->gpu_usable = gpu_upload(*fGpu, fLayers, maxw);
            fGpu->weights_dirty = false;
        }
        if (fGpu->gpu_usable &&
            gpu_forward(*fGpu, fLayers, x.data(), fParm.fInNodes,
                        fOut, fParm.fOutNodes)) {
            fOnGPU = true;
            if (fPlotter) {
                bool good = true;
                if (out != 0) good = out[0] > fParm.fThreshold;
                fPlotter->AddTestSample(fOut[0], good);
            }
            return fOut;
        }
        fOnGPU = false;   // fall through to CPU
    }
#endif

    // CPU fallback
    std::vector<double> xd(fParm.fInNodes);
    for (int i = 0; i < fParm.fInNodes; ++i) xd[i] = x[i];
    ForwardFloat(xd.data(), fOut);
    if (fPlotter) {
        bool good = true;
        if (out != 0) good = out[0] > fParm.fThreshold;
        fPlotter->AddTestSample(fOut[0], good);
    }
    return fOut;
}

double TXMLPGPU::Train(NNO_INTYPE* in, NNO_OUTTYPE* trout)
{
    fShouldSave = true;
    InvalidateGpu(); // float weights changed

    std::vector<double> x(fParm.fInNodes);
    for (int i = 0; i < fParm.fInNodes; ++i)
        x[i] = static_cast<double>(in[i]) * fParm.fInScale;

    std::vector<std::vector<double>> acts(fLayers.size() + 1);
    acts[0] = x;
    for (size_t L = 0; L < fLayers.size(); ++L) {
        const LayerF& ly = fLayers[L];
        acts[L + 1].assign(ly.fOutNodes, 0.0);
        for (int o = 0; o < ly.fOutNodes; ++o) {
            double sum = ly.fBFloat[o];
            const double* w = &ly.fWFloat[static_cast<size_t>(o) * ly.fInNodes];
            const double* xi = acts[L].data();
            for (int i = 0; i < ly.fInNodes; ++i) sum += w[i] * xi[i];
            if (ly.fReLU) sum = relu(sum);
            acts[L + 1][o] = sum;
        }
    }
    for (int o = 0; o < fParm.fOutNodes; ++o) fOut[o] = acts.back()[o];

    double S_Err = 0.0;
    std::vector<double> delta(fParm.fOutNodes);
    for (int o = 0; o < fParm.fOutNodes; ++o) {
        double diff = static_cast<double>(trout[o]) - fOut[o];
        S_Err += diff * diff;
        double d = diff;
        if (fLayers.back().fReLU) d *= relu_deriv(fOut[o]);
        delta[o] = d;
    }

    for (int L = static_cast<int>(fLayers.size()) - 1; L >= 0; --L) {
        LayerF& ly = fLayers[L];
        const double lr = ly.fLearnStep;
        std::vector<double> delta_prev(ly.fInNodes, 0.0);
        std::vector<double> wcopy = ly.fWFloat;

        for (int o = 0; o < ly.fOutNodes; ++o) {
            double d = delta[o];
            if (!std::isfinite(d)) d = 0.0;
            if (d > 5.0) d = 5.0;
            if (d < -5.0) d = -5.0;

            double nb = ly.fBFloat[o] + d * lr;
            if (std::isfinite(nb)) ly.fBFloat[o] = nb;

            double* w = &ly.fWFloat[static_cast<size_t>(o) * ly.fInNodes];
            const double* w0 = &wcopy[static_cast<size_t>(o) * ly.fInNodes];
            const double* xi = acts[L].data();
            for (int i = 0; i < ly.fInNodes; ++i) {
                double nw = w[i] + xi[i] * d * lr;
                if (std::isfinite(nw)) {
                    if (nw > 10.0) nw = 10.0;
                    if (nw < -10.0) nw = -10.0;
                    w[i] = nw;
                }
                delta_prev[i] += w0[i] * d;
            }
        }

        if (L > 0) {
            if (fLayers[L - 1].fReLU) {
                for (int i = 0; i < ly.fInNodes; ++i)
                    delta_prev[i] *= relu_deriv(acts[L][i]);
            }
            delta.swap(delta_prev);
        }
    }

    if (fPlotter) fPlotter->AddTrainSample(trout[0], trout[0] > fParm.fThreshold);
    return S_Err;
}

double TXMLPGPU::CompareCpuVsGpu(NNO_INTYPE* in)
{
    std::vector<double> xf(fParm.fInNodes), y_cpu(fParm.fOutNodes),
        y_gpu(fParm.fOutNodes);
    for (int i = 0; i < fParm.fInNodes; ++i)
        xf[i] = static_cast<double>(in[i]) * fParm.fInScale;
    ForwardFloat(xf.data(), y_cpu.data());

    // Force upload+GPU recall
    bool prev = fDisableGPU;
    fDisableGPU = false;
    InvalidateGpu();
    Recall(in, nullptr);
    for (int o = 0; o < fParm.fOutNodes; ++o) y_gpu[o] = fOut[o];
    fDisableGPU = prev;

    double md = 0.0;
    for (int o = 0; o < fParm.fOutNodes; ++o)
        md = std::max(md, std::fabs(y_cpu[o] - y_gpu[o]));
    return md;
}

// ----- persistence ----------------------------------------------------------

void TXMLPGPU::WriteText()
{
    fprintf(fFile, "layers    %i\n", fParm.fLayers);
    fprintf(fFile, "in_scale  %le\n", fParm.fInScale);
    fprintf(fFile, "in_nodes  %i\n", fParm.fInNodes);
    fprintf(fFile, "out_nodes %i\n", fParm.fOutNodes);
    fprintf(fFile, "relu_out  %i\n", fReLUOutput ? 1 : 0);
    fprintf(fFile, "momentum  %le\n", fMomentum);
    for (int L = 0; L < fParm.fLayers; ++L) {
        const LayerF& ly = fLayers[L];
        fprintf(fFile, "\nlayer %i\n", L);
        fprintf(fFile, "innodes     %i\n", ly.fInNodes);
        fprintf(fFile, "outnodes    %i\n", ly.fOutNodes);
        fprintf(fFile, "learn_step  %le\n", ly.fLearnStep);
        fprintf(fFile, "relu        %i\n", ly.fReLU ? 1 : 0);
        fprintf(fFile, "weights_f\n");
        for (double w : ly.fWFloat) fprintf(fFile, "%le\n", w);
        fprintf(fFile, "bias_f\n");
        for (double b : ly.fBFloat) fprintf(fFile, "%le\n", b);
    }
}

void TXMLPGPU::ReadText()
{
    int layers = 0, relu_out = 0;
    double scale = 1.0, mom = 0.0;
    int in_nodes = 0, out_nodes = 0;
    fscanf(fFile, "layers    %i\n", &layers);
    fscanf(fFile, "in_scale  %le\n", &scale);
    fscanf(fFile, "in_nodes  %i\n", &in_nodes);
    fscanf(fFile, "out_nodes %i\n", &out_nodes);
    fscanf(fFile, "relu_out  %i\n", &relu_out);
    fscanf(fFile, "momentum  %le\n", &mom);
    fParm.fLayers = layers;
    fParm.fInScale = scale;
    fParm.fInNodes = in_nodes;
    fParm.fOutNodes = out_nodes;
    fReLUOutput = relu_out != 0;
    fMomentum = mom;

    fLayers.clear();
    fLayers.resize(layers);
    for (int L = 0; L < layers; ++L) {
        int lid = 0, inn = 0, outn = 0, relu = 0;
        double ls = 0.01;
        fscanf(fFile, "\nlayer %i\n", &lid);
        fscanf(fFile, "innodes     %i\n", &inn);
        fscanf(fFile, "outnodes    %i\n", &outn);
        fscanf(fFile, "learn_step  %le\n", &ls);
        fscanf(fFile, "relu        %i\n", &relu);
        LayerF& ly = fLayers[L];
        ly.fInNodes = inn;
        ly.fOutNodes = outn;
        ly.fLearnStep = ls;
        ly.fReLU = relu != 0;
        ly.fWFloat.assign(static_cast<size_t>(outn) * inn, 0.0);
        ly.fBFloat.assign(outn, 0.0);
        fscanf(fFile, "weights_f\n");
        for (size_t k = 0; k < ly.fWFloat.size(); ++k)
            fscanf(fFile, "%le\n", &ly.fWFloat[k]);
        fscanf(fFile, "bias_f\n");
        for (int o = 0; o < outn; ++o)
            fscanf(fFile, "%le\n", &ly.fBFloat[o]);
    }
    if (fOut) delete[] fOut;
    fOut = new double[fParm.fOutNodes];
    InvalidateGpu();
}

void TXMLPGPU::WriteBinary() { WriteText(); }
void TXMLPGPU::ReadBinary()  { ReadText(); }
