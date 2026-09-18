// TXMLPInt8
//
// Int8-quantized Multi-Layer Perceptron with ReLU
// Part of the Neural Network Objects package (NNO)

#include "TXMLPInt8.h"
#include "VNeuralNetPlotter.h"
#include "TPerceptron.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

inline double relu(double x) { return x > 0.0 ? x : 0.0; }
inline double relu_deriv(double y) { return y > 0.0 ? 1.0 : 0.0; }

} // namespace

int8_t TXMLPInt8::QuantizeValue(double v, float scale)
{
    if (!(scale > 0.f) || !std::isfinite(v)) return 0;
    double q = std::round(v / static_cast<double>(scale));
    if (q > 127.0) q = 127.0;
    if (q < -127.0) q = -127.0; // keep -128 free for some kernels; symmetric
    return static_cast<int8_t>(q);
}

float TXMLPInt8::SymmetricScale(const double* data, int n, double eps)
{
    double amax = 0.0;
    for (int i = 0; i < n; ++i) {
        double a = std::fabs(data[i]);
        if (a > amax) amax = a;
    }
    if (amax < eps) amax = eps;
    return static_cast<float>(amax / 127.0);
}

void TXMLPInt8::BuildLayers(int layers, int innodes,
                            const std::vector<int>& nodes,
                            const std::vector<double>& steps,
                            bool reluOutput)
{
    if (layers < 1) Errorf((char*)"(TXMLPInt8) at least one layer necessary");
    if ((int)nodes.size() < layers) Errorf((char*)"(TXMLPInt8) nodes size < layers");
    if ((int)steps.size() < layers) Errorf((char*)"(TXMLPInt8) learnSteps size < layers");

    fParm.fLayers = layers;
    fParm.fInNodes = innodes;
    fParm.fOutNodes = nodes[layers - 1];
    fParm.fTransferId = TNeuralNetParameters::TR_RELU;
    fReLUOutput = reluOutput;
    fQuantized = false;

    fLayers.clear();
    fLayers.resize(layers);

    int prev = innodes;
    for (int L = 0; L < layers; ++L) {
        LayerQ& ly = fLayers[L];
        ly.fInNodes = prev;
        ly.fOutNodes = nodes[L];
        ly.fLearnStep = steps[L];
        ly.fReLU = (L < layers - 1) ? true : reluOutput;
        ly.fWFloat.assign(static_cast<size_t>(ly.fOutNodes) * ly.fInNodes, 0.0);
        ly.fBFloat.assign(ly.fOutNodes, 0.0);
        ly.fW.assign(ly.fWFloat.size(), 0);
        ly.fB.assign(ly.fOutNodes, 0);
        ly.fWScale = 1.f;
        ly.fXScale = 1.f;
        prev = ly.fOutNodes;
    }

    fOut = new double[fParm.fOutNodes];
    TestPointer(fOut);
    fAct.assign(std::max(innodes, fParm.fOutNodes) + 8, 0.0);
    // scratch sized to max layer width
    int maxw = innodes;
    for (int L = 0; L < layers; ++L) maxw = std::max(maxw, fLayers[L].fOutNodes);
    fAct.assign(maxw + 8, 0.0);
    fDelta.assign(maxw + 8, 0.0);
    fActQ.assign(maxw + 8, 0);
}

TXMLPInt8::TXMLPInt8(int layers, double inputRange, std::string netFile,
                     int innodes, const std::vector<int>& nodes,
                     const std::vector<double>& learnSteps, bool reluOutput)
    : VSupervisedNet("XMLP8", innodes, 0, netFile)
{
    fParm.fInScale = 1.0 / inputRange;
    BuildLayers(layers, innodes, nodes, learnSteps, reluOutput);
    InitNet();
}

TXMLPInt8::TXMLPInt8(int layers, double inputRange, std::string netFile,
                     int innodes, int n0, int n1, int n2,
                     double s0, double s1, double s2, bool reluOutput)
    : VSupervisedNet("XMLP8", innodes, 0, netFile)
{
    if (layers != 3) Errorf((char*)"(TXMLPInt8) 3-arg node ctor needs layers==3");
    fParm.fInScale = 1.0 / inputRange;
    BuildLayers(3, innodes, {n0, n1, n2}, {s0, s1, s2}, reluOutput);
    InitNet();
}

TXMLPInt8::TXMLPInt8(std::string netFile) : VSupervisedNet(netFile)
{
    ReadNet("XMLP8");
}

TXMLPInt8::TXMLPInt8(TXMLP& src, std::string netFile)
    : VSupervisedNet("XMLP8", src.GetParameters().fInNodes, src.GetParameters().fOutNodes, netFile)
{
    const TNeuralNetParameters& sp = src.GetParameters();
    fParm.fInScale = sp.fInScale;
    fParm.fLayers = sp.fLayers;
    fParm.fInNodes = sp.fInNodes;
    fParm.fOutNodes = sp.fOutNodes;
    fParm.fTransferId = TNeuralNetParameters::TR_RELU;
    fReLUOutput = false;
    fQuantized = false;

    std::vector<int> nodes(sp.fLayers);
    std::vector<double> steps(sp.fLayers, 0.01);
    for (int L = 0; L < sp.fLayers; ++L) {
        TPerceptron* p = src.GetPerceptron(L);
        nodes[L] = p->GetParameters().fOutNodes;
        steps[L] = p->GetParameters().fLearnStep;
    }
    BuildLayers(sp.fLayers, sp.fInNodes, nodes, steps, /*reluOutput*/ false);

    // Copy float weights from TPerceptron units
    for (int L = 0; L < sp.fLayers; ++L) {
        TPerceptron* p = src.GetPerceptron(L);
        LayerQ& ly = fLayers[L];
        ly.fReLU = (L < sp.fLayers - 1) ||
                   (p->GetParameters().fTransferId == TNeuralNetParameters::TR_RELU);
        if (L == sp.fLayers - 1) fReLUOutput = ly.fReLU;

        for (int o = 0; o < ly.fOutNodes; ++o) {
            PerceptronUnit* up = &p->fU[o];
            ly.fBFloat[o] = -up->fThreshold; // TPerceptron uses sum - threshold
            for (int i = 0; i < ly.fInNodes; ++i)
                ly.fWFloat[static_cast<size_t>(o) * ly.fInNodes + i] = up->fVector[i];
        }
    }
    Quantize();
}

TXMLPInt8::~TXMLPInt8()
{
    if (fFilename != "" && fShouldSave) Save();
}

void TXMLPInt8::AllocNet()
{
    // parameters already in fParm from file header
    std::vector<int> nodes(fParm.fLayers, 0);
    std::vector<double> steps(fParm.fLayers, 0.01);
    // Real sizes filled in Read*; placeholder Build then overwritten
    for (int i = 0; i < fParm.fLayers; ++i) nodes[i] = 1;
    BuildLayers(fParm.fLayers, fParm.fInNodes, nodes, steps, fReLUOutput);
}

void TXMLPInt8::InitNet()
{
    // He initialization for ReLU: N(0, sqrt(2/fan_in))
    for (auto& ly : fLayers) {
        const double sigma = std::sqrt(2.0 / std::max(1, ly.fInNodes));
        for (size_t k = 0; k < ly.fWFloat.size(); ++k) {
            double u1 = 0.0, u2 = 0.0;
            // open interval (0,1]
            while (u1 < 1e-12) u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
            u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
            double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            ly.fWFloat[k] = z * sigma;
        }
        for (int o = 0; o < ly.fOutNodes; ++o)
            ly.fBFloat[o] = 0.0;
    }
    fQuantized = false;
}

void TXMLPInt8::SetMomentumTerm(double f)
{
    fMomentum = f;
    fParm.fMu = f;
}

void TXMLPInt8::Quantize()
{
    // Input scale from network input range (after fInScale, |x| ~ 1 typically)
    float xscale = static_cast<float>(1.0 / 127.0); // default if all zero
    // Use |scaled input| max ~1.0 when fInScale applied to data in [-inputRange,inputRange]
    xscale = 1.f / 127.f;

    for (size_t L = 0; L < fLayers.size(); ++L) {
        LayerQ& ly = fLayers[L];
        ly.fXScale = xscale;
        ly.fWScale = SymmetricScale(ly.fWFloat.data(), static_cast<int>(ly.fWFloat.size()));

        ly.fW.resize(ly.fWFloat.size());
        for (size_t k = 0; k < ly.fWFloat.size(); ++k)
            ly.fW[k] = QuantizeValue(ly.fWFloat[k], ly.fWScale);

        // Bias in accumulator domain: b_real ≈ b_q * (w_scale * x_scale)
        const double b_scale = static_cast<double>(ly.fWScale) * ly.fXScale;
        ly.fB.resize(ly.fOutNodes);
        for (int o = 0; o < ly.fOutNodes; ++o) {
            double bq = (b_scale > 0.0) ? std::round(ly.fBFloat[o] / b_scale) : 0.0;
            if (bq > static_cast<double>(std::numeric_limits<int32_t>::max()))
                bq = static_cast<double>(std::numeric_limits<int32_t>::max());
            if (bq < static_cast<double>(std::numeric_limits<int32_t>::min()))
                bq = static_cast<double>(std::numeric_limits<int32_t>::min());
            ly.fB[o] = static_cast<int32_t>(bq);
        }

        // Estimate next activation scale from a synthetic pass is hard offline;
        // use weight*input range bound: |y| <= in*127*127*w_scale*x_scale roughly
        // Practical: scan float forward weights max activation magnitude via |W|*|x|_max
        double amax = 0.0;
        for (int o = 0; o < ly.fOutNodes; ++o) {
            double s = std::fabs(ly.fBFloat[o]);
            for (int i = 0; i < ly.fInNodes; ++i)
                s += std::fabs(ly.fWFloat[static_cast<size_t>(o) * ly.fInNodes + i]) *
                     (127.0 * ly.fXScale);
            if (s > amax) amax = s;
        }
        if (amax < 1e-8) amax = 1e-8;
        xscale = static_cast<float>(amax / 127.0);
        // Keep xscale for next layer
        if (L + 1 < fLayers.size()) {
            // will set next fXScale at loop head via xscale
        }
    }
    // Fix next-layer x scales properly in a second pass using amax chain
    xscale = 1.f / 127.f;
    for (size_t L = 0; L < fLayers.size(); ++L) {
        fLayers[L].fXScale = xscale;
        double amax = 0.0;
        LayerQ& ly = fLayers[L];
        for (int o = 0; o < ly.fOutNodes; ++o) {
            double s = std::fabs(ly.fBFloat[o]);
            for (int i = 0; i < ly.fInNodes; ++i)
                s += std::fabs(ly.fWFloat[static_cast<size_t>(o) * ly.fInNodes + i]) *
                     (127.0 * ly.fXScale);
            if (ly.fReLU && s < 0) s = 0; // not applicable to bound
            if (s > amax) amax = s;
        }
        if (amax < 1e-8) amax = 1e-8;
        xscale = static_cast<float>(amax / 127.0);

        // Recompute bias quant with final x scale
        const double b_scale = static_cast<double>(ly.fWScale) * ly.fXScale;
        for (int o = 0; o < ly.fOutNodes; ++o) {
            double bq = (b_scale > 0.0) ? std::round(ly.fBFloat[o] / b_scale) : 0.0;
            ly.fB[o] = static_cast<int32_t>(std::max(
                static_cast<double>(std::numeric_limits<int32_t>::min()),
                std::min(static_cast<double>(std::numeric_limits<int32_t>::max()), bq)));
        }
    }

    fQuantized = true;
}

void TXMLPInt8::ForwardFloat(const double* inScaled, double* out)
{
    // fAct holds layer input
    std::vector<double> cur(inScaled, inScaled + fParm.fInNodes);
    std::vector<double> nxt;

    for (size_t L = 0; L < fLayers.size(); ++L) {
        const LayerQ& ly = fLayers[L];
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

void TXMLPInt8::ForwardInt8(const NNO_INTYPE* in, double* out)
{
    if (!fQuantized) Quantize();

    // Quantize input
    std::vector<int8_t> cur(fParm.fInNodes);
    const float in_scale = fLayers.empty() ? (1.f / 127.f) : fLayers[0].fXScale;
    for (int i = 0; i < fParm.fInNodes; ++i) {
        double xs = static_cast<double>(in[i]) * fParm.fInScale;
        cur[i] = QuantizeValue(xs, in_scale);
    }

    std::vector<int8_t> nxt;
    std::vector<double> yf;

    for (size_t L = 0; L < fLayers.size(); ++L) {
        const LayerQ& ly = fLayers[L];
        // Ensure input uses this layer's x scale (requantize if needed)
        if (L == 0) {
            // already using fLayers[0].fXScale
        }
        yf.resize(ly.fOutNodes);
        const double deq = static_cast<double>(ly.fWScale) * ly.fXScale;
        for (int o = 0; o < ly.fOutNodes; ++o) {
            int32_t acc = ly.fB[o];
            const int8_t* w = &ly.fW[static_cast<size_t>(o) * ly.fInNodes];
            for (int i = 0; i < ly.fInNodes; ++i)
                acc += static_cast<int32_t>(w[i]) * static_cast<int32_t>(cur[i]);
            double y = static_cast<double>(acc) * deq;
            if (ly.fReLU) y = relu(y);
            yf[o] = y;
        }

        if (L + 1 < fLayers.size()) {
            const float next_xs = fLayers[L + 1].fXScale;
            nxt.resize(ly.fOutNodes);
            for (int o = 0; o < ly.fOutNodes; ++o)
                nxt[o] = QuantizeValue(yf[o], next_xs);
            cur.swap(nxt);
        } else {
            for (int o = 0; o < ly.fOutNodes; ++o) out[o] = yf[o];
        }
    }
}

double* TXMLPInt8::Recall(NNO_INTYPE* in, NNO_OUTTYPE* out)
{
    ForwardInt8(in, fOut);

    if (fPlotter) {
        bool good = true;
        if (out != 0) good = out[0] > fParm.fThreshold;
        fPlotter->AddTestSample(fOut[0], good);
    }
    return fOut;
}

double TXMLPInt8::Train(NNO_INTYPE* in, NNO_OUTTYPE* trout)
{
    fShouldSave = true;
    fQuantized = false; // float weights changed

    // Scale input
    std::vector<double> x(fParm.fInNodes);
    for (int i = 0; i < fParm.fInNodes; ++i)
        x[i] = static_cast<double>(in[i]) * fParm.fInScale;

    // Forward float, cache post-activation per layer (acts[0]=input)
    std::vector<std::vector<double>> acts(fLayers.size() + 1);
    acts[0] = x;
    for (size_t L = 0; L < fLayers.size(); ++L) {
        const LayerQ& ly = fLayers[L];
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

    // Squared error + output delta (matches TXMLP: delta = target - output)
    double S_Err = 0.0;
    std::vector<double> delta(fParm.fOutNodes);
    for (int o = 0; o < fParm.fOutNodes; ++o) {
        double diff = static_cast<double>(trout[o]) - fOut[o];
        S_Err += diff * diff;
        double d = diff;
        if (fLayers.back().fReLU) d *= relu_deriv(fOut[o]);
        delta[o] = d;
    }

    // Backprop: snapshot weights, update, propagate with pre-update W
    for (int L = static_cast<int>(fLayers.size()) - 1; L >= 0; --L) {
        LayerQ& ly = fLayers[L];
        const double lr = ly.fLearnStep;
        std::vector<double> delta_prev(ly.fInNodes, 0.0);
        std::vector<double> wcopy = ly.fWFloat;

        for (int o = 0; o < ly.fOutNodes; ++o) {
            double d = delta[o];
            if (!std::isfinite(d)) d = 0.0;
            // clip delta
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
                    // soft weight clip
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

double TXMLPInt8::CompareFloatVsInt8(NNO_INTYPE* in)
{
    std::vector<double> xf(fParm.fInNodes), yf(fParm.fOutNodes), yq(fParm.fOutNodes);
    for (int i = 0; i < fParm.fInNodes; ++i)
        xf[i] = static_cast<double>(in[i]) * fParm.fInScale;
    ForwardFloat(xf.data(), yf.data());
    bool was = fQuantized;
    Quantize();
    ForwardInt8(in, yq.data());
    fQuantized = was;
    double maxd = 0.0;
    for (int o = 0; o < fParm.fOutNodes; ++o)
        maxd = std::max(maxd, std::fabs(yf[o] - yq[o]));
    return maxd;
}

void TXMLPInt8::WriteText()
{
    fprintf(fFile, "layers    %i\n", fParm.fLayers);
    fprintf(fFile, "in_scale  %le\n", fParm.fInScale);
    fprintf(fFile, "in_nodes  %i\n", fParm.fInNodes);
    fprintf(fFile, "out_nodes %i\n", fParm.fOutNodes);
    fprintf(fFile, "relu_out  %i\n", fReLUOutput ? 1 : 0);
    fprintf(fFile, "momentum  %le\n", fMomentum);
    if (!fQuantized) Quantize();
    for (int L = 0; L < fParm.fLayers; ++L) {
        const LayerQ& ly = fLayers[L];
        fprintf(fFile, "\nlayer %i\n", L);
        fprintf(fFile, "innodes     %i\n", ly.fInNodes);
        fprintf(fFile, "outnodes    %i\n", ly.fOutNodes);
        fprintf(fFile, "learn_step  %le\n", ly.fLearnStep);
        fprintf(fFile, "relu        %i\n", ly.fReLU ? 1 : 0);
        fprintf(fFile, "w_scale     %a\n", ly.fWScale);
        fprintf(fFile, "x_scale     %a\n", ly.fXScale);
        fprintf(fFile, "weights_q\n");
        for (int o = 0; o < ly.fOutNodes; ++o) {
            for (int i = 0; i < ly.fInNodes; ++i)
                fprintf(fFile, "%d\n", static_cast<int>(ly.fW[static_cast<size_t>(o) * ly.fInNodes + i]));
        }
        fprintf(fFile, "bias_q\n");
        for (int o = 0; o < ly.fOutNodes; ++o) fprintf(fFile, "%d\n", ly.fB[o]);
        fprintf(fFile, "weights_f\n");
        for (double w : ly.fWFloat) fprintf(fFile, "%le\n", w);
        fprintf(fFile, "bias_f\n");
        for (double b : ly.fBFloat) fprintf(fFile, "%le\n", b);
    }
}

void TXMLPInt8::ReadText()
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

    std::vector<int> nodes(layers);
    std::vector<double> steps(layers, 0.01);
    // Peek layer headers by reading into temp structure twice — single pass build
    fLayers.clear();
    fLayers.resize(layers);
    int prev = in_nodes;
    for (int L = 0; L < layers; ++L) {
        int lid = 0, inn = 0, outn = 0, relu = 0;
        double ls = 0.01;
        float ws = 1.f, xs = 1.f;
        fscanf(fFile, "\nlayer %i\n", &lid);
        fscanf(fFile, "innodes     %i\n", &inn);
        fscanf(fFile, "outnodes    %i\n", &outn);
        fscanf(fFile, "learn_step  %le\n", &ls);
        fscanf(fFile, "relu        %i\n", &relu);
        fscanf(fFile, "w_scale     %a\n", &ws);
        fscanf(fFile, "x_scale     %a\n", &xs);
        LayerQ& ly = fLayers[L];
        ly.fInNodes = inn;
        ly.fOutNodes = outn;
        ly.fLearnStep = ls;
        ly.fReLU = relu != 0;
        ly.fWScale = ws;
        ly.fXScale = xs;
        ly.fW.assign(static_cast<size_t>(outn) * inn, 0);
        ly.fB.assign(outn, 0);
        ly.fWFloat.assign(ly.fW.size(), 0.0);
        ly.fBFloat.assign(outn, 0.0);
        fscanf(fFile, "weights_q\n");
        for (int o = 0; o < outn; ++o)
            for (int i = 0; i < inn; ++i) {
                int v = 0;
                fscanf(fFile, "%d\n", &v);
                ly.fW[static_cast<size_t>(o) * inn + i] = static_cast<int8_t>(v);
            }
        fscanf(fFile, "bias_q\n");
        for (int o = 0; o < outn; ++o) {
            int v = 0;
            fscanf(fFile, "%d\n", &v);
            ly.fB[o] = v;
        }
        fscanf(fFile, "weights_f\n");
        for (size_t k = 0; k < ly.fWFloat.size(); ++k) fscanf(fFile, "%le\n", &ly.fWFloat[k]);
        fscanf(fFile, "bias_f\n");
        for (int o = 0; o < outn; ++o) fscanf(fFile, "%le\n", &ly.fBFloat[o]);
        nodes[L] = outn;
        steps[L] = ls;
        prev = outn;
        (void)prev;
    }
    if (fOut) delete[] fOut;
    fOut = new double[fParm.fOutNodes];
    fQuantized = true;
    int maxw = fParm.fInNodes;
    for (auto& ly : fLayers) maxw = std::max(maxw, ly.fOutNodes);
    fAct.assign(maxw + 8, 0.0);
    fDelta.assign(maxw + 8, 0.0);
    fActQ.assign(maxw + 8, 0);
}

void TXMLPInt8::WriteBinary()
{
    if (!fQuantized) Quantize();
    fwrite(&fParm, sizeof(TNeuralNetParameters), 1, fFile);
    int ro = fReLUOutput ? 1 : 0;
    fwrite(&ro, sizeof(int), 1, fFile);
    fwrite(&fMomentum, sizeof(double), 1, fFile);
    for (int L = 0; L < fParm.fLayers; ++L) {
        LayerQ& ly = fLayers[L];
        fwrite(&ly.fInNodes, sizeof(int), 1, fFile);
        fwrite(&ly.fOutNodes, sizeof(int), 1, fFile);
        fwrite(&ly.fLearnStep, sizeof(double), 1, fFile);
        int relu = ly.fReLU ? 1 : 0;
        fwrite(&relu, sizeof(int), 1, fFile);
        fwrite(&ly.fWScale, sizeof(float), 1, fFile);
        fwrite(&ly.fXScale, sizeof(float), 1, fFile);
        fwrite(ly.fW.data(), sizeof(int8_t), ly.fW.size(), fFile);
        fwrite(ly.fB.data(), sizeof(int32_t), ly.fB.size(), fFile);
        fwrite(ly.fWFloat.data(), sizeof(double), ly.fWFloat.size(), fFile);
        fwrite(ly.fBFloat.data(), sizeof(double), ly.fBFloat.size(), fFile);
    }
}

void TXMLPInt8::ReadBinary()
{
    fread(&fParm, sizeof(TNeuralNetParameters), 1, fFile);
    int ro = 0;
    fread(&ro, sizeof(int), 1, fFile);
    fReLUOutput = ro != 0;
    fread(&fMomentum, sizeof(double), 1, fFile);
    fLayers.resize(fParm.fLayers);
    for (int L = 0; L < fParm.fLayers; ++L) {
        LayerQ& ly = fLayers[L];
        fread(&ly.fInNodes, sizeof(int), 1, fFile);
        fread(&ly.fOutNodes, sizeof(int), 1, fFile);
        fread(&ly.fLearnStep, sizeof(double), 1, fFile);
        int relu = 0;
        fread(&relu, sizeof(int), 1, fFile);
        ly.fReLU = relu != 0;
        fread(&ly.fWScale, sizeof(float), 1, fFile);
        fread(&ly.fXScale, sizeof(float), 1, fFile);
        ly.fW.resize(static_cast<size_t>(ly.fOutNodes) * ly.fInNodes);
        ly.fB.resize(ly.fOutNodes);
        ly.fWFloat.resize(ly.fW.size());
        ly.fBFloat.resize(ly.fOutNodes);
        fread(ly.fW.data(), sizeof(int8_t), ly.fW.size(), fFile);
        fread(ly.fB.data(), sizeof(int32_t), ly.fB.size(), fFile);
        fread(ly.fWFloat.data(), sizeof(double), ly.fWFloat.size(), fFile);
        fread(ly.fBFloat.data(), sizeof(double), ly.fBFloat.size(), fFile);
    }
    if (fOut) delete[] fOut;
    fOut = new double[fParm.fOutNodes];
    fQuantized = true;
    int maxw = fParm.fInNodes;
    for (auto& ly : fLayers) maxw = std::max(maxw, ly.fOutNodes);
    fAct.assign(maxw + 8, 0.0);
    fDelta.assign(maxw + 8, 0.0);
    fActQ.assign(maxw + 8, 0);
}
