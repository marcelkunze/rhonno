// TXMLPNPU
//
// Float-shadow-trained MLP whose Recall() runs on AMD Ryzen XDNA NPU via
// ONNX Runtime + VitisAI EP. Training happens in float on CPU and is
// numerically equal to TXMLPInt8's shadow path.
//
// Part of the Neural Network Objects package (NNO)

#include "TXMLPNPU.h"
#include "VNeuralNetPlotter.h"
#include "TPerceptron.h"

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

struct TXMLPNPUOrtState; // defined in TXMLPNPU.cxx (Pimpl)

// Define it here (after ORT include) so the header's unique_ptr sees a complete dtor.
struct TXMLPNPUOrtState {
    std::unique_ptr<Ort::Env>         env;
    std::unique_ptr<Ort::Session>     session;
    std::unique_ptr<Ort::MemoryInfo>  mem_info;
};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

inline double relu(double x) { return x > 0.0 ? x : 0.0; }
inline double relu_deriv(double y) { return y > 0.0 ? 1.0 : 0.0; }

// --- Minimal protobuf writer for onnx.ModelProto ---------------------------
// onnx.proto3 stable field numbers needed for a Gemm/Relu MLP graph:
//   ModelProto: ir_version=1 (int64), producer_name=2, graph=7, opset_import=8
//   GraphProto: node=1 (repeated), name=2, initializer=5, input=11, output=12
//   NodeProto:  input=1 (repeated string), output=2 (repeated), name=3, op_type=4
//   TensorProto: dims=1 (repeated int64), data_type=2 (int32), float_data=4 (packed fixed32),
//                name=8, raw_data=9
//   ValueInfoProto: name=1, type=2
//   TypeProto: tensor_type=1 { elem_type=1, shape=2 { dim=1 { dim_value=1 } } }
//   OperatorSetIdProto: domain=1, version=2
// We keep floats in float_data (tag 4, packed) — simpler than raw_data.
namespace pb {

inline void varint_len(std::string& s, size_t n);

inline void tag(std::string& s, int field, int wire) {
    uint64_t v = (static_cast<uint64_t>(field) << 3) | static_cast<uint64_t>(wire);
    while (true) {
        unsigned char b = v & 0x7f;
        v >>= 7;
        if (v) s.push_back(static_cast<char>(b | 0x80));
        else { s.push_back(static_cast<char>(b)); break; }
    }
}
inline void varint(std::string& s, int field, uint64_t v) {
    tag(s, field, 0);
    while (true) {
        unsigned char b = v & 0x7f;
        v >>= 7;
        if (v) s.push_back(static_cast<char>(b | 0x80));
        else { s.push_back(static_cast<char>(b)); break; }
    }
}
inline void bytes(std::string& s, int field, const std::string& b) {
    tag(s, field, 2);
    varint_len(s, b.size());
    s += b;
}
inline void varint_len(std::string& s, size_t n) {
    uint64_t v = n;
    while (true) {
        unsigned char b = v & 0x7f;
        v >>= 7;
        if (v) s.push_back(static_cast<char>(b | 0x80));
        else { s.push_back(static_cast<char>(b)); break; }
    }
}
inline void str_field(std::string& s, int field, const char* v) { bytes(s, field, std::string(v)); }
inline void packed_f32(std::string& s, int field, const float* d, size_t n) {
    tag(s, field, 2);
    varint_len(s, n * 4);
    size_t off = s.size();
    s.resize(off + n * 4);
    std::memcpy(&s[off], d, n * 4);
}

std::string make_dim_value(int64_t v) {
    std::string d;
    varint(d, 1, static_cast<uint64_t>(v)); // dim_value
    return d;
}
std::string make_tensor_shape(int64_t d0, int64_t d1) {
    std::string sh;
    bytes(sh, 1, make_dim_value(d0));       // dim (repeated)
    if (d1 >= 0) bytes(sh, 1, make_dim_value(d1));
    return sh;
}
std::string make_tensor_type(int elem, int64_t d0, int64_t d1) {
    std::string tt;
    varint(tt, 1, static_cast<uint64_t>(elem)); // elem_type (1 = FLOAT)
    bytes(tt, 2, make_tensor_shape(d0, d1));    // shape
    return tt;
}
std::string make_value_info(const std::string& name, int elem, int64_t d0, int64_t d1) {
    std::string vi;
    bytes(vi, 1, name);      // name
    std::string ty;
    bytes(ty, 1, make_tensor_type(elem, d0, d1)); // tensor_type
    bytes(vi, 2, ty);        // type
    return vi;
}
std::string make_tensor_f32(const std::string& name, int64_t rows, int64_t cols,
                            const std::vector<float>& data) {
    std::string t;
    varint(t, 1, static_cast<uint64_t>(rows));                 // dims
    if (cols >= 0) varint(t, 1, static_cast<uint64_t>(cols));  // dims
    varint(t, 2, 1);                                           // data_type = FLOAT
    packed_f32(t, 4, data.data(), data.size());                // float_data
    bytes(t, 8, name);                                         // name
    return t;
}
std::string make_node(const char* op, const std::vector<std::string>& in,
                      const std::vector<std::string>& out, const std::string& name) {
    std::string n;
    for (const auto& i : in) bytes(n, 1, i);
    for (const auto& o : out) bytes(n, 2, o);
    bytes(n, 3, name);
    bytes(n, 4, std::string(op));
    return n;
}
// Gemm attributes: alpha (f=1), beta (f=2), transB (i=4)
std::string make_attr_float(const char* name, float v) {
    std::string a;
    bytes(a, 1, std::string(name));   // name
    tag(a, 2, 5);                     // f (fixed32, wire 5)
    size_t off = a.size();
    a.resize(off + 4);
    std::memcpy(&a[off], &v, 4);
    varint(a, 20, 1);                 // type = FLOAT (AttributeProto.AttributeType)
    return a;
}
std::string make_attr_int(const char* name, int64_t v) {
    std::string a;
    bytes(a, 1, std::string(name));
    varint(a, 3, static_cast<uint64_t>(v)); // i
    varint(a, 20, 2);                        // type = INT
    return a;
}
std::string make_gemm(const std::string& x, const std::string& w, const std::string& b,
                      const std::string& y, const std::string& name) {
    std::string n;
    bytes(n, 1, x); bytes(n, 1, w); bytes(n, 1, b);       // inputs
    bytes(n, 2, y);                                        // output
    bytes(n, 3, name);                                     // name
    bytes(n, 4, std::string("Gemm"));                      // op_type
    bytes(n, 5, make_attr_float("alpha", 1.0f));           // attribute
    bytes(n, 5, make_attr_float("beta", 1.0f));
    bytes(n, 5, make_attr_int("transB", 1));               // W is [out][in] → transB=1
    return n;
}
std::string make_relu(const std::string& in, const std::string& out, const std::string& name) {
    std::string n;
    bytes(n, 1, in);
    bytes(n, 2, out);
    bytes(n, 3, name);
    bytes(n, 4, std::string("Relu"));
    return n;
}
std::string make_opset(int64_t version) {
    std::string o;
    varint(o, 2, static_cast<uint64_t>(version)); // version (domain "" omitted)
    return o;
}

} // namespace pb
} // anonymous namespace

// Query the vaIP compile cache for a VAIML partition summary and report whether
// any subgraph was actually placed on the AIE/NPU. Tiny models often end up
// fully on CPU EP even when VitisAI accepts the graph.
static bool DetectVaimlNpuOffload(const std::string& cacheDir,
                                  const std::string& cacheKey,
                                  std::string* detail)
{
    std::string base = cacheDir;
    if (!base.empty() && base.back() != '/') base += '/';
    base += cacheKey;
    const std::string rai = base + "/" + cacheKey + ".rai";

    // The .rai JSON manifest lists device per subgraph. We count VAIML/AIE entries.
    int vaimlCount = 0;
    int totalDev = 0;
    if (FILE* f = fopen(rai.c_str(), "rb")) {
        fseek(f, 0, SEEK_END);
        long sz = ftell(f);
        fseek(f, 0, SEEK_SET);
        if (sz > 0 && sz < (32 << 20)) {
            std::string buf(sz, '\0');
            fread(&buf[0], 1, sz, f);
            fclose(f);
            size_t pos = 0;
            while ((pos = buf.find("\"device\"", pos)) != std::string::npos) {
                size_t colon = buf.find(':', pos);
                size_t q1 = buf.find('"', colon + 1);
                size_t q2 = q1 == std::string::npos ? q1 : buf.find('"', q1 + 1);
                if (q1 != std::string::npos && q2 != std::string::npos) {
                    std::string dev = buf.substr(q1 + 1, q2 - q1 - 1);
                    totalDev++;
                    if (dev.find("VAIML") != std::string::npos ||
                        dev.find("AIE") != std::string::npos) vaimlCount++;
                }
                pos = q2 != std::string::npos ? q2 + 1 : pos + 7;
            }
        } else if (f) fclose(f);
    }

    std::ostringstream os;
    os << "VAIML_devices=" << vaimlCount << "/" << totalDev << " (rai manifest)";
    if (detail) *detail = os.str();
    return vaimlCount > 0;
}

// ---------------------------------------------------------------------------

void TXMLPNPU::BuildLayers(int layers, int innodes,
                           const std::vector<int>& nodes,
                           const std::vector<double>& steps,
                           bool reluOutput)
{
    if (layers < 1) Errorf((char*)"(TXMLPNPU) at least one layer necessary");
    if ((int)nodes.size() < layers) Errorf((char*)"(TXMLPNPU) nodes size < layers");
    if ((int)steps.size() < layers) Errorf((char*)"(TXMLPNPU) learnSteps size < layers");

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

TXMLPNPU::TXMLPNPU(int layers, double inputRange, std::string netFile,
                   int innodes, const std::vector<int>& nodes,
                   const std::vector<double>& learnSteps, bool reluOutput)
    : VSupervisedNet("XMLPNPU", innodes, 0, netFile)
{
    fParm.fInScale = 1.0 / inputRange;
    BuildLayers(layers, innodes, nodes, learnSteps, reluOutput);
    InitNet();
}

TXMLPNPU::TXMLPNPU(int layers, double inputRange, std::string netFile,
                   int innodes, int n0, int n1, int n2,
                   double s0, double s1, double s2, bool reluOutput)
    : VSupervisedNet("XMLPNPU", innodes, 0, netFile)
{
    if (layers != 3) Errorf((char*)"(TXMLPNPU) 3-arg node ctor needs layers==3");
    fParm.fInScale = 1.0 / inputRange;
    BuildLayers(3, innodes, {n0, n1, n2}, {s0, s1, s2}, reluOutput);
    InitNet();
}

TXMLPNPU::TXMLPNPU(std::string netFile) : VSupervisedNet(netFile)
{
    ReadNet("XMLPNPU");
}

TXMLPNPU::TXMLPNPU(TXMLP& src, std::string netFile)
    : VSupervisedNet("XMLPNPU", src.GetParameters().fInNodes, src.GetParameters().fOutNodes, netFile)
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
    BuildLayers(sp.fLayers, sp.fInNodes, nodes, steps, /*reluOutput*/ false);

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

TXMLPNPU::TXMLPNPU()
{
    // out-of-line so unique_ptr<TXMLPNPUOrtState> sees a complete type
}

TXMLPNPU::~TXMLPNPU()
{
    if (fFilename != "" && fShouldSave) Save();
}

void TXMLPNPU::AllocNet()
{
    std::vector<int> nodes(fParm.fLayers, 1);
    std::vector<double> steps(fParm.fLayers, 0.01);
    BuildLayers(fParm.fLayers, fParm.fInNodes, nodes, steps, fReLUOutput);
}

void TXMLPNPU::InitNet()
{
    // He init for ReLU
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
    InvalidateSession();
}

void TXMLPNPU::SetMomentumTerm(double f)
{
    fMomentum = f;
    fParm.fMu = f;
}

void TXMLPNPU::SetDisableNPU(bool v)
{
    if (fDisableNPU == v) return;
    fDisableNPU = v;
    InvalidateSession();
}

bool TXMLPNPU::DumpOnnx(const char* path) const
{
    std::string m = BuildOnnxModel();
    FILE* f = fopen(path, "wb");
    if (!f) return false;
    fwrite(m.data(), 1, m.size(), f);
    fclose(f);
    return true;
}

void TXMLPNPU::InvalidateSession()
{
    if (fOrt) fOrt->session.reset();
    fSessionValid = false;
    fOnNPU = false;
}

std::string TXMLPNPU::BuildOnnxModel() const
{
    // Graph: input "X" [1, in]
    // For each layer L: initializer W_L [out,in], B_L [out]
    //   y = Gemm(x, W, B) transB=1 → [1, out]
    //   if ReLU: a = Relu(y)
    // Output "Y" [1, out_nodes]
    std::string graph;
    const std::string X = "X";
    std::string cur = X;

    for (int L = 0; L < fParm.fLayers; ++L) {
        const LayerF& ly = fLayers[L];
        std::string wname = "W" + std::to_string(L);
        std::string bname = "B" + std::to_string(L);
        std::string yname = "Y" + std::to_string(L);
        std::string aname = (L + 1 < fParm.fLayers || ly.fReLU) ? ("A" + std::to_string(L)) : "";

        // Initializers: convert double → float32
        std::vector<float> wf(ly.fWFloat.begin(), ly.fWFloat.end());
        std::vector<float> bf(ly.fBFloat.begin(), ly.fBFloat.end());
        pb::bytes(graph, 5, pb::make_tensor_f32(wname, ly.fOutNodes, ly.fInNodes, wf));
        pb::bytes(graph, 5, pb::make_tensor_f32(bname, ly.fOutNodes, -1, bf));

        std::string outname = ly.fReLU ? yname : ((L + 1 == fParm.fLayers) ? "Y" : yname);
        pb::bytes(graph, 1, pb::make_gemm(cur, wname, bname, outname,
                                          "Gemm" + std::to_string(L)));
        if (ly.fReLU) {
            std::string relu_out = (L + 1 == fParm.fLayers) ? "Y" : aname;
            pb::bytes(graph, 1, pb::make_relu(yname, relu_out, "Relu" + std::to_string(L)));
            cur = relu_out;
        } else {
            cur = yname;
        }
    }
    // Make sure final tensor is named "Y"
    if (cur != "Y") {
        // Identity rename in a robust way: emit a final Relu-less layer is wrong,
        // instead we made the last Gemm output "Y" by choice above. When last layer has
        // ReLU true and L==layers-1 we already output "Y" in the relu step.
        // So nothing to do here.
    }

    // Defensive shape metadata for intermediates — the VitisAI/VAIML backend of
    // ORT 1.27 mis-infers Gemm output shapes in a multi-Gemm chain when the
    // intermediate ValueInfo entries are absent. Providing them explicitly is
    // harmless for CPU EP and required for VAIML partitioning.
    for (int L = 0; L < fParm.fLayers - 1; ++L) {
        pb::bytes(graph, 13, pb::make_value_info("Y" + std::to_string(L), 1, 1,
                                                  fLayers[L].fOutNodes)); // value_info
        pb::bytes(graph, 13, pb::make_value_info("A" + std::to_string(L), 1, 1,
                                                  fLayers[L].fOutNodes));
    }

    pb::str_field(graph, 2, "TXMLPNPU"); // name
    pb::bytes(graph, 11, pb::make_value_info("X", 1, 1, fParm.fInNodes));   // input
    pb::bytes(graph, 12, pb::make_value_info("Y", 1, 1, fParm.fOutNodes));  // output

    std::string model;
    pb::varint(model, 1, 8);                    // ir_version = 8 (ORT 1.27 supports ir<=10)
    pb::str_field(model, 2, "rhonno");          // producer_name
    pb::bytes(model, 7, graph);                 // graph
    pb::bytes(model, 8, pb::make_opset(17));    // opset_import default 17
    return model;
}

void TXMLPNPU::EnsureSession()
{
    if (fSessionValid) return;

    if (!fOrt) fOrt = std::make_unique<TXMLPNPUOrtState>();
    if (!fOrt->env) {
        fOrt->env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "rhonno_txmlpnpu");
        fOrt->mem_info = std::make_unique<Ort::MemoryInfo>(
            Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault));
    }

    Ort::SessionOptions so;
    so.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

    bool npu_ok = false;
    if (!fDisableNPU) {
        // VitisAI EP: needs config file path; vaip_config.json is part of the
        // Ryzen AI C++ install. Location may be overridden via env.
        const char* cfg = std::getenv("RHONNO_VAIP_CONFIG");
        std::string cfg_path;
        if (cfg && *cfg) {
            cfg_path = cfg;
        } else {
            const char* home = std::getenv("HOME");
            cfg_path = std::string(home ? home : "") +
                "/ryzen-ai-cpp/voe-4.0-linux_x86_64/vaip_config.json";
        }
        try {
            std::unordered_map<std::string, std::string> kv;
            kv["config_file"] = cfg_path;
            // VAIML compile cache. Key depends on the model fingerprint because
            // vaip reuses whatever was compiled under the key; otherwise we hit
            // stale images after Train() or when the CPU-fallback and NPU path
            // share the same cache.
            const char* cache_dir = std::getenv("RHONNO_VAIP_CACHE");
            const std::string cacheDir = (cache_dir && *cache_dir)
                ? std::string(cache_dir)
                : std::string(std::getenv("HOME") ? std::getenv("HOME") : ".")
                  + "/.cache/txmlpnpu_vaip";
            std::string model = BuildOnnxModel();
            size_t h = std::hash<std::string>{}(model);
            char keybuf[32];
            snprintf(keybuf, sizeof(keybuf), "txmlpnpu_%zx", h);
            kv["cacheDir"] = cacheDir;
            kv["cacheKey"] = keybuf;
            so.AppendExecutionProvider_VitisAI(kv);
            // Trial session to detect whether VitisAI actually accepts
            Ort::Session trial(*fOrt->env, model.data(), model.size(), so);
            fOrt->session = std::make_unique<Ort::Session>(std::move(trial));

            // EP-acceptance is not the same as actual NPU offload. Look at the
            // partition summary inside the VAIP cache.
            std::string det;
            bool offloaded = DetectVaimlNpuOffload(cacheDir, keybuf, &det);
            std::cerr << "TXMLPNPU: VAIML partition: " << det << std::endl;
            if (offloaded) {
                npu_ok = true;
            } else {
                std::cerr << "TXMLPNPU: no AIE subgraph — inference stays on CPU "
                             "despite VitisAI EP (expected for tiny/topologically "
                             "unsupported models)"
                          << std::endl;
                npu_ok = false;
            }
        } catch (const Ort::Exception& e) {
            std::cerr << "TXMLPNPU: VitisAI EP unavailable (" << e.what()
                      << ") — falling back to CPU EP" << std::endl;
            fOrt->session.reset();
            npu_ok = false;
        }
    }

    if (!fOrt->session) {
        Ort::SessionOptions so_cpu;
        so_cpu.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
        std::string model = BuildOnnxModel();
        fOrt->session = std::make_unique<Ort::Session>(*fOrt->env, model.data(), model.size(), so_cpu);
    }

    // Cache input/output names
    {
        Ort::AllocatorWithDefaultOptions alloc;
        auto in = fOrt->session->GetInputNameAllocated(0, alloc);
        auto out = fOrt->session->GetOutputNameAllocated(0, alloc);
        fInputName = in.get();
        fOutputName = out.get();
    }

    fOnNPU = npu_ok;
    fSessionValid = true;
}

void TXMLPNPU::ForwardFloat(const double* inScaled, double* out) const
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

double* TXMLPNPU::Recall(NNO_INTYPE* in, NNO_OUTTYPE* out)
{
    EnsureSession();

    // Scale input to net range (float32)
    std::vector<float> x(fParm.fInNodes);
    for (int i = 0; i < fParm.fInNodes; ++i)
        x[i] = static_cast<float>(static_cast<double>(in[i]) * fParm.fInScale);

    std::array<int64_t, 2> shape{1, static_cast<int64_t>(fParm.fInNodes)};
    Ort::Value input = Ort::Value::CreateTensor<float>(
        *fOrt->mem_info, x.data(), x.size(), shape.data(), shape.size());

    const char* in_names[] = { fInputName.c_str() };
    const char* out_names[] = { fOutputName.c_str() };

    auto outs = fOrt->session->Run(Ort::RunOptions{}, in_names, &input, 1,
                              out_names, 1);
    float* y = outs[0].GetTensorMutableData<float>();
    for (int o = 0; o < fParm.fOutNodes; ++o) fOut[o] = static_cast<double>(y[o]);

    if (fPlotter) {
        bool good = true;
        if (out != 0) good = out[0] > fParm.fThreshold;
        fPlotter->AddTestSample(fOut[0], good);
    }
    return fOut;
}

double TXMLPNPU::Train(NNO_INTYPE* in, NNO_OUTTYPE* trout)
{
    fShouldSave = true;
    InvalidateSession(); // float weights changed

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

// --- persistence (text only; binary is aliased to text for simplicity) -----

void TXMLPNPU::WriteText()
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

void TXMLPNPU::ReadText()
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
    InvalidateSession();
}

void TXMLPNPU::WriteBinary() { WriteText(); }
void TXMLPNPU::ReadBinary()  { ReadText(); }
