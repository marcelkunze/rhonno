#ifndef TXMLPNPU_H
#define TXMLPNPU_H
// TXMLPNPU
//
// Multi-Layer Perceptron that trains on CPU (float shadow weights, backprop,
// same dynamics as TXMLPInt8) and executes Recall() on the AMD Ryzen XDNA NPU
// via ONNX Runtime with the VitisAI Execution Provider.
//
// The ONNX graph is built in-memory from the float shadow weights
// (no external onnx library or disk export needed) and cached; it is rebuilt
// whenever weights change (after Train()) or the EP disappears.
//
// If VitisAI EP / XRT / NPU are missing, code falls back to ORT
// CPUExecutionProvider automatically and sets IsNPU() = false. Behaviour of
// Recall() stays numerically identical (float32 MLP, ReLU on hidden layers).
//
// Part of the Neural Network Objects package (NNO)

#include "VSupervisedNet.h"
#include "TXMLP.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Pimpl holder — defined in TXMLPNPU.cxx (needs ORT headers at definition site)
struct TXMLPNPUOrtState;

class TXMLPNPU : public VSupervisedNet {
public:
    struct LayerF {
        int fInNodes = 0;
        int fOutNodes = 0;
        std::vector<double> fWFloat; // row-major [out][in]
        std::vector<double> fBFloat; // [out]
        bool fReLU = true;           // hidden layers: true
        double fLearnStep = 0.01;
    };

    TXMLPNPU();
    /// layers >= 1, nodes[0..layers-1] gives the layer widths
    TXMLPNPU(int layers, double inputRange, std::string netFile,
             int innodes, const std::vector<int>& nodes,
             const std::vector<double>& learnSteps,
             bool reluOutput = false);
    /// Convenience 3-layer ctor
    TXMLPNPU(int layers, double inputRange, std::string netFile,
             int innodes, int n0, int n1, int n2,
             double s0, double s1, double s2,
             bool reluOutput = false);
    /// Load from net file
    TXMLPNPU(std::string netFile);
    /// Build NPU net from a trained float TXMLP
    explicit TXMLPNPU(TXMLP& src, std::string netFile = "xmlpnpu.net");

    virtual ~TXMLPNPU();

    double Train(NNO_INTYPE* in, NNO_OUTTYPE* out) override;
    double* Recall(NNO_INTYPE* in, NNO_OUTTYPE* out = 0) override;

    void SetMomentumTerm(double f) override;

    /// True when the ORT session currently runs on VitisAI EP (NPU),
    /// false when it fell back to CPU.
    bool IsNPU() const { return fOnNPU; }

    /// Force CPU path even when NPU is available (testing/parity checks)
    void SetDisableNPU(bool v);

    int GetNLayers() const { return fParm.fLayers; }
    const LayerF& GetLayer(int i) const { return fLayers[i]; }
    LayerF& GetLayer(int i) { return fLayers[i]; }

    /// Build the ONNX model bytes from current float weights (exposed for tests)
    std::string BuildOnnxModel() const;
    /// Dump last built ONNX model to file (for offline inspection)
    bool DumpOnnx(const char* path) const;

protected:
    void AllocNet() override;
    void InitNet() override;
    void WriteText() override;
    void WriteBinary() override;
    void ReadText() override;
    void ReadBinary() override;

private:
    std::vector<LayerF> fLayers;
    bool fReLUOutput = false;
    double fMomentum = 0.0;
    bool fDisableNPU = false;

    // ORT state — built lazily on first Recall and whenever weights change
    std::unique_ptr<TXMLPNPUOrtState> fOrt;
    bool fSessionValid = false;
    bool fOnNPU = false;
    std::string fInputName;
    std::string fOutputName;

    void BuildLayers(int layers, int innodes, const std::vector<int>& nodes,
                     const std::vector<double>& steps, bool reluOutput);
    void ForwardFloat(const double* inScaled, double* out) const;
    void InvalidateSession();
    void EnsureSession();
};

#endif
