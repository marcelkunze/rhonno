#ifndef TXMLPGPU_H
#define TXMLPGPU_H
// TXMLPGPU
//
// Multi-Layer Perceptron that trains on CPU (float shadow weights, backprop,
// same dynamics as TXMLPNPU / TXMLPInt8) and executes Recall() on the AMD
// GPU (Radeon 890M, gfx1103) via HIP + hipBLAS.
//
// Design mirrors TXMLPNPU:
//   - training happens on CPU in float (identical semantics)
//   - Recall() runs layer-by-layer as hipblasSgemm + elementwise ReLU kernel
//   - GPU weights are uploaded lazily on first Recall and re-uploaded after
//     any Train() call (weights dirty flag)
//   - if HIP is unavailable at runtime, we silently fall back to the CPU
//     forward and SetGPU(false) becomes a no-op
//
// Part of the Neural Network Objects package (NNO)

#include "VSupervisedNet.h"
#include "TXMLP.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Pimpl holder — defined in TXMLPGPU.cxx (needs HIP at definition site)
struct TXMLPGPUGpuState;

class TXMLPGPU : public VSupervisedNet {
public:
    struct LayerF {
        int fInNodes = 0;
        int fOutNodes = 0;
        std::vector<double> fWFloat; // row-major [out][in]
        std::vector<double> fBFloat; // [out]
        bool fReLU = true;           // hidden layers: true
        double fLearnStep = 0.01;
    };

    TXMLPGPU();
    /// layers >= 1, nodes[0..layers-1] gives the layer widths
    TXMLPGPU(int layers, double inputRange, std::string netFile,
             int innodes, const std::vector<int>& nodes,
             const std::vector<double>& learnSteps,
             bool reluOutput = false);
    /// Convenience 3-layer ctor
    TXMLPGPU(int layers, double inputRange, std::string netFile,
             int innodes, int n0, int n1, int n2,
             double s0, double s1, double s2,
             bool reluOutput = false);
    /// Load from net file
    TXMLPGPU(std::string netFile);
    /// Build GPU net from a trained float TXMLP
    explicit TXMLPGPU(TXMLP& src, std::string netFile = "xmlpgpu.net");

    virtual ~TXMLPGPU();

    double Train(NNO_INTYPE* in, NNO_OUTTYPE* out) override;
    double* Recall(NNO_INTYPE* in, NNO_OUTTYPE* out = 0) override;

    void SetMomentumTerm(double f) override;

    /// True when the last Recall actually ran on GPU (weights uploaded and
    /// hipblas accepted the graph). False when we fell back to CPU.
    bool IsGPU() const { return fOnGPU; }

    /// Force CPU path even when HIP is available (testing/parity checks)
    void SetDisableGPU(bool v);

    /// True when compiled with HIP support
    static bool HasGpuSupport();

    int GetNLayers() const { return fParm.fLayers; }
    const LayerF& GetLayer(int i) const { return fLayers[i]; }
    LayerF& GetLayer(int i) { return fLayers[i]; }

    /// Max abs error of float CPU forward vs GPU forward on one sample
    double CompareCpuVsGpu(NNO_INTYPE* in);

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
    bool fDisableGPU = false;

    std::unique_ptr<TXMLPGPUGpuState> fGpu;
    bool fOnGPU = false;

    void BuildLayers(int layers, int innodes, const std::vector<int>& nodes,
                     const std::vector<double>& steps, bool reluOutput);
    void ForwardFloat(const double* inScaled, double* out) const;
    void InvalidateGpu();
};

#endif
