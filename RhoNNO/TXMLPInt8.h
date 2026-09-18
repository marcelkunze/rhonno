#ifndef TXMLPINT8_H
#define TXMLPINT8_H
// TXMLPInt8
//
// Int8-quantized Multi-Layer Perceptron with ReLU (hidden) / Linear (output).
// Part of the Neural Network Objects package (NNO)
//
// Training uses float shadow weights (same dynamics as TXMLP + TR_RELU).
// Call Quantize() (also done automatically after Train) to pack int8 weights
// and per-layer scales for fast int8 Recall().
//
// Forward (inference, after Quantize):
//   acc_i32 = sum_j W_q[i,j] * x_q[j] + B_i32
//   y = (acc_i32) * (w_scale * x_scale)
//   y = ReLU(y) for hidden layers
//   next x_q = round(y / next_x_scale) clamped to int8
//

#include "VSupervisedNet.h"
#include "TXMLP.h"

#include <cstdint>
#include <string>
#include <vector>

class TXMLPInt8 : public VSupervisedNet {
public:
    struct LayerQ {
        int fInNodes = 0;
        int fOutNodes = 0;
        std::vector<int8_t> fW;      // row-major [out][in]
        std::vector<int32_t> fB;     // bias in accumulator domain
        std::vector<double> fWFloat; // shadow float weights [out*in]
        std::vector<double> fBFloat; // shadow float biases [out]
        float fWScale = 1.f;         // real_w ≈ w_q * fWScale
        float fXScale = 1.f;         // real_x ≈ x_q * fXScale
        bool fReLU = true;           // hidden: true, last layer often false
        double fLearnStep = 0.01;
    };

    TXMLPInt8() {}
    /// layers >= 1. nodes[0..layers-1] sizes; last layer is output.
    /// Hidden layers use ReLU; output uses linear (no ReLU) by default.
    TXMLPInt8(int layers, double inputRange, std::string netFile,
              int innodes, const std::vector<int>& nodes,
              const std::vector<double>& learnSteps,
              bool reluOutput = false);
    /// Convenience 3-layer ctor (hidden-hidden-out or in practice 2 hidden + out)
    TXMLPInt8(int layers, double inputRange, std::string netFile,
              int innodes, int n0, int n1, int n2,
              double s0, double s1, double s2,
              bool reluOutput = false);
    /// Load from text/binary net file
    TXMLPInt8(std::string netFile);
    /// Build int8 net by quantizing a float TXMLP (ReLU recommended on float net)
    explicit TXMLPInt8(TXMLP& src, std::string netFile = "xmlpint8.net");

    virtual ~TXMLPInt8();

    double Train(NNO_INTYPE* in, NNO_OUTTYPE* out) override;
    double* Recall(NNO_INTYPE* in, NNO_OUTTYPE* out = 0) override;

    void SetMomentumTerm(double f) override;

    /// Pack float shadow weights into int8 + scales (symmetric per-layer)
    void Quantize();
    /// True after successful Quantize()
    bool IsQuantized() const { return fQuantized; }

    int GetNLayers() const { return fParm.fLayers; }
    const LayerQ& GetLayer(int i) const { return fLayers[i]; }
    LayerQ& GetLayer(int i) { return fLayers[i]; }

    /// Max abs error of float vs int8 forward on one sample (debug)
    double CompareFloatVsInt8(NNO_INTYPE* in);

protected:
    void AllocNet() override;
    void InitNet() override;
    void WriteText() override;
    void WriteBinary() override;
    void ReadText() override;
    void ReadBinary() override;

private:
    std::vector<LayerQ> fLayers;
    std::vector<double> fAct;      // float activations scratch
    std::vector<double> fDelta;    // backprop scratch
    std::vector<int8_t> fActQ;     // int8 activations scratch
    bool fQuantized = false;
    bool fReLUOutput = false;
    double fMomentum = 0.0;

    void BuildLayers(int layers, int innodes, const std::vector<int>& nodes,
                     const std::vector<double>& steps, bool reluOutput);
    void ForwardFloat(const double* inScaled, double* out);
    void ForwardInt8(const NNO_INTYPE* in, double* out);
    static int8_t QuantizeValue(double v, float scale);
    static float SymmetricScale(const double* data, int n, double eps = 1e-8);

    // No ClassDef: int8/STL members break ROOT auto-streamers; I/O via WriteText/Binary.
};

#endif
