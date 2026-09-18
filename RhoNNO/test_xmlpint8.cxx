// test_xmlpint8 — smoke test for int8 TXMLP with ReLU
// Learns a simple linear target (ReLU-friendly), quantizes, checks Recall.

#include "TROOT.h"
#include "TXMLPInt8.h"
#include "TXMLP.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

// Target: y = 0.25*(x0+x1+x2+x3) - 0.5  in [-0.5, 0.5] for binary inputs
static float target_of(const NNO_INTYPE* x, int n)
{
    float s = 0.f;
    for (int i = 0; i < n; ++i) s += x[i];
    return 0.25f * s - 0.5f;
}

int main()
{
    const int IN = 4;
    const int EPOCHS = 120;
    const int SAMPLES = 200;

    cout << "TXMLPInt8 ReLU smoke test" << endl;
    srand(42);

    TXMLPInt8* qnet = new TXMLPInt8(
        3, 1.0, "xmlpint8_relu.net",
        IN, 16, 8, 1,
        0.02, 0.01, 0.005,
        /*reluOutput=*/false);
    qnet->SetThreshold(0.0);

    // Float TXMLP with ReLU hidden for conversion path
    TXMLP* fnet = new TXMLP(
        3, 1.0, "xmlp_relu_float.net",
        IN, 16, 8, 1,
        0.01, 0.005, 0.002,
        TNeuralNetParameters::TR_RELU,
        TNeuralNetParameters::TR_RELU,
        TNeuralNetParameters::TR_LINEAR);
    fnet->SetThreshold(0.0);

    NNO_INTYPE in[IN];
    NNO_OUTTYPE out[1];

    auto fill_random = [&]() {
        for (int i = 0; i < IN; ++i) in[i] = (rand() % 2) ? 1.f : 0.f;
        out[0] = target_of(in, IN);
    };

    double err_q = 0;
    for (int e = 0; e < EPOCHS; ++e) {
        err_q = 0;
        double err_f = 0;
        for (int s = 0; s < SAMPLES; ++s) {
            fill_random();
            err_q += qnet->Train(in, out);
            err_f += fnet->Train(in, out);
        }
        if ((e + 1) % 30 == 0)
            cout << "epoch " << (e + 1)
                 << "  int8_shadow_mse=" << err_q / SAMPLES
                 << "  float_mse=" << err_f / SAMPLES << endl;
    }

    qnet->Quantize();
    if (!qnet->IsQuantized()) {
        cerr << "FAIL: Quantize() did not set IsQuantized" << endl;
        return 1;
    }

    TXMLPInt8 from_float(*fnet, "xmlpint8_from_float.net");

    // Evaluate MSE on all 16 binary patterns
    double mse_q = 0, mse_f = 0, mse_c = 0, max_fq = 0;
    int npat = 1 << IN;
    for (int s = 0; s < npat; ++s) {
        for (int i = 0; i < IN; ++i) in[i] = ((s >> i) & 1) ? 1.f : 0.f;
        out[0] = target_of(in, IN);

        double* yq = qnet->Recall(in, out);
        double* yf = fnet->Recall(in, out);
        double* yc = from_float.Recall(in, out);

        mse_q += (yq[0] - out[0]) * (yq[0] - out[0]);
        mse_f += (yf[0] - out[0]) * (yf[0] - out[0]);
        mse_c += (yc[0] - out[0]) * (yc[0] - out[0]);

        double d = qnet->CompareFloatVsInt8(in);
        if (d > max_fq) max_fq = d;
    }
    mse_q /= npat;
    mse_f /= npat;
    mse_c /= npat;

    cout << "MSE TXMLPInt8 (int8 recall): " << mse_q << endl;
    cout << "MSE float TXMLP(ReLU):       " << mse_f << endl;
    cout << "MSE from quantized TXMLP:    " << mse_c << endl;
    cout << "max |float-int8| shadow:     " << max_fq << endl;
    cout << "layers=" << qnet->GetNLayers()
         << " L0 w_scale=" << qnet->GetLayer(0).fWScale
         << " x_scale=" << qnet->GetLayer(0).fXScale
         << " relu_h=" << qnet->GetLayer(0).fReLU
         << " relu_out=" << qnet->GetLayer(2).fReLU << endl;

    // -ffast-math can break std::isfinite; use x!=x for NaN
    if (mse_q != mse_q || mse_q > 0.02) {
        cerr << "FAIL: TXMLPInt8 MSE bad (" << mse_q << ")" << endl;
        return 2;
    }
    if (mse_f != mse_f) {
        cerr << "WARN: float ReLU TXMLP MSE is NaN (known fragile with -ffast-math / high LR)" << endl;
    } else if (mse_f > 0.08) {
        cerr << "WARN: float ReLU TXMLP MSE high (" << mse_f << ")" << endl;
    }
    if (max_fq != max_fq || max_fq > 0.15) {
        cerr << "FAIL: float vs int8 mismatch too large (" << max_fq << ")" << endl;
        return 4;
    }

    cout << "PASS" << endl;
    delete qnet;
    delete fnet;
    return 0;
}
