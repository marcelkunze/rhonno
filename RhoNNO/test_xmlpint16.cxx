// test_xmlpint16 — int16 TXMLP + ReLU + ROOT streamer round-trip

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TClass.h"
#include "TXMLPInt16.h"
#include "TXMLP.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

static float target_of(const NNO_INTYPE* x, int n)
{
    float s = 0.f;
    for (int i = 0; i < n; ++i) s += x[i];
    return 0.25f * s - 0.5f;
}

int main()
{
    const int IN = 4;
    const int EPOCHS = 100;
    const int SAMPLES = 200;
    cout << "TXMLPInt16 ReLU + ROOT streamer smoke test" << endl;
    srand(7);

    // Dictionary must be loaded with lib
    if (!TClass::GetClass("TXMLPInt16")) {
        cerr << "FAIL: TClass TXMLPInt16 not found (dictionary missing)" << endl;
        return 10;
    }
    if (!TClass::GetClass("TXMLPInt16Layer")) {
        cerr << "FAIL: TClass TXMLPInt16Layer not found" << endl;
        return 11;
    }
    cout << "ROOT dict: TXMLPInt16 + TXMLPInt16Layer OK" << endl;

    TXMLPInt16* net = new TXMLPInt16(
        3, 1.0, "xmlpint16_relu.net",
        IN, 16, 8, 1,
        0.02, 0.01, 0.005,
        kFALSE);
    // Regression smoke test: linear output (not FERMI classification)
    net->SetOutputTransfer(TXMLPInt16::kLinear);
    net->SetThreshold(0.0);

    NNO_INTYPE in[IN];
    NNO_OUTTYPE out[1];
    auto fill = [&]() {
        for (int i = 0; i < IN; ++i) in[i] = (rand() % 2) ? 1.f : 0.f;
        out[0] = target_of(in, IN);
    };

    for (int e = 0; e < EPOCHS; ++e) {
        double err = 0;
        for (int s = 0; s < SAMPLES; ++s) {
            fill();
            err += net->Train(in, out);
        }
        if ((e + 1) % 25 == 0)
            cout << "epoch " << (e + 1) << " mse=" << err / SAMPLES << endl;
    }
    net->Quantize();
    if (!net->IsQuantized()) {
        cerr << "FAIL: not quantized" << endl;
        return 1;
    }

    // Evaluate before I/O
    double mse0 = 0;
    const int npat = 1 << IN;
    double y_before[16];
    for (int s = 0; s < npat; ++s) {
        for (int i = 0; i < IN; ++i) in[i] = ((s >> i) & 1) ? 1.f : 0.f;
        out[0] = target_of(in, IN);
        double* y = net->Recall(in, out);
        y_before[s] = y[0];
        mse0 += (y[0] - out[0]) * (y[0] - out[0]);
    }
    mse0 /= npat;
    cout << "MSE before ROOT I/O: " << mse0 << endl;
    if (mse0 != mse0 || mse0 > 0.02) {
        cerr << "FAIL: MSE too high before I/O" << endl;
        return 2;
    }

    // ROOT streamer round-trip
    const char* rpath = "xmlpint16_roundtrip.root";
    if (!net->WriteToRootFile(rpath, "net16")) {
        cerr << "FAIL: WriteToRootFile" << endl;
        return 3;
    }
    cout << "Wrote " << rpath << endl;

    TXMLPInt16* loaded = TXMLPInt16::ReadFromRootFile(rpath, "net16");
    if (!loaded) {
        cerr << "FAIL: ReadFromRootFile returned null" << endl;
        return 4;
    }
    if (loaded->GetNLayers() != net->GetNLayers()) {
        cerr << "FAIL: layer count mismatch " << loaded->GetNLayers()
             << " vs " << net->GetNLayers() << endl;
        return 5;
    }
    auto* L0 = loaded->GetLayer(0);
    if (!L0 || L0->fW.GetSize() < 1) {
        cerr << "FAIL: loaded layer0 weights empty" << endl;
        return 6;
    }
    cout << "Loaded layers=" << loaded->GetNLayers()
         << " L0 Wsize=" << L0->fW.GetSize()
         << " w_scale=" << L0->fWScale
         << " transfer=" << L0->fTransfer << endl;

    double mse1 = 0, maxdiff = 0;
    for (int s = 0; s < npat; ++s) {
        for (int i = 0; i < IN; ++i) in[i] = ((s >> i) & 1) ? 1.f : 0.f;
        out[0] = target_of(in, IN);
        double* y = loaded->Recall(in, out);
        mse1 += (y[0] - out[0]) * (y[0] - out[0]);
        double d = fabs(y[0] - y_before[s]);
        if (d > maxdiff) maxdiff = d;
    }
    mse1 /= npat;
    cout << "MSE after ROOT I/O:  " << mse1 << endl;
    cout << "max |y_load - y0|:   " << maxdiff << endl;

    if (mse1 != mse1 || mse1 > 0.02) {
        cerr << "FAIL: MSE after I/O bad" << endl;
        return 7;
    }
    if (maxdiff > 1e-6) {
        cerr << "FAIL: streamer changed outputs (diff=" << maxdiff << ")" << endl;
        return 8;
    }

    // Also exercise TFile direct Write of TObject
    {
        TFile f("xmlpint16_direct.root", "RECREATE");
        net->Write("direct");
        f.Close();
        TFile f2("xmlpint16_direct.root", "READ");
        auto* n2 = dynamic_cast<TXMLPInt16*>(f2.Get("direct"));
        if (!n2) {
            cerr << "FAIL: direct TFile::Get" << endl;
            return 9;
        }
        cout << "direct TFile GetNLayers=" << n2->GetNLayers() << endl;
        f2.Close();
    }

    gSystem->Unlink(rpath);
    gSystem->Unlink("xmlpint16_direct.root");

    cout << "PASS" << endl;
    delete loaded;
    delete net;
    return 0;
}
