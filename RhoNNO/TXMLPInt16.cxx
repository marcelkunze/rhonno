// TXMLPInt16 — int16 MLP, ReLU hidden / FERMI output, ROOT streamers
// Part of the Neural Network Objects package (NNO)

#include "TXMLPInt16.h"
#include "VNeuralNetPlotter.h"
#include "TPerceptron.h"

#include "TFile.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ClassImp(TXMLPInt16Layer)
ClassImp(TXMLPInt16)

namespace {
constexpr Double_t kQMax = 32767.0;

inline Double_t clampf(Double_t x, Double_t lo, Double_t hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}
} // namespace

// ----- activation -----------------------------------------------------------

Double_t TXMLPInt16::Activate(Int_t transfer, Double_t x)
{
    switch (transfer) {
    case kReLU:
        return x > 0.0 ? x : 0.0;
    case kFermi: {
        x = clampf(x, -20.0, 20.0);
        return 1.0 / (1.0 + std::exp(-x));
    }
    case kLinear:
    default:
        return x;
    }
}

Double_t TXMLPInt16::ActivateDeriv(Int_t transfer, Double_t y)
{
    // y is post-activation value
    switch (transfer) {
    case kReLU:
        return y > 0.0 ? 1.0 : 0.0;
    case kFermi:
        return y * (1.0 - y);
    case kLinear:
    default:
        return 1.0;
    }
}

// ----- TXMLPInt16Layer -------------------------------------------------------

TXMLPInt16Layer::TXMLPInt16Layer()
    : TObject(), fInNodes(0), fOutNodes(0), fWScale(1.f), fXScale(1.f),
      fTransfer(TXMLPInt16::kReLU), fLearnStep(0.01)
{
}

TXMLPInt16Layer::TXMLPInt16Layer(Int_t inNodes, Int_t outNodes,
                                 Double_t learnStep, Int_t transfer)
    : TObject(), fInNodes(0), fOutNodes(0), fWScale(1.f), fXScale(1.f),
      fTransfer(transfer), fLearnStep(learnStep)
{
    Resize(inNodes, outNodes);
}

TXMLPInt16Layer::TXMLPInt16Layer(const TXMLPInt16Layer& o)
    : TObject(o), fInNodes(o.fInNodes), fOutNodes(o.fOutNodes),
      fW(o.fW), fB(o.fB), fWFloat(o.fWFloat), fBFloat(o.fBFloat),
      fDW(o.fDW), fDB(o.fDB), fWScale(o.fWScale), fXScale(o.fXScale),
      fTransfer(o.fTransfer), fLearnStep(o.fLearnStep)
{
}

TXMLPInt16Layer& TXMLPInt16Layer::operator=(const TXMLPInt16Layer& o)
{
    if (this == &o) return *this;
    TObject::operator=(o);
    fInNodes = o.fInNodes;
    fOutNodes = o.fOutNodes;
    fW = o.fW;
    fB = o.fB;
    fWFloat = o.fWFloat;
    fBFloat = o.fBFloat;
    fDW = o.fDW;
    fDB = o.fDB;
    fWScale = o.fWScale;
    fXScale = o.fXScale;
    fTransfer = o.fTransfer;
    fLearnStep = o.fLearnStep;
    return *this;
}

void TXMLPInt16Layer::Resize(Int_t inNodes, Int_t outNodes)
{
    fInNodes = inNodes;
    fOutNodes = outNodes;
    const Int_t nw = inNodes * outNodes;
    fW.Set(nw);
    fB.Set(outNodes);
    fWFloat.Set(nw);
    fBFloat.Set(outNodes);
    fDW.Set(nw);
    fDB.Set(outNodes);
    ClearWeights();
}

void TXMLPInt16Layer::ClearWeights()
{
    if (fW.fArray) std::memset(fW.fArray, 0, sizeof(Short_t) * fW.GetSize());
    if (fB.fArray) std::memset(fB.fArray, 0, sizeof(Int_t) * fB.GetSize());
    auto zeron = [](TArrayD& a) {
        if (a.fArray)
            for (Int_t i = 0; i < a.GetSize(); ++i) a.fArray[i] = 0.0;
    };
    zeron(fWFloat);
    zeron(fBFloat);
    zeron(fDW);
    zeron(fDB);
}

// ----- helpers ---------------------------------------------------------------

Short_t TXMLPInt16::QuantizeValue(Double_t v, Float_t scale)
{
    if (!(scale > 0.f) || !(v == v) || !std::isfinite(v)) return 0;
    Double_t q = std::round(v / static_cast<Double_t>(scale));
    if (q > kQMax) q = kQMax;
    if (q < -kQMax) q = -kQMax;
    return static_cast<Short_t>(q);
}

Float_t TXMLPInt16::SymmetricScale(const Double_t* data, Int_t n, Double_t eps)
{
    Double_t amax = 0.0;
    for (Int_t i = 0; i < n; ++i) {
        Double_t a = std::fabs(data[i]);
        if (a > amax) amax = a;
    }
    if (amax < eps) amax = eps;
    return static_cast<Float_t>(amax / kQMax);
}

void TXMLPInt16::ClearLayers()
{
    fLayers.Delete();
    fLayers.Clear();
}

void TXMLPInt16::EnsureScratch()
{
    Int_t maxw = fParm.fInNodes;
    for (Int_t L = 0; L < fLayers.GetEntriesFast(); ++L) {
        auto* ly = static_cast<TXMLPInt16Layer*>(fLayers.At(L));
        if (ly && ly->fOutNodes > maxw) maxw = ly->fOutNodes;
    }
    if (fAct.GetSize() < maxw + 8) fAct.Set(maxw + 8);
    if (fDelta.GetSize() < maxw + 8) fDelta.Set(maxw + 8);
    if (fActQ.GetSize() < maxw + 8) fActQ.Set(maxw + 8);
}

void TXMLPInt16::BuildLayers(Int_t layers, Int_t innodes,
                             const Int_t* nodes, const Double_t* steps,
                             Bool_t reluOutput)
{
    if (layers < 1) Errorf((char*)"(TXMLPInt16) at least one layer necessary");
    if (!nodes || !steps) Errorf((char*)"(TXMLPInt16) null nodes/steps");

    ClearLayers();
    fLayers.SetOwner(kTRUE);

    fParm.fLayers = layers;
    fParm.fInNodes = innodes;
    fParm.fOutNodes = nodes[layers - 1];
    fParm.fTransferId = TNeuralNetParameters::TR_RELU;
    fReLUOutput = reluOutput;
    fQuantized = kFALSE;
    // Classification default: FERMI output (threshold 0.5). ReLU-out only if requested.
    fOutTransfer = reluOutput ? kReLU : kFermi;

    Int_t prev = innodes;
    for (Int_t L = 0; L < layers; ++L) {
        Int_t tr = kReLU;
        if (L == layers - 1) tr = fOutTransfer;
        auto* ly = new TXMLPInt16Layer(prev, nodes[L], steps[L], tr);
        fLayers.Add(ly);
        prev = nodes[L];
    }

    if (fOut) {
        delete[] fOut;
        fOut = nullptr;
    }
    fOut = new Double_t[fParm.fOutNodes];
    TestPointer(fOut);
    EnsureScratch();
}

void TXMLPInt16::SetOutputTransfer(Int_t t)
{
    fOutTransfer = t;
    fReLUOutput = (t == kReLU);
    const Int_t nL = fLayers.GetEntriesFast();
    if (nL > 0) GetLayer(nL - 1)->fTransfer = t;
}

// ----- ctors -----------------------------------------------------------------

TXMLPInt16::TXMLPInt16()
    : VSupervisedNet(), fLayers(0), fQuantized(kFALSE), fReLUOutput(kFALSE),
      fUseInt16Recall(kFALSE), fOutTransfer(kFermi), fMomentum(0.0)
{
    fLayers.SetOwner(kTRUE);
}

TXMLPInt16::TXMLPInt16(Int_t layers, Double_t inputRange, std::string netFile,
                       Int_t innodes, Int_t n0, Int_t n1, Int_t n2,
                       Double_t s0, Double_t s1, Double_t s2, Bool_t reluOutput)
    : VSupervisedNet("XMLP16", innodes, 0, netFile), fLayers(0),
      fQuantized(kFALSE), fReLUOutput(reluOutput), fUseInt16Recall(kFALSE),
      fOutTransfer(reluOutput ? kReLU : kFermi), fMomentum(0.0)
{
    if (layers != 3) Errorf((char*)"(TXMLPInt16) 3-node ctor needs layers==3");
    fParm.fInScale = 1.0 / inputRange;
    fLayers.SetOwner(kTRUE);
    const Int_t nodes[3] = {n0, n1, n2};
    const Double_t steps[3] = {s0, s1, s2};
    BuildLayers(3, innodes, nodes, steps, reluOutput);
    InitNet();
}

TXMLPInt16::TXMLPInt16(Int_t layers, Double_t inputRange, std::string netFile,
                       Int_t innodes, const Int_t* nodes, const Double_t* learnSteps,
                       Bool_t reluOutput)
    : VSupervisedNet("XMLP16", innodes, 0, netFile), fLayers(0),
      fQuantized(kFALSE), fReLUOutput(reluOutput), fUseInt16Recall(kFALSE),
      fOutTransfer(reluOutput ? kReLU : kFermi), fMomentum(0.0)
{
    fParm.fInScale = 1.0 / inputRange;
    fLayers.SetOwner(kTRUE);
    BuildLayers(layers, innodes, nodes, learnSteps, reluOutput);
    InitNet();
}

TXMLPInt16::TXMLPInt16(std::string netFile)
    : VSupervisedNet(netFile), fLayers(0), fQuantized(kFALSE),
      fReLUOutput(kFALSE), fUseInt16Recall(kFALSE), fOutTransfer(kFermi),
      fMomentum(0.0)
{
    fLayers.SetOwner(kTRUE);
    ReadNet("XMLP16");
}

TXMLPInt16::TXMLPInt16(TXMLP& src, std::string netFile)
    : VSupervisedNet("XMLP16", src.GetParameters().fInNodes,
                     src.GetParameters().fOutNodes, netFile),
      fLayers(0), fQuantized(kFALSE), fReLUOutput(kFALSE),
      fUseInt16Recall(kFALSE), fOutTransfer(kFermi), fMomentum(0.0)
{
    fLayers.SetOwner(kTRUE);
    const TNeuralNetParameters& sp = src.GetParameters();
    fParm.fInScale = sp.fInScale;
    fParm.fLayers = sp.fLayers;
    fParm.fInNodes = sp.fInNodes;
    fParm.fOutNodes = sp.fOutNodes;
    fParm.fTransferId = TNeuralNetParameters::TR_RELU;

    std::vector<Int_t> nodes(sp.fLayers);
    std::vector<Double_t> steps(sp.fLayers, 0.01);
    for (Int_t L = 0; L < sp.fLayers; ++L) {
        TPerceptron* p = src.GetPerceptron(L);
        nodes[L] = p->GetParameters().fOutNodes;
        steps[L] = p->GetParameters().fLearnStep;
    }
    BuildLayers(sp.fLayers, sp.fInNodes, nodes.data(), steps.data(), kFALSE);

    for (Int_t L = 0; L < sp.fLayers; ++L) {
        TPerceptron* p = src.GetPerceptron(L);
        auto* ly = GetLayer(L);
        const auto tid = p->GetParameters().fTransferId;
        if (tid == TNeuralNetParameters::TR_RELU) ly->fTransfer = kReLU;
        else if (tid == TNeuralNetParameters::TR_FERMI) ly->fTransfer = kFermi;
        else if (tid == TNeuralNetParameters::TR_LINEAR) ly->fTransfer = kLinear;
        else ly->fTransfer = (L == sp.fLayers - 1) ? kFermi : kReLU;
        if (L == sp.fLayers - 1) {
            fOutTransfer = ly->fTransfer;
            fReLUOutput = (ly->fTransfer == kReLU);
        }
        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            PerceptronUnit* up = &p->fU[o];
            ly->fBFloat.fArray[o] = -up->fThreshold;
            for (Int_t i = 0; i < ly->fInNodes; ++i)
                ly->fWFloat.fArray[o * ly->fInNodes + i] = up->fVector[i];
        }
    }
    Quantize();
}

TXMLPInt16::TXMLPInt16(const TXMLPInt16& o)
    : VSupervisedNet("XMLP16", o.fParm.fInNodes, o.fParm.fOutNodes, o.fFilename),
      fLayers(0), fQuantized(o.fQuantized), fReLUOutput(o.fReLUOutput),
      fUseInt16Recall(o.fUseInt16Recall), fOutTransfer(o.fOutTransfer),
      fMomentum(o.fMomentum)
{
    fLayers.SetOwner(kTRUE);
    fParm = o.fParm;
    fFilename = o.fFilename;
    for (Int_t i = 0; i < o.fLayers.GetEntriesFast(); ++i) {
        auto* src = static_cast<TXMLPInt16Layer*>(o.fLayers.At(i));
        fLayers.Add(new TXMLPInt16Layer(*src));
    }
    if (fParm.fOutNodes > 0) {
        fOut = new Double_t[fParm.fOutNodes];
        if (o.fOut) std::memcpy(fOut, o.fOut, sizeof(Double_t) * fParm.fOutNodes);
    }
    EnsureScratch();
}

TXMLPInt16& TXMLPInt16::operator=(const TXMLPInt16& o)
{
    if (this == &o) return *this;
    fParm = o.fParm;
    fFilename = o.fFilename;
    fQuantized = o.fQuantized;
    fReLUOutput = o.fReLUOutput;
    fUseInt16Recall = o.fUseInt16Recall;
    fOutTransfer = o.fOutTransfer;
    fMomentum = o.fMomentum;
    ClearLayers();
    fLayers.SetOwner(kTRUE);
    for (Int_t i = 0; i < o.fLayers.GetEntriesFast(); ++i) {
        auto* src = static_cast<TXMLPInt16Layer*>(o.fLayers.At(i));
        fLayers.Add(new TXMLPInt16Layer(*src));
    }
    if (fOut) {
        delete[] fOut;
        fOut = nullptr;
    }
    if (fParm.fOutNodes > 0) {
        fOut = new Double_t[fParm.fOutNodes];
        if (o.fOut) std::memcpy(fOut, o.fOut, sizeof(Double_t) * fParm.fOutNodes);
    }
    EnsureScratch();
    return *this;
}

TXMLPInt16::~TXMLPInt16()
{
    if (fFilename != "" && fShouldSave) Save();
    ClearLayers();
}

TXMLPInt16Layer* TXMLPInt16::GetLayer(Int_t i)
{
    return static_cast<TXMLPInt16Layer*>(fLayers.At(i));
}

const TXMLPInt16Layer* TXMLPInt16::GetLayer(Int_t i) const
{
    return static_cast<const TXMLPInt16Layer*>(fLayers.At(i));
}

void TXMLPInt16::AllocNet()
{
    std::vector<Int_t> nodes(std::max(1, fParm.fLayers), 1);
    std::vector<Double_t> steps(std::max(1, fParm.fLayers), 0.01);
    if (fParm.fLayers < 1) fParm.fLayers = 1;
    if (fParm.fInNodes < 1) fParm.fInNodes = 1;
    BuildLayers(fParm.fLayers, fParm.fInNodes, nodes.data(), steps.data(), fReLUOutput);
}

void TXMLPInt16::InitNet()
{
    for (Int_t L = 0; L < fLayers.GetEntriesFast(); ++L) {
        auto* ly = GetLayer(L);
        // He for ReLU, Xavier-ish for FERMI/linear
        const Double_t sigma =
            (ly->fTransfer == kReLU)
                ? std::sqrt(2.0 / std::max(1, ly->fInNodes))
                : std::sqrt(1.0 / std::max(1, ly->fInNodes));
        for (Int_t k = 0; k < ly->fWFloat.GetSize(); ++k) {
            Double_t u1 = 0.0, u2 = 0.0;
            while (u1 < 1e-12) u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
            u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
            Double_t z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            ly->fWFloat.fArray[k] = z * sigma;
        }
        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            ly->fBFloat.fArray[o] = 0.0;
            ly->fDB.fArray[o] = 0.0;
        }
        for (Int_t k = 0; k < ly->fDW.GetSize(); ++k) ly->fDW.fArray[k] = 0.0;
    }
    fQuantized = kFALSE;
}

void TXMLPInt16::SetMomentumTerm(Double_t f)
{
    fMomentum = f;
    fParm.fMu = f;
}

// ----- forward ---------------------------------------------------------------

void TXMLPInt16::Quantize()
{
    // Input assumed roughly in [-1,1] after fInScale / autoscale
    Float_t xscale = static_cast<Float_t>(1.0 / kQMax);

    for (Int_t L = 0; L < fLayers.GetEntriesFast(); ++L) {
        auto* ly = GetLayer(L);
        ly->fXScale = xscale;
        ly->fWScale = SymmetricScale(ly->fWFloat.fArray, ly->fWFloat.GetSize());

        if (ly->fW.GetSize() != ly->fWFloat.GetSize()) ly->fW.Set(ly->fWFloat.GetSize());
        for (Int_t k = 0; k < ly->fWFloat.GetSize(); ++k)
            ly->fW.fArray[k] = QuantizeValue(ly->fWFloat.fArray[k], ly->fWScale);

        // Bound next activation by |b| + sum|w| * xmax_in  (xmax_in = kQMax*xscale)
        const Double_t xmax_in = kQMax * ly->fXScale;
        Double_t amax = 0.0;
        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            Double_t s = std::fabs(ly->fBFloat.fArray[o]);
            for (Int_t i = 0; i < ly->fInNodes; ++i)
                s += std::fabs(ly->fWFloat.fArray[o * ly->fInNodes + i]) * xmax_in;
            // FERMI outputs in (0,1): bound next scale by 1
            if (ly->fTransfer == kFermi) s = std::min(s, 1.0);
            if (s > amax) amax = s;
        }
        if (amax < 1e-8) amax = 1e-8;
        xscale = static_cast<Float_t>(amax / kQMax);

        if (ly->fB.GetSize() != ly->fOutNodes) ly->fB.Set(ly->fOutNodes);
        const Double_t b_scale = static_cast<Double_t>(ly->fWScale) * ly->fXScale;
        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            Double_t bq = (b_scale > 0.0) ? std::round(ly->fBFloat.fArray[o] / b_scale) : 0.0;
            bq = clampf(bq,
                        static_cast<Double_t>(std::numeric_limits<Int_t>::min()),
                        static_cast<Double_t>(std::numeric_limits<Int_t>::max()));
            ly->fB.fArray[o] = static_cast<Int_t>(bq);
        }
    }
    fQuantized = kTRUE;
}

void TXMLPInt16::ForwardFloat(const Double_t* inScaled, Double_t* out)
{
    EnsureScratch();
    std::vector<Double_t> cur(inScaled, inScaled + fParm.fInNodes);
    std::vector<Double_t> nxt;

    for (Int_t L = 0; L < fLayers.GetEntriesFast(); ++L) {
        const auto* ly = GetLayer(L);
        nxt.assign(ly->fOutNodes, 0.0);
        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            Double_t sum = ly->fBFloat.fArray[o];
            const Double_t* w = &ly->fWFloat.fArray[o * ly->fInNodes];
            for (Int_t i = 0; i < ly->fInNodes; ++i) sum += w[i] * cur[i];
            nxt[o] = Activate(ly->fTransfer, sum);
        }
        cur.swap(nxt);
    }
    for (Int_t o = 0; o < fParm.fOutNodes; ++o) out[o] = cur[o];
}

void TXMLPInt16::ForwardInt16(const NNO_INTYPE* in, Double_t* out)
{
    if (!fQuantized) Quantize();
    EnsureScratch();

    const Int_t nL = fLayers.GetEntriesFast();
    if (nL < 1) return;

    const auto* l0 = GetLayer(0);
    std::vector<Short_t> cur(fParm.fInNodes);
    for (Int_t i = 0; i < fParm.fInNodes; ++i) {
        Double_t xs = static_cast<Double_t>(in[i]) * fParm.fInScale;
        cur[i] = QuantizeValue(xs, l0->fXScale);
    }

    std::vector<Short_t> nxt;
    std::vector<Double_t> yf;

    for (Int_t L = 0; L < nL; ++L) {
        const auto* ly = GetLayer(L);
        yf.resize(ly->fOutNodes);
        const Double_t deq = static_cast<Double_t>(ly->fWScale) * ly->fXScale;
        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            Long64_t acc = ly->fB.fArray[o];
            const Short_t* w = &ly->fW.fArray[o * ly->fInNodes];
            for (Int_t i = 0; i < ly->fInNodes; ++i)
                acc += static_cast<Long64_t>(w[i]) * static_cast<Long64_t>(cur[i]);
            Double_t y = Activate(ly->fTransfer, static_cast<Double_t>(acc) * deq);
            yf[o] = y;
        }
        if (L + 1 < nL) {
            const Float_t next_xs = GetLayer(L + 1)->fXScale;
            nxt.resize(ly->fOutNodes);
            for (Int_t o = 0; o < ly->fOutNodes; ++o)
                nxt[o] = QuantizeValue(yf[o], next_xs);
            cur.swap(nxt);
        } else {
            for (Int_t o = 0; o < ly->fOutNodes; ++o) out[o] = yf[o];
        }
    }
}

Double_t* TXMLPInt16::Recall(NNO_INTYPE* in, NNO_OUTTYPE* out)
{
    if (fUseInt16Recall) {
        ForwardInt16(in, fOut);
    } else {
        std::vector<Double_t> x(fParm.fInNodes);
        for (Int_t i = 0; i < fParm.fInNodes; ++i)
            x[i] = static_cast<Double_t>(in[i]) * fParm.fInScale;
        ForwardFloat(x.data(), fOut);
    }
    if (fPlotter) {
        Bool_t good = kTRUE;
        if (out != 0) good = out[0] > fParm.fThreshold;
        fPlotter->AddTestSample(fOut[0], good);
    }
    return fOut;
}

Double_t* TXMLPInt16::RecallInt16(NNO_INTYPE* in, NNO_OUTTYPE* out)
{
    ForwardInt16(in, fOut);
    if (fPlotter) {
        Bool_t good = kTRUE;
        if (out != 0) good = out[0] > fParm.fThreshold;
        fPlotter->AddTestSample(fOut[0], good);
    }
    return fOut;
}

// ----- train -----------------------------------------------------------------

Double_t TXMLPInt16::Train(NNO_INTYPE* in, NNO_OUTTYPE* trout)
{
    fShouldSave = kTRUE;
    fQuantized = kFALSE;

    const Int_t nL = fLayers.GetEntriesFast();
    std::vector<Double_t> x(fParm.fInNodes);
    for (Int_t i = 0; i < fParm.fInNodes; ++i)
        x[i] = static_cast<Double_t>(in[i]) * fParm.fInScale;

    // acts[L] = layer L input; acts[L+1] = layer L post-activation
    std::vector<std::vector<Double_t>> acts(nL + 1);
    acts[0] = x;
    for (Int_t L = 0; L < nL; ++L) {
        auto* ly = GetLayer(L);
        acts[L + 1].assign(ly->fOutNodes, 0.0);
        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            Double_t sum = ly->fBFloat.fArray[o];
            const Double_t* w = &ly->fWFloat.fArray[o * ly->fInNodes];
            const Double_t* xi = acts[L].data();
            for (Int_t i = 0; i < ly->fInNodes; ++i) sum += w[i] * xi[i];
            acts[L + 1][o] = Activate(ly->fTransfer, sum);
        }
    }

    for (Int_t o = 0; o < fParm.fOutNodes; ++o) fOut[o] = acts.back()[o];

    Double_t S_Err = 0.0;
    std::vector<Double_t> delta(fParm.fOutNodes);
    // Match classic TXMLP: last-layer fDiffSrc is overwritten with (t-y) ONLY
    // (no output transfer'). Hidden layers still multiply by transfer' in backprop.
    // Using FERMI' on the output crushed gradients and left scores in a [0,0.5] box.
    for (Int_t o = 0; o < fParm.fOutNodes; ++o) {
        Double_t diff = static_cast<Double_t>(trout[o]) - fOut[o];
        S_Err += diff * diff;
        delta[o] = diff;
    }

    const Double_t mu = fMomentum;

    for (Int_t L = nL - 1; L >= 0; --L) {
        auto* ly = GetLayer(L);
        const Double_t lr = ly->fLearnStep;
        std::vector<Double_t> delta_prev(ly->fInNodes, 0.0);
        // snapshot W for backprop (TPerceptron uses pre-update weights)
        std::vector<Double_t> wcopy(ly->fWFloat.fArray,
                                    ly->fWFloat.fArray + ly->fWFloat.GetSize());

        for (Int_t o = 0; o < ly->fOutNodes; ++o) {
            Double_t d = delta[o];
            if (!(d == d) || !std::isfinite(d)) d = 0.0;
            d = clampf(d, -5.0, 5.0);

            // momentum on bias (TPerceptron: threshold -= d*lr  <=>  bias += d*lr)
            Double_t db = d + mu * ly->fDB.fArray[o];
            ly->fDB.fArray[o] = db;
            ly->fBFloat.fArray[o] += db * lr;

            Double_t* w = &ly->fWFloat.fArray[o * ly->fInNodes];
            Double_t* dw = &ly->fDW.fArray[o * ly->fInNodes];
            const Double_t* w0 = &wcopy[o * ly->fInNodes];
            const Double_t* xi = acts[L].data();
            for (Int_t i = 0; i < ly->fInNodes; ++i) {
                Double_t g = xi[i] * d;
                Double_t m = g + mu * dw[i];
                dw[i] = m;
                w[i] += m * lr;
                w[i] = clampf(w[i], -50.0, 50.0);
                delta_prev[i] += w0[i] * d;
            }
        }

        if (L > 0) {
            // Like TPerceptron: prev DiffSrc was transfer'(act); multiply by backprop sum
            auto* prev = GetLayer(L - 1);
            for (Int_t i = 0; i < ly->fInNodes; ++i)
                delta_prev[i] *= ActivateDeriv(prev->fTransfer, acts[L][i]);
            delta.swap(delta_prev);
        }
    }

    // Classic NNO train pad: fill TARGET labels (0 or 1) → sharp peaks at 0 and 1
    // (score distribution is shown on the validation pad via Recall → AddTestSample)
    if (fPlotter)
        fPlotter->AddTrainSample(trout[0], trout[0] > fParm.fThreshold);
    return S_Err;
}

Double_t TXMLPInt16::CompareFloatVsInt16(NNO_INTYPE* in)
{
    std::vector<Double_t> xf(fParm.fInNodes), yf(fParm.fOutNodes), yq(fParm.fOutNodes);
    for (Int_t i = 0; i < fParm.fInNodes; ++i)
        xf[i] = static_cast<Double_t>(in[i]) * fParm.fInScale;
    ForwardFloat(xf.data(), yf.data());
    Bool_t was = fQuantized;
    Quantize();
    ForwardInt16(in, yq.data());
    fQuantized = was;
    Double_t maxd = 0.0;
    for (Int_t o = 0; o < fParm.fOutNodes; ++o)
        maxd = std::max(maxd, std::fabs(yf[o] - yq[o]));
    return maxd;
}

Bool_t TXMLPInt16::WriteToRootFile(const char* path, const char* name)
{
    if (!fQuantized) Quantize();
    TFile f(path, "RECREATE");
    if (f.IsZombie()) return kFALSE;
    Int_t nbytes = Write(name, TObject::kOverwrite);
    f.Write();
    f.Close();
    return nbytes > 0;
}

TXMLPInt16* TXMLPInt16::ReadFromRootFile(const char* path, const char* name)
{
    TFile f(path, "READ");
    if (f.IsZombie()) return nullptr;
    auto* net = dynamic_cast<TXMLPInt16*>(f.Get(name));
    if (!net) {
        f.Close();
        return nullptr;
    }
    auto* copy = static_cast<TXMLPInt16*>(net->Clone(name));
    f.Close();
    if (copy) {
        copy->EnsureScratch();
        if (!copy->fOut && copy->fParm.fOutNodes > 0)
            copy->fOut = new Double_t[copy->fParm.fOutNodes];
        copy->fLayers.SetOwner(kTRUE);
    }
    return copy;
}

// ----- text/binary I/O -------------------------------------------------------

void TXMLPInt16::WriteText()
{
    fprintf(fFile, "layers    %i\n", fParm.fLayers);
    fprintf(fFile, "in_scale  %le\n", fParm.fInScale);
    fprintf(fFile, "in_nodes  %i\n", fParm.fInNodes);
    fprintf(fFile, "out_nodes %i\n", fParm.fOutNodes);
    fprintf(fFile, "out_xfer  %i\n", fOutTransfer);
    fprintf(fFile, "momentum  %le\n", fMomentum);
    if (!fQuantized) Quantize();
    for (Int_t L = 0; L < fLayers.GetEntriesFast(); ++L) {
        const auto* ly = GetLayer(L);
        fprintf(fFile, "\nlayer %i\n", L);
        fprintf(fFile, "innodes     %i\n", ly->fInNodes);
        fprintf(fFile, "outnodes    %i\n", ly->fOutNodes);
        fprintf(fFile, "learn_step  %le\n", ly->fLearnStep);
        fprintf(fFile, "transfer    %i\n", ly->fTransfer);
        fprintf(fFile, "w_scale     %a\n", ly->fWScale);
        fprintf(fFile, "x_scale     %a\n", ly->fXScale);
        fprintf(fFile, "weights_q\n");
        for (Int_t k = 0; k < ly->fW.GetSize(); ++k)
            fprintf(fFile, "%d\n", static_cast<int>(ly->fW.fArray[k]));
        fprintf(fFile, "bias_q\n");
        for (Int_t o = 0; o < ly->fOutNodes; ++o)
            fprintf(fFile, "%d\n", ly->fB.fArray[o]);
        fprintf(fFile, "weights_f\n");
        for (Int_t k = 0; k < ly->fWFloat.GetSize(); ++k)
            fprintf(fFile, "%le\n", ly->fWFloat.fArray[k]);
        fprintf(fFile, "bias_f\n");
        for (Int_t o = 0; o < ly->fOutNodes; ++o)
            fprintf(fFile, "%le\n", ly->fBFloat.fArray[o]);
    }
}

void TXMLPInt16::ReadText()
{
    Int_t layers = 0, in_nodes = 0, out_nodes = 0, out_xfer = kFermi;
    Double_t scale = 1.0, mom = 0.0;
    fscanf(fFile, "layers    %i\n", &layers);
    fscanf(fFile, "in_scale  %le\n", &scale);
    fscanf(fFile, "in_nodes  %i\n", &in_nodes);
    fscanf(fFile, "out_nodes %i\n", &out_nodes);
    // backward compatible: old files had relu_out
    long pos = ftell(fFile);
    if (fscanf(fFile, "out_xfer  %i\n", &out_xfer) != 1) {
        fseek(fFile, pos, SEEK_SET);
        Int_t relu_out = 0;
        fscanf(fFile, "relu_out  %i\n", &relu_out);
        out_xfer = relu_out ? kReLU : kFermi;
    }
    fscanf(fFile, "momentum  %le\n", &mom);
    fParm.fLayers = layers;
    fParm.fInScale = scale;
    fParm.fInNodes = in_nodes;
    fParm.fOutNodes = out_nodes;
    fOutTransfer = out_xfer;
    fReLUOutput = (out_xfer == kReLU);
    fMomentum = mom;

    ClearLayers();
    fLayers.SetOwner(kTRUE);
    for (Int_t L = 0; L < layers; ++L) {
        Int_t lid = 0, inn = 0, outn = 0, tr = kReLU;
        Double_t ls = 0.01;
        Float_t ws = 1.f, xs = 1.f;
        fscanf(fFile, "\nlayer %i\n", &lid);
        fscanf(fFile, "innodes     %i\n", &inn);
        fscanf(fFile, "outnodes    %i\n", &outn);
        fscanf(fFile, "learn_step  %le\n", &ls);
        pos = ftell(fFile);
        if (fscanf(fFile, "transfer    %i\n", &tr) != 1) {
            fseek(fFile, pos, SEEK_SET);
            Int_t relu = 1;
            fscanf(fFile, "relu        %i\n", &relu);
            tr = relu ? kReLU : kLinear;
            if (L == layers - 1) tr = out_xfer;
        }
        fscanf(fFile, "w_scale     %a\n", &ws);
        fscanf(fFile, "x_scale     %a\n", &xs);
        auto* ly = new TXMLPInt16Layer(inn, outn, ls, tr);
        ly->fWScale = ws;
        ly->fXScale = xs;
        fscanf(fFile, "weights_q\n");
        for (Int_t k = 0; k < ly->fW.GetSize(); ++k) {
            int v = 0;
            fscanf(fFile, "%d\n", &v);
            ly->fW.fArray[k] = static_cast<Short_t>(v);
        }
        fscanf(fFile, "bias_q\n");
        for (Int_t o = 0; o < outn; ++o) {
            int v = 0;
            fscanf(fFile, "%d\n", &v);
            ly->fB.fArray[o] = v;
        }
        fscanf(fFile, "weights_f\n");
        for (Int_t k = 0; k < ly->fWFloat.GetSize(); ++k)
            fscanf(fFile, "%le\n", &ly->fWFloat.fArray[k]);
        fscanf(fFile, "bias_f\n");
        for (Int_t o = 0; o < outn; ++o)
            fscanf(fFile, "%le\n", &ly->fBFloat.fArray[o]);
        fLayers.Add(ly);
    }
    if (fOut) {
        delete[] fOut;
        fOut = nullptr;
    }
    fOut = new Double_t[fParm.fOutNodes];
    fQuantized = kTRUE;
    EnsureScratch();
}

void TXMLPInt16::WriteBinary()
{
    if (!fQuantized) Quantize();
    fwrite(&fParm, sizeof(TNeuralNetParameters), 1, fFile);
    Int_t ot = fOutTransfer;
    Int_t nL = fLayers.GetEntriesFast();
    fwrite(&ot, sizeof(Int_t), 1, fFile);
    fwrite(&fMomentum, sizeof(Double_t), 1, fFile);
    fwrite(&nL, sizeof(Int_t), 1, fFile);
    for (Int_t L = 0; L < nL; ++L) {
        auto* ly = GetLayer(L);
        fwrite(&ly->fInNodes, sizeof(Int_t), 1, fFile);
        fwrite(&ly->fOutNodes, sizeof(Int_t), 1, fFile);
        fwrite(&ly->fLearnStep, sizeof(Double_t), 1, fFile);
        fwrite(&ly->fTransfer, sizeof(Int_t), 1, fFile);
        fwrite(&ly->fWScale, sizeof(Float_t), 1, fFile);
        fwrite(&ly->fXScale, sizeof(Float_t), 1, fFile);
        Int_t nw = ly->fW.GetSize();
        fwrite(&nw, sizeof(Int_t), 1, fFile);
        fwrite(ly->fW.fArray, sizeof(Short_t), nw, fFile);
        fwrite(ly->fB.fArray, sizeof(Int_t), ly->fOutNodes, fFile);
        fwrite(ly->fWFloat.fArray, sizeof(Double_t), nw, fFile);
        fwrite(ly->fBFloat.fArray, sizeof(Double_t), ly->fOutNodes, fFile);
    }
}

void TXMLPInt16::ReadBinary()
{
    fread(&fParm, sizeof(TNeuralNetParameters), 1, fFile);
    Int_t ot = kFermi, nL = 0;
    fread(&ot, sizeof(Int_t), 1, fFile);
    // old format stored relu_out 0/1 — values 0,1,2 are valid transfers
    fOutTransfer = ot;
    fReLUOutput = (ot == kReLU);
    fread(&fMomentum, sizeof(Double_t), 1, fFile);
    fread(&nL, sizeof(Int_t), 1, fFile);
    ClearLayers();
    fLayers.SetOwner(kTRUE);
    for (Int_t L = 0; L < nL; ++L) {
        Int_t inn = 0, outn = 0, tr = kReLU, nw = 0;
        Double_t ls = 0.01;
        Float_t ws = 1.f, xs = 1.f;
        fread(&inn, sizeof(Int_t), 1, fFile);
        fread(&outn, sizeof(Int_t), 1, fFile);
        fread(&ls, sizeof(Double_t), 1, fFile);
        fread(&tr, sizeof(Int_t), 1, fFile);
        fread(&ws, sizeof(Float_t), 1, fFile);
        fread(&xs, sizeof(Float_t), 1, fFile);
        fread(&nw, sizeof(Int_t), 1, fFile);
        // old files stored relu as 0/1 in tr slot — map 0→linear/fermi, 1→relu
        if (tr == 0 && L == nL - 1) tr = fOutTransfer;
        auto* ly = new TXMLPInt16Layer(inn, outn, ls, tr);
        ly->fWScale = ws;
        ly->fXScale = xs;
        if (ly->fW.GetSize() != nw) ly->fW.Set(nw);
        if (ly->fWFloat.GetSize() != nw) ly->fWFloat.Set(nw);
        fread(ly->fW.fArray, sizeof(Short_t), nw, fFile);
        fread(ly->fB.fArray, sizeof(Int_t), outn, fFile);
        fread(ly->fWFloat.fArray, sizeof(Double_t), nw, fFile);
        fread(ly->fBFloat.fArray, sizeof(Double_t), outn, fFile);
        fLayers.Add(ly);
    }
    if (fOut) {
        delete[] fOut;
        fOut = nullptr;
    }
    fOut = new Double_t[fParm.fOutNodes];
    fQuantized = kTRUE;
    EnsureScratch();
}
