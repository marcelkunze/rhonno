#ifndef TXMLPINT16_H
#define TXMLPINT16_H
// TXMLPInt16
//
// Int16-quantized Multi-Layer Perceptron:
//   hidden: ReLU, output: FERMI (sigmoid) by default for classification.
// ROOT streamers: TObject layers + TArray{S,I,D}.
//
// Training & Recall() use float shadow weights (same as TXMLP behaviour).
// Quantize() packs int16; RecallInt16() runs the quantized graph.
//

#include "VSupervisedNet.h"
#include "TXMLP.h"

#include "TArrayD.h"
#include "TArrayI.h"
#include "TArrayS.h"
#include "TObject.h"
#include "TObjArray.h"

#include <string>

/// One MLP layer — ROOT-streamable
class TXMLPInt16Layer : public TObject {
public:
    Int_t    fInNodes;
    Int_t    fOutNodes;
    TArrayS  fW;           // Short_t weights [out*in]
    TArrayI  fB;           // Int_t bias (accumulator domain)
    TArrayD  fWFloat;      // float shadow weights
    TArrayD  fBFloat;      // float shadow biases
    TArrayD  fDW;          // momentum buffer weights
    TArrayD  fDB;          // momentum buffer biases
    Float_t  fWScale;
    Float_t  fXScale;
    Int_t    fTransfer;    // 0=linear, 1=ReLU, 2=FERMI
    Double_t fLearnStep;

    TXMLPInt16Layer();
    TXMLPInt16Layer(Int_t inNodes, Int_t outNodes, Double_t learnStep, Int_t transfer);
    TXMLPInt16Layer(const TXMLPInt16Layer& o);
    TXMLPInt16Layer& operator=(const TXMLPInt16Layer& o);
    virtual ~TXMLPInt16Layer() {}

    void Resize(Int_t inNodes, Int_t outNodes);
    void ClearWeights();

    ClassDef(TXMLPInt16Layer, 2) // Int16 MLP layer
};

class TXMLPInt16 : public VSupervisedNet {
public:
    enum ETransfer { kLinear = 0, kReLU = 1, kFermi = 2 };

    TXMLPInt16();
    TXMLPInt16(Int_t layers, Double_t inputRange, std::string netFile,
               Int_t innodes, Int_t n0, Int_t n1, Int_t n2,
               Double_t s0, Double_t s1, Double_t s2,
               Bool_t reluOutput = kFALSE);
    TXMLPInt16(Int_t layers, Double_t inputRange, std::string netFile,
               Int_t innodes,
               const Int_t* nodes, const Double_t* learnSteps,
               Bool_t reluOutput = kFALSE);
    TXMLPInt16(std::string netFile);
    explicit TXMLPInt16(TXMLP& src, std::string netFile = "xmlpint16.net");
    TXMLPInt16(const TXMLPInt16& o);
    TXMLPInt16& operator=(const TXMLPInt16& o);
    virtual ~TXMLPInt16();

    Double_t Train(NNO_INTYPE* in, NNO_OUTTYPE* out) override;
    /// Float forward (training / classError). Matches TXMLP Recall semantics.
    Double_t* Recall(NNO_INTYPE* in, NNO_OUTTYPE* out = 0) override;
    /// Explicit int16 quantized forward (after Quantize).
    Double_t* RecallInt16(NNO_INTYPE* in, NNO_OUTTYPE* out = 0);
    void SetMomentumTerm(Double_t f) override;

    /// Output transfer: default FERMI for 0/1 classification (threshold 0.5).
    void SetOutputTransfer(Int_t t);
    Int_t GetOutputTransfer() const { return fOutTransfer; }
    /// If true, Recall() uses int16 path (default false).
    void SetUseInt16Recall(Bool_t v) { fUseInt16Recall = v; }
    Bool_t GetUseInt16Recall() const { return fUseInt16Recall; }

    void Quantize();
    Bool_t IsQuantized() const { return fQuantized; }

    Int_t GetNLayers() const { return fLayers.GetEntriesFast(); }
    TXMLPInt16Layer* GetLayer(Int_t i);
    const TXMLPInt16Layer* GetLayer(Int_t i) const;

    Double_t CompareFloatVsInt16(NNO_INTYPE* in);

    Bool_t WriteToRootFile(const char* path, const char* name = "xmlpint16");
    static TXMLPInt16* ReadFromRootFile(const char* path, const char* name = "xmlpint16");

protected:
    void AllocNet() override;
    void InitNet() override;
    void WriteText() override;
    void WriteBinary() override;
    void ReadText() override;
    void ReadBinary() override;

private:
    TObjArray fLayers;
    Bool_t    fQuantized;
    Bool_t    fReLUOutput;       // if true, last layer ReLU; else FERMI (classification)
    Bool_t    fUseInt16Recall;   // Recall() dispatches to int16 when true
    Int_t     fOutTransfer;      // kLinear / kReLU / kFermi
    Double_t  fMomentum;

    TArrayD fAct;   //!
    TArrayD fDelta; //!
    TArrayS fActQ;  //!

    void ClearLayers();
    void BuildLayers(Int_t layers, Int_t innodes,
                     const Int_t* nodes, const Double_t* steps,
                     Bool_t reluOutput);
    void EnsureScratch();
    void ForwardFloat(const Double_t* inScaled, Double_t* out);
    void ForwardInt16(const NNO_INTYPE* in, Double_t* out);

    static Double_t Activate(Int_t transfer, Double_t x);
    static Double_t ActivateDeriv(Int_t transfer, Double_t y /*post-act*/);
    static Short_t QuantizeValue(Double_t v, Float_t scale);
    static Float_t SymmetricScale(const Double_t* data, Int_t n, Double_t eps = 1e-8);

    ClassDef(TXMLPInt16, 2) // Int16 MLP ReLU/FERMI (ROOT I/O)
};

#endif
