#ifndef VNEURALNETPLOTTER_H
#define VNEURALNETPLOTTER_H
//////////////////////////////////////////////////////////////////////////
//									//
// VNeuralNetPlotter							//
//									//
// Base classes for network plotters					//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// M.Kunze, Bochum University						//
// (C) Copyright Marcel Kunze 2001, Ruhr-University Bochum.		//
//									//
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TString.h"

class TH1D;
class TGraph;
class TCanvas;

// Base class of all network plotterss
class VNeuralNetPlotter : public TNamed {
public:
    VNeuralNetPlotter(const char *name);
    virtual ~VNeuralNetPlotter() {}
    virtual void Initialize() = 0;
    virtual void Reset() = 0;
    virtual void AddTrainSample(Double_t trn, Bool_t good = kTRUE) = 0;
    virtual void AddTestSample(Double_t tst, Bool_t good = kTRUE) = 0;
    virtual void AddTrainGraph(Double_t trn) = 0;
    virtual void AddTestGraph(Double_t tst) = 0;
    virtual void ShowPlots() = 0;
    virtual void DrawT(Text_t *text, Float_t x, Float_t y, Float_t angle=0., Int_t color=1);

    ClassDef(VNeuralNetPlotter,1)  // Base class of all network plotters
};

// Default plotter
class TSimpleNeuralNetPlotter : public VNeuralNetPlotter {
protected:
    Bool_t	    fPlots;	    //!Show plots
    TString	    fTrnPlot;	    //!
    TString	    fTstPlot;	    //!
    Int_t	    fIndex;	    //!the index number
    TH1D*	    fTrnHistTrue;   //!Training histogram
    TH1D*	    fTrnHistFalse;  //!Training histogram
    TH1D*	    fTstHistTrue;   //!Test histogram
    TH1D*	    fTstHistFalse;  //!Test histogram
    TGraph*	    fTrnCurve;	    //!Error plot
    Int_t	    fNtrn;	    //!Number of points
    Double_t*	    fXtrn;	    //!
    Double_t*	    fYtrn;	    //!
    TGraph*	    fTstCurve;	    //!Test plot
    Int_t	    fNtst;	    //!Number of points
    Double_t*	    fXtst;	    //!
    Double_t*	    fYtst;	    //!
    TCanvas*	    fCanvas;	    //!Drawing

public:
    TSimpleNeuralNetPlotter(const char *name);
    virtual ~TSimpleNeuralNetPlotter();
    void Initialize();
    void Reset();
    void AddTrainSample(Double_t trn, Bool_t good = kTRUE);
    void AddTestSample(Double_t tst, Bool_t good = kTRUE);
    void AddTrainGraph(Double_t trn);
    void AddTestGraph(Double_t tst);
    void ShowPlots();

    ClassDef(TSimpleNeuralNetPlotter,1)  // Base class of all network plotters
};
    
#endif
