#ifndef VNEURALNETPLOTTER_H
#define VNEURALNETPLOTTER_H
// VNeuralNetPlotter
// TSimpleNeuralNetPlotter
//
// Base classes for network plotters
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

class TH1D;
class TGraph;
class TCanvas;

// Base class of all network plotters
#include <TNamed.h>
class VNeuralNetPlotter : public TNamed {
public:
    VNeuralNetPlotter(std::string name);
    virtual ~VNeuralNetPlotter() {}
    virtual void Initialize() = 0;
    virtual void Reset() = 0;
    virtual void AddTrainSample(double trn, bool good = true) = 0;
    virtual void AddTestSample(double tst, bool good = true) = 0;
    virtual void AddTrainGraph(double trn) = 0;
    virtual void AddTestGraph(double tst) = 0;
    virtual void ShowPlots() = 0;
    virtual void DrawT(Text_t *text, float x, float y, float angle=0., int color=1);

    ClassDef(VNeuralNetPlotter,1)  // Base class of all network plotters
};

// Default plotter
class TSimpleNeuralNetPlotter : public VNeuralNetPlotter {
protected:
    bool	    fPlots;	    //!Show plots
    std::string fTrnPlot;	//!
    std::string fTstPlot;	//!
    int	    fIndex;	    //!the index number
    TCanvas*    fCanvas;    //!Drawing
    TH1D*	    fTrnHistTrue;   //!Training histogram
    TH1D*	    fTrnHistFalse;  //!Training histogram
    TH1D*	    fTstHistTrue;   //!Test histogram
    TH1D*	    fTstHistFalse;  //!Test histogram
    TGraph*	    fTrnCurve;	    //!Error plot
    int	    fNtrn;	    //!Number of points
    double*	fXtrn;	    //!
    double*	fYtrn;	    //!
    TGraph*	    fTstCurve;	//!Test plot
    int	    fNtst;	    //!Number of points
    double*	fXtst;	    //!
    double*	fYtst;	    //!

public:
    TSimpleNeuralNetPlotter(std::string name);
    virtual ~TSimpleNeuralNetPlotter();
    void Initialize();
    void Reset();
    void AddTrainSample(double trn, bool good = true);
    void AddTestSample(double tst, bool good = true);
    void AddTrainGraph(double trn);
    void AddTestGraph(double tst);
    void ShowPlots();

    ClassDef(TSimpleNeuralNetPlotter,1)  // Base class of all network plotters
};
    
#endif
