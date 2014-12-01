#ifndef VNEURALNET_H
#define VNEURALNET_H
//////////////////////////////////////////////////////////////////////////
//									//
// VNeuralNet								//
//									//
// Base classes for unsupervised and supervised networks		//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TString.h"

// Define the precision
typedef Float_t NNO_INTYPE;
typedef Float_t NNO_OUTTYPE;


class VNeuralNetPlotter;
class TDataServe;
class TTree;
class TGraph;
class TCanvas;

// Base struct of all networks
class TNeuralNetParameters : public TObject {
public:
    // ID of transferfunction
    enum TRANSFER {TR_USER,TR_FERMI,TR_LINEAR,TR_LINEAR_BEND,TR_SIGMOID};

    char fNetId[9];
    Int_t    fLayers;	    // number of perceptron layers
    Double_t fInScale;	    // scale input vector
    Int_t    fInNodes;	    // number of input nodes
    Int_t    fOutNodes;	    // number of output nodes
    Double_t fLearnStep;    // learning step
    Double_t fMu;	    // momentum term
    Double_t fFse;	    // flat spot elimination
    TRANSFER fTransferId;   // transfer function
    Int_t fPerceptronId;    // ID of perceptron
    Double_t fThreshold;    // Threshold for output
public:
    TNeuralNetParameters();
    virtual ~TNeuralNetParameters() {}

    ClassDef(TNeuralNetParameters,1)  // Parameters for all supervised networks
};

// Base class of all networks
class VNeuralNet : public TNamed {
public:
    // Abstract interface for all networks

    virtual void AllocNet() = 0;
    virtual void InitNet() = 0;
    virtual void WriteText() = 0;
    virtual void WriteBinary() = 0;
    virtual void ReadText() = 0;
    virtual void ReadBinary() = 0;
    virtual Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out=0) = 0;
    virtual Double_t Train(NNO_INTYPE* in,NNO_OUTTYPE* out=0) = 0;    // returns squared error

    // For backwards compatibility
    Double_t*	Recallstep(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { return Recall(in,out); }
    Double_t	Learnstep(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { return Train(in,out); }

    // Training and testing
    Double_t	TrainEpoch(TDataServe *server, Int_t nEpoch=1);
    Double_t	TestEpoch(TDataServe *server);
    Double_t	TrainEpoch(const char *file, Int_t nEpoch=1);
    Double_t	TestEpoch(const char *file);
    Double_t	Test(NNO_INTYPE* in,NNO_OUTTYPE* trn);        // returns squared error
    void	BalanceSamples(Bool_t yesNo = kTRUE) { fBalance = yesNo; }
    virtual void SetMomentumTerm(Double_t f);
    virtual void SetFlatSpotElimination(Double_t f);

    enum FILE_TYPE {FILE_BINARY,FILE_TEXT,FILE_ROOT};

private:
    Bool_t  fBalance;		    //!Balance positive and negative samples
    Bool_t  fOwnPlotter;	    //!Ownership of plotter

    void    WriteNetText();
    void    WriteNetBinary();
    void    ReadNetText();
    void    ReadNetBinary();
    void    WriteNet();
    UInt_t  BalancedTrnIndex(TDataServe *server);
    UInt_t  BalancedTstIndex(TDataServe *server);
    
protected:
    TString	    fFilename;	    // Name of network file
    FILE_TYPE	    fFiletype;	    // Type of file
    TNeuralNetParameters  fParm;    // Topology of network
    VNeuralNetPlotter* fPlotter;    //!Show plots
    Double_t*	    fOut;	    //!Output nodes
    Bool_t	    fShouldSave;    //!Save the network
    FILE*	    fFile;	    //!File pointer

    // Help functions
    void freadvar(NNO_INTYPE Var) { fread(&Var,sizeof(Var),1,fFile); }
    void fwritevar(NNO_OUTTYPE Var) { fwrite(&Var,sizeof(Var),1,fFile); }
    void ReadNet(const char* netID);	
    void TestPointer(void* ptr);
    void Errorf(char* format,...);
    void Warningf(FILE* f,char* format,...);
    void Messagef(FILE* f,char* format,...);
    
public:
    VNeuralNet();
    VNeuralNet(const char* netID,Int_t innodes,Int_t outnodes,const char* netFile);
    VNeuralNet(const char* netFile);
    virtual ~VNeuralNet();
    void Save(char* file);
    void Save();

    // Getter functions
    TNeuralNetParameters& GetParameters() { return fParm; }
    VNeuralNetPlotter& GetPlotter() const { return *fPlotter; }
    const char *    GetFilename() const { return fFilename; }
    Double_t*	    GetOutput() const { return fOut; }
    Double_t	    GetThreshold() const { return fParm.fThreshold; }
    const char *    GetNetID() const { return fParm.fNetId; }
    static Double_t Random(void);

    // Setter functions
    void SetThreshold(Double_t t) { fParm.fThreshold = t; }
    void SetPlotter(VNeuralNetPlotter *plotter) { fPlotter = plotter; }

    // Plotter functions
    void SetupPlots(VNeuralNetPlotter *plotter=0);
    void FillPlots(Double_t trn=0.0, Double_t tst=0.0);
    void ShowPlots();

    ClassDef(VNeuralNet,1)  // Base class of all networks
};
    
#endif
