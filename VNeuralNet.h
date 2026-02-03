#ifndef VNEURALNET_H
#define VNEURALNET_H
// VNeuralNet
//
// Base classes for unsupervised and supervised networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TNamed.h"
#include <string>

// Define the precision
typedef float NNO_INTYPE;
typedef float NNO_OUTTYPE;


class VNeuralNetPlotter;
class TDataServe;
class TTree;
class TGraph;
class TCanvas;

// Base struct of all networks
class TNeuralNetParameters : public TObject {
public:
    // ID of transferfunction
    enum TRANSFER {TR_USER,TR_FERMI,TR_LINEAR,TR_LINEAR_BEND,TR_SIGMOID,TR_RELU};
    
    char fNetId[9];
    int    fLayers;	    // number of perceptron layers
    double fInScale;	    // scale input vector
    int    fInNodes;	    // number of input nodes
    int    fOutNodes;	    // number of output nodes
    double fLearnStep;    // learning step
    double fMu;	    // momentum term
    double fFse;	    // flat spot elimination
    TRANSFER fTransferId;   // transfer function
    int fPerceptronId;    // ID of perceptron
    double fThreshold;    // Threshold for output
public:
    TNeuralNetParameters();
    virtual ~TNeuralNetParameters() {}

    ClassDef(TNeuralNetParameters,1)  // Parameters for all supervised networks#endif
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
    virtual double* Inference(NNO_INTYPE* in,NNO_OUTTYPE* out=0) = 0;
    virtual double Train(NNO_INTYPE* in,NNO_OUTTYPE* out=0) = 0;    // returns squared error
    
    // For backwards compatibility
    double*	Inferencestep(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { return Inference(in,out); }
    double	Learnstep(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { return Train(in,out); }
    
    // Training and testing
    double	TrainEpoch(TDataServe *server, int nEpoch=1);
    double	TestEpoch(TDataServe *server);
    double	TrainEpoch(std::string file, int nEpoch=1);
    double	TestEpoch(std::string file);
    double	Test(NNO_INTYPE* in,NNO_OUTTYPE* trn);        // returns squared error
    void	BalanceSamples(bool yesNo = true) { fBalance = yesNo; }
    virtual void SetMomentumTerm(double f);
    virtual void SetFlatSpotElimination(double f);
    
    enum FILE_TYPE {FILE_BINARY,FILE_TEXT,FILE_ROOT};
    
private:
    bool  fBalance;		    //!Balance positive and negative samples
    bool  fOwnPlotter;	    //!Ownership of plotter
    
    void    WriteNetText();
    void    WriteNetBinary();
    void    ReadNetText();
    void    ReadNetBinary();
    void    WriteNet();
    unsigned int  BalancedTrnIndex(TDataServe *server);
    unsigned int  BalancedTstIndex(TDataServe *server);
    
protected:
    std::string	    fFilename;	    // Name of network file
    FILE_TYPE	    fFiletype;	    // Type of file
    TNeuralNetParameters  fParm;    // Topology of network
    VNeuralNetPlotter* fPlotter;    //!Show plots
    double*	    fOut;	    //!Output nodes
    bool	    fShouldSave;    //!Save the network
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
    VNeuralNet(std::string netID,int innodes,int outnodes,std::string netFile);
    VNeuralNet(std::string netFile);
    virtual ~VNeuralNet();
    void Save(std::string file);
    void Save();
    
    // Getter functions
    TNeuralNetParameters& GetParameters() { return fParm; }
    VNeuralNetPlotter& GetPlotter() const { return *fPlotter; }
    std::string     GetFilename() const { return fFilename; }
    double*	    GetOutput() const { return fOut; }
    double	    GetThreshold() const { return fParm.fThreshold; }
    std::string     GetNetID() const { return fParm.fNetId; }
    static double Random(void);
    
    // Setter functions
    void SetThreshold(double t) { fParm.fThreshold = t; }
    void SetPlotter(VNeuralNetPlotter *plotter) { fPlotter = plotter; }
    
    // Plotter functions
    void SetupPlots(VNeuralNetPlotter *plotter=0);
    void FillPlots(double trn=0.0, double tst=0.0);
    void ShowPlots();

    ClassDef(VNeuralNet,1)  // Base class of all networks
};

#endif
