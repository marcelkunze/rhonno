// NetworkTrainer
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// M.Kunze, Bochum University
// (C) Copyright M.Kunze 1997-2001, Ruhr-University Bochum.

#ifndef _NETWORKTRAINER_H_
#define _NETWORKTRAINER_H_

#include "TObject.h"
#include "TMath.h"
#include "string.h"
#include "TList.h"

#include "RhoNNO/VNeuralNet.h"

class TFile;
class TDataServe;

// Boundary conditions
#define NNODIMENSION 100
#define THRESHOLD 0.5

// Parameters for adaptation of learning rate
const Long_t adjust_steps	= 20000;  // Decay of learning rate
const Double_t adaptThreshold	= 3500.;
const Double_t adaptRate	= 0.99;


class NetworkTrainer : public TObject
{
private:
    Int_t	fStartEpoch;	    // first epoch and
    Int_t	fStopEpoch;	        // last epoch to be trained
    Int_t	fVectorsEpoch;	    // number of vectors per epoch
    TFile      *fFile;		    // input file
    std::string	fNetworkFile;	// name of the newtork file
    std::string	fDataPath;	    // path to data directory
    std::string	fNetworkPath;	// path to network directory
    VNeuralNet *fNet;		    // pointer to actual network
    Double_t	fMomentum;	    // momentum term
    TDataServe *fPidDataServer,	// input data set
    *fTrainingServer;           // training data set
    Int_t	fTrnMax,	        // number of training vectors
    fTstMax;	                // number of test vectors
    std::string	fModel;		    // name of network model
    Int_t	fInNodes,	        // size of input layer
    fHid1Nodes,	                // size of first hidden layer
    fHid2Nodes,	                // size of second hidden layer
    fOutNodes,	                // size of output layer
    fCells;		                // size growing net
    std::string	fInBranch[NNODIMENSION],  // names of input branches
    fOutBranch[NNODIMENSION];   // names of output branches
    Bool_t	fBalance;	        // take equal number of pro and con samples
    Bool_t	fPlots;		        // produce plots
    std::string	fTree;		        // name of tree
    std::string	fInput;		    // name of input file
    std::string	fOutput;	    // name of output file
    Double_t	fScale;		    // global input scale
    Bool_t	fAutoScale;	        // determine scale
    Double_t	fInMean[NNODIMENSION],	// mean of inputs
    fOutMean[NNODIMENSION];	    // mean of outputs
    Double_t	fInScale[NNODIMENSION],	// scale of inputs
    fOutScale[NNODIMENSION];    // scale of outputs
    Float_t	fInVector[NNODIMENSION],// input vector
    fOutVector[NNODIMENSION];   // output vector
    std::string	fCut;		    // cut formula
    TList	fAll;		        // list of input files (pro and con)
    TList	fPro;		        // list of input files (pro)
    TList	fCon;		        // list of input files (con)
    TNeuralNetParameters::TRANSFER fTransfer; // transfer function
    
public:
    NetworkTrainer() {}
    NetworkTrainer(std::string file,Int_t start=1,Int_t end=100);
    virtual ~NetworkTrainer();
    Bool_t	ReadSteeringFile(std::string file);
    void	SetupDataServer(std::string file="input.root");
    void	SetupNetworks();
    void	WriteSourceCode(std::string file);
    Double_t	Train();
    Double_t	Test();
    void	PrintOn();
    void	SetDataPath(std::string path) { fDataPath = path;}
    void	SetNetworkPath(std::string path) { fNetworkPath = path;}
    std::string     Makename(Int_t z,  std::string fNetworkPath, std::string name);
    void	AutoScale(Bool_t yesNo=kTRUE) { fAutoScale = yesNo; }
    void	BalanceSamples(Bool_t yesNo=kTRUE) { fBalance = yesNo; }
    void	ShowControlPlots(Bool_t yesNo=kTRUE) { fPlots = yesNo; }
    std::string GetModel() { return fModel; }
};

#endif

