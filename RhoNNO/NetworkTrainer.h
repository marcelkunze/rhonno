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
const long adjust_steps	= 20000;  // Decay of learning rate
const double adaptThreshold	= 3500.;
const double adaptRate	= 0.99;


class NetworkTrainer : public TObject
{
private:
    TDataServe *fPidDataServer, // input data set
    *fTrainingServer;           // training data set
    int	fStartEpoch;	    // first epoch and
    int	fStopEpoch;	        // last epoch to be trained
    int	fVectorsEpoch;	    // number of vectors per epoch
    TFile      *fFile;		    // input file
    std::string	fNetworkFile;	// name of the newtork file
    std::string	fDataPath;	    // path to data directory
    std::string	fNetworkPath;	// path to network directory
    VNeuralNet *fNet;		    // pointer to actual network
    double	fMomentum;	    // momentum term
    int	fTrnMax,	        // number of training vectors
    fTstMax;	                // number of test vectors
    std::string	fModel;		    // name of network model
    int	fInNodes,	        // size of input layer
    fHid1Nodes,	                // size of first hidden layer
    fHid2Nodes,	                // size of second hidden layer
    fOutNodes,	                // size of output layer
    fCells;		                // size growing net
    std::string	fInBranch[NNODIMENSION],  // names of input branches
    fOutBranch[NNODIMENSION];   // names of output branches
    bool	fBalance;	        // take equal number of pro and con samples
    bool	fPlots;		        // produce plots
    std::string	fTree;		        // name of tree
    std::string	fInput;		    // name of input file
    std::string	fOutput;	    // name of output file
    double	fScale;		    // global input scale
    bool	fAutoScale;	        // determine scale
    double	fInMean[NNODIMENSION],	// mean of inputs
    fOutMean[NNODIMENSION];	    // mean of outputs
    double	fInScale[NNODIMENSION],	// scale of inputs
    fOutScale[NNODIMENSION];    // scale of outputs
    float	fInVector[NNODIMENSION],// input vector
    fOutVector[NNODIMENSION];   // output vector
    std::string	fCut;		    // cut formula
    TList	fAll;		        // list of input files (pro and con)
    TList	fPro;		        // list of input files (pro)
    TList	fCon;		        // list of input files (con)
    TNeuralNetParameters::TRANSFER fTransfer; // transfer function
    
public:
    NetworkTrainer() {}
    NetworkTrainer(std::string file,int start=1,int end=100);
    virtual ~NetworkTrainer();
    bool	ReadSteeringFile(std::string file);
    void	SetupDataServer(std::string file="input.root");
    void	SetupNetworks();
    void	WriteSourceCode(std::string file);
    double	Train();
    double	Test();
    void	PrintOn();
    void	SetDataPath(std::string path) { fDataPath = path;}
    void	SetNetworkPath(std::string path) { fNetworkPath = path;}
    std::string     Makename(int z,  std::string fNetworkPath, std::string name);
    void	AutoScale(bool yesNo=kTRUE) { fAutoScale = yesNo; }
    void	BalanceSamples(bool yesNo=kTRUE) { fBalance = yesNo; }
    void	ShowControlPlots(bool yesNo=kTRUE) { fPlots = yesNo; }
    std::string GetModel() { return fModel; }
};

#endif

