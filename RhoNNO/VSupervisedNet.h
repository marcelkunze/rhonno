#ifndef VSUPERVISEDNET_H
#define VSUPERVISEDNET_H
// VSupervisedNet
//
// Base classes for supervised learning
// Abstract base class of all unsupervised networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TNtuple.h"

#include "VNeuralNet.h"

class VSupervisedNet : public VNeuralNet {
protected:
    TNtuple *fTuple;  // Training data
public:
    VSupervisedNet() : VNeuralNet() {};
    VSupervisedNet(std::string netID,int innodes,int outnodes,std::string netFile) :
    VNeuralNet(netID,innodes,outnodes,netFile) {}
    VSupervisedNet(std::string netFile) :
    VNeuralNet(netFile) {}
    virtual ~VSupervisedNet() {};
    virtual long TrainEpoch (TNtuple *tuple, bool random=true);    // learn the hits from the ntuple
    
    ClassDef(VSupervisedNet,1) // Supervised training
};

#endif
