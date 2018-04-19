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
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include <TNtuple.h>

#include "RhoNNO/VNeuralNet.h"

class VSupervisedNet : public VNeuralNet {
protected:
    TNtuple *fTuple;  // Training data
public:
    VSupervisedNet() : VNeuralNet() {};
    VSupervisedNet(std::string netID,Int_t innodes,Int_t outnodes,std::string netFile) :
    VNeuralNet(netID,innodes,outnodes,netFile) {}
    VSupervisedNet(std::string netFile) :
    VNeuralNet(netFile) {}
    virtual ~VSupervisedNet() {};
    virtual Long_t TrainEpoch (TNtuple *tuple, Bool_t random=kTRUE);    // learn the hits from the ntuple
    
    ClassDef(VSupervisedNet,1) // Supervised training
};

#endif
