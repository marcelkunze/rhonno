#ifndef VSUPERVISEDNET_H
#define VSUPERVISEDNET_H
//////////////////////////////////////////////////////////////////////////
//									                                    //
// VSupervisedNet							                            //
//									                                    //
// Base classes for Supervised Learning					                //
// Abstract base class of all supervised networks			            //
// Part of the Neural Network Objects package (NNO)			            //
//									                                    //
// Author List:								                            //
// Johannes Steffens, Bochum University					                //
// M.Kunze, Bochum University						                    //
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	    //
//									                                    //
//////////////////////////////////////////////////////////////////////////

#include <TNtuple.h>

#include "RhoNNO/VNeuralNet.h"

class VSupervisedNet : public VNeuralNet {
protected:
    TNtuple *fTuple;  // Training data
public:
    VSupervisedNet() : VNeuralNet() {};
    VSupervisedNet(const char* netID,Int_t innodes,Int_t outnodes,const char* netFile) :
    VNeuralNet(netID,innodes,outnodes,netFile) {}
    VSupervisedNet(const char* netFile) :
    VNeuralNet(netFile) {}
    virtual ~VSupervisedNet() {};
    virtual Long_t TrainEpoch (TNtuple *tuple, Bool_t random=kTRUE);    // learn the hits from the ntuple
    
    ClassDef(VSupervisedNet,1) // Supervised training
};

#endif
