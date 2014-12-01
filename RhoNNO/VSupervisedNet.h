#ifndef VSUPERVISEDNET_H
#define VSUPERVISEDNET_H
//////////////////////////////////////////////////////////////////////////
//									//
// VSupervisedNet							//
//									//
// Base classes for Supervised Learning					//
// Abstract base class of all supervised networks			//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VNeuralNet.h"

class VSupervisedNet : public VNeuralNet {
public:
    VSupervisedNet() : VNeuralNet() {};
    VSupervisedNet(const char* netID,Int_t innodes,Int_t outnodes,const char* netFile) : 
      VNeuralNet(netID,innodes,outnodes,netFile) {}
    VSupervisedNet(const char* netFile) : 
      VNeuralNet(netFile) {}
    virtual ~VSupervisedNet() {};

    ClassDef(VSupervisedNet,1) // Supervised training
};

#endif
