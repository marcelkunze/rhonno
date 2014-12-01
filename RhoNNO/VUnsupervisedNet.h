#ifndef VUNSUPERVISEDNET_H
#define VUNSUPERVISEDNET_H
//////////////////////////////////////////////////////////////////////////
//									//
// VUnsupervisedNet							//
//									//
// Base classes for unsupervised learning				//
// Abstract base class of all unsupervised networks			//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VNeuralNet.h"
#include "RhoNNO/TNeuralNetCell.h"

class VUnsupervisedNet : public VNeuralNet {
protected:
    TNeuralNetCell* fU;	    //! Temp. unit
    TNeuralNetCellParameters fXB;  //Network parameters

public:
    VUnsupervisedNet()  : fXB(), VNeuralNet() {};
    VUnsupervisedNet(const char* netID,Int_t innodes,Int_t outnodes,const char* netFile) : 
      fXB(), VNeuralNet(netID,innodes,outnodes,netFile) {}
    VUnsupervisedNet(const char* netFile) : 
      fXB(), VNeuralNet(netFile) {}
    TNeuralNetCellParameters &GetParameters() { return fXB; }
    Int_t GetNumberOfCells() const { return fXB.fCells; }
    virtual Int_t GetWinnerCell(NNO_INTYPE* in) = 0;
    virtual Long_t TrainEpoch (FILE* trainFile);  // returns number of records in LearnFile
    
    ClassDef(VUnsupervisedNet,1) // Unsupervised training
};

#endif
