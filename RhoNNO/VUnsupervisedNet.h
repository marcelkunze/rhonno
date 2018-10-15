#ifndef VUNSUPERVISEDNET_H
#define VUNSUPERVISEDNET_H
// VUnsupervisedNet
//
// Base classes for unsupervised learning
// Abstract base class of all unsupervised networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TNtuple.h"

#include "VNeuralNet.h"
#include "TNeuralNetCell.h"

class VUnsupervisedNet : public VNeuralNet {
protected:
    TNeuralNetCell* fU;	    //! Temp. unit
    TNeuralNetCellParameters fXB;  //Network parameters
    TNtuple *fTuple;  // Training data
    
public:
    VUnsupervisedNet()  :  VNeuralNet(), fXB() {};
    VUnsupervisedNet(std::string netID,int innodes,int outnodes,std::string netFile) :
    VNeuralNet(netID,innodes,outnodes,netFile), fXB() {}
    VUnsupervisedNet(std::string netFile) :
    VNeuralNet(netFile), fXB() {}
    TNeuralNetCellParameters &GetParameters() { return fXB; }
    int GetNumberOfCells() const { return fXB.fCells; }
    const TNeuralNetCell* GetCell(int n) const { if (n<fXB.fCells) return &fU[n]; else return NULL;};
    virtual int GetWinnerCell(NNO_INTYPE* in) = 0;
    virtual long TrainEpoch (FILE* trainFile);  // returns number of records in LearnFile
    virtual long TrainEpoch (TNtuple *tuple, bool random=true);    // learn the hits from the ntuple
    virtual void Draw (Option_t *option="");
    virtual void Print (Option_t *option="") const;
    
    ClassDef(VUnsupervisedNet,1) // Unsupervised training
};

#endif
