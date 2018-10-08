#ifndef TLVQ_H
#define TLVQ_H
// TLVQ
//
// Implementation of the LEARNING-VECTOR-QUANTISATION (LVQ)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "VUnsupervisedNet.h"

class TLVQ : public VUnsupervisedNet {
public:
    TLVQ() {}
    TLVQ(int innodes,int cells,double winStep,std::string netFile);
    TLVQ(std::string netFile) : VUnsupervisedNet(netFile) {ReadNet("LVQ");};
    TLVQ(const TLVQ& lvq,std::string netFile); // copy constructor
    
    virtual ~TLVQ();                          //destructor of network  (File will be saved)
private:
    TNeuralNetCell*  fUbound;	//! Temp. unit
    TNeuralNetCell*  fUwin;	//! Temp. unit
    
    virtual void AllocNet(void);
    virtual void InitNet(void);
    virtual void WriteText(void);
    virtual void WriteBinary(void);
    virtual void ReadText(void);
    virtual void ReadBinary(void);
    
public:
    virtual double  Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    virtual double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* =0) { GetWinnerCell(in); return fUwin->fVector; }
    virtual int GetWinnerCell(NNO_INTYPE* in);
    
    void Implant(int nr,int Class,NNO_INTYPE* in);
    void CopyData(const TLVQ& lvq); // copies data from another lvq network
    
    ClassDef(TLVQ,1)	// Learning Vector Quantization
};

#endif
