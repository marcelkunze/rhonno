#ifndef TLVQ_H
#define TLVQ_H
//////////////////////////////////////////////////////////////////////////
//									//
// TLVQ									//
//									//
// Implementation of the LEARNING-VECTOR-QUANTISATION (LVQ)		//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VUnsupervisedNet.h"

class TLVQ : public VUnsupervisedNet {
public:
    TLVQ() {}
    TLVQ(Int_t innodes,Int_t cells,Double_t winStep,const char* netFile);
    TLVQ(const char* netFile) : VUnsupervisedNet(netFile) {ReadNet("LVQ");};
    TLVQ(const TLVQ& lvq,const char* netFile); // copy constructor
    
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
    virtual Double_t  Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    virtual Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { GetWinnerCell(in); return fUwin->fVector; }
    virtual Int_t GetWinnerCell(NNO_INTYPE* in);
        
    void Implant(Int_t nr,Int_t Class,NNO_INTYPE* in);
    void CopyData(const TLVQ& lvq); // copies data from another lvq network
    
    ClassDef(TLVQ,1)	// Learning Vector Quantization
};

#endif
