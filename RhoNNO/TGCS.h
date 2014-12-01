#ifndef TGCS_H
#define TGCS_H
//////////////////////////////////////////////////////////////////////////
//									//
// TGCS									//
//									//
// Implementation of the GROWING-CELL-STRUCTURE (GCS)			//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VUnsupervisedNet.h"

class TGCS : public VUnsupervisedNet {
public:
    TGCS() {}
    TGCS(Int_t innodes,Int_t cells,Int_t maxCells,Double_t winStep,Double_t neiStep,
	Double_t aWinCount,Int_t connectors,
	Long_t insertStep,Long_t deleteStep,const char* netFile);
    TGCS(const char* netFile) : VUnsupervisedNet(netFile) {ReadNet("GCS");};
    TGCS(const TGCS& gcs,const char* netFile); // copy constructor
    virtual ~TGCS();                             // destructor of network  (File will be saved)
private:
    TNeuralNetCell*  fUbound;	//! Temp. unit
    TNeuralNetCell*  fUwin;	//! Temp. unit
    
    virtual void AllocNet(void);
    virtual void InitNet(void);
    virtual void WriteText(void);
    virtual void WriteBinary(void);
    virtual void ReadText(void);
    virtual void ReadBinary(void);
    
    void Remove(TNeuralNetCell* Urem);
    
public:
    virtual Double_t  Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    virtual Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { GetWinnerCell(in); return fUwin->fVector; }
    virtual Int_t GetWinnerCell(NNO_INTYPE* in);
        
    void Deviation(void); // calculates s_dev of all cells
    void CopyData(const TGCS& GCS); // copies data from another gcs network
    
    Int_t Insert(void);  // this function will be called automatically if insert_step>0
    //Int_t Insert(NNO_INTYPE* In);
    Int_t Prune(void);  // this function will be called automatically if delete_step>0
    
    ClassDef(TGCS,1)	// Unsupervised Growing Cellstructure
};

#endif
