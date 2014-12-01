#ifndef TGNG_H
#define TGNG_H
//////////////////////////////////////////////////////////////////////////
//									//
// TGNG									//
//									//
// Implementation of the GROWING-NEURAL-GAS (GNG)			//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VUnsupervisedNet.h"

class TGNG : public VUnsupervisedNet {
public:
    TGNG() {}
    TGNG(Int_t innodes,
	Int_t maxCells,
	Double_t winStep,
	Double_t neiStep,
	Double_t aWinCount,
	Double_t aEdgeCount,
	Double_t minCount,
	Int_t  connectors,
	Long_t insertStep,
	Long_t deleteStep,
	const char* netFile);
    
    TGNG(const char* netFile) : VUnsupervisedNet(netFile) {ReadNet("GNG");};
    TGNG(const TGNG& gng,const char* netFile); // copy constructor
    
    virtual ~TGNG();                          //destructor of network  (File will be saved)

protected:
    TNeuralNetCell*  fUbound;	//! Temp. unit
    TNeuralNetCell*  fUwin1;	//! Temp. unit
    TNeuralNetCell*  fUwin2;	//! Temp. unit
    Double_t fMinDistSquare1;
    Double_t fMinDistSquare2;
    
    virtual void AllocNet(void);
    virtual void InitNet(void);
    virtual void WriteText(void);
    virtual void WriteBinary(void);
    virtual void ReadText(void);
    virtual void ReadBinary(void);
    
    Int_t  CondDisconnect(TNeuralNetCell* up1,TNeuralNetCell* up2);
    void Connect(TNeuralNetCell* up1,TNeuralNetCell* up2);
    virtual void UpdateConnector(TNeuralNetCell* up1,TNeuralNetCell* up2);
    
public:
    virtual Double_t  Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    virtual Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out=0) { GetWinnerCell(in); return fUwin1->fVector; }
    virtual Int_t GetWinnerCell(NNO_INTYPE* in);

    void Deviation(void); // calculates s_dev of all cells
    void CopyData(const TGNG& gng); // copies data from another gcs network
    
    Int_t Insert(void);  // this function will be called automatically if insert_step>0
    void Prune(void);  // this function will be called automatically if delete_step>0
    
    ClassDef(TGNG,1)	// Unsupervised Growing Neural Gas
};


#endif

