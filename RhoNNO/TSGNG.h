#ifndef TSGNG_H
#define TSGNG_H
//////////////////////////////////////////////////////////////////////////
//									//
// TSGNG								//
//									//
// Implementation of the SUPERVISED-GROWING-NEURAL-GAS (SGNG)		//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/TNeuralNetCell.h"
#include "RhoNNO/VSupervisedNet.h"

class TSGNG : public VSupervisedNet {
public:
    TSGNG() {};
    TSGNG(Int_t innodes,
	Int_t outnodes,
	Int_t maxCells,
	Double_t winStep,
	Double_t neiStep,
	Double_t neuStep,
	Double_t aErrCount,
	Double_t aEdgeCount,
	Double_t bSDev,
	Double_t minCount,
	Int_t  connectors,
	Long_t insertStep,
	Long_t deleteStep,
	const char* netFile);
    TSGNG(const char* netFile) : fXB(), VSupervisedNet(netFile) {ReadNet("SGNG");};
    TSGNG(const TSGNG& sgng,const char* netFile); // copy constructor
    virtual ~TSGNG();                          //destructor of network  (File will be saved)
       
private:
    TNeuralNetCell* fU;		//! Temp. Cell
    TNeuralNetCellParameters fXB;//  Network parameters
    TNeuralNetCell*  fUbound;	//! Temp. Cell
    TNeuralNetCell*  fUwin1;	//! Temp. Cell
    TNeuralNetCell*  fUwin2;	//! Temp. Cell
    Double_t fMinDistSquare1;
    Double_t fMinDistSquare2;
    
    void AllocNet(void);
    void InitNet(void);
    void WriteText(void);
    void WriteBinary(void);
    void ReadText(void);
    void ReadBinary(void);
    
    Int_t  CondDisconnect(TNeuralNetCell* Up1,TNeuralNetCell* Up2);
    void Connect(TNeuralNetCell* Up1,TNeuralNetCell* Up2);
    void UpdateConnector(TNeuralNetCell* Up1,TNeuralNetCell* Up2);
    
public:
    TNeuralNetCellParameters &GetParameters() { return fXB; }
    Int_t GetNumberOfCells() const { return fXB.fCells; }
    Double_t Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    void   CopyData(const TSGNG& SGNG); // copies data from another sgng network
    
    Int_t Insert(void);  // this function will be called automatically if insert_step>0
    void Prune(void);  // this function will be called automatically if delete_step>0
    
    ClassDef(TSGNG,1)	// Supervised Growing Neural Gas
};

#endif
