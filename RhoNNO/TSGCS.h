#ifndef TSGCS_H
#define TSGCS_H
//////////////////////////////////////////////////////////////////////////
//									//
// TSGCS								//
//									//
// Implementation of the SUPERVISED-GROWING-CELL-STRUCTURE (SGCS)	//
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

class TSGCS : public VSupervisedNet {
public:
    TSGCS() {};
    TSGCS(Int_t   innodes,
	Int_t    cells,
	Int_t    outnodes,
	Int_t    maxCells,
	Double_t winStep,
	Double_t neiStep,
	Double_t neuStep,
	Double_t aErrCount,
	Double_t bSDev,
	Int_t    connectors,
	long   insertStep,
	long   deleteStep,
	const char*  netFile);
    
    TSGCS(const char* netFile) : fXB(), VSupervisedNet(netFile) {ReadNet("SGCS");};
    TSGCS(const TSGCS& sgcs,const char* netFile); // copy constructor
    virtual ~TSGCS();                          //destructor of network  (File will be saved)
       
private:
    TNeuralNetCell* fU;		//! Temp.Cell
    TNeuralNetCellParameters fXB;//  Network parameters
    TNeuralNetCell*  fUbound;	//! Temporary cell
    TNeuralNetCell*  fUwin;	//! Temp. Cell
    
    void AllocNet(void);
    void InitNet(void);
    void WriteText(void);
    void WriteBinary(void);
    void ReadText(void);
    void ReadBinary(void);
    
    void Remove(TNeuralNetCell* Urem);
    
public:
    TNeuralNetCellParameters &GetParameters() { return fXB; }
    Int_t GetNumberOfCells() const { return fXB.fCells; }
    Double_t Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    
    void CopyData(const TSGCS& sgcs); // copies data from another sgcs network
    
    Int_t Insert(void); // this function will be called automatically if insert_step>0
    Int_t Prune(void);  // this function will be called automatically if delete_step>0
    
    ClassDef(TSGCS,1) // Supervised Growing Cell Structure
};


#endif

