#ifndef TSGNG_H
#define TSGNG_H
// TSGNG
//
// Implementation of the SUPERVISED-GROWING-NEURAL-GAS (SGNG)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "TNeuralNetCell.h"
#include "VSupervisedNet.h"

class TSGNG : public VSupervisedNet {
public:
    TSGNG() {};
    TSGNG(int innodes,
	int outnodes,
	int maxCells,
	double winStep,
	double neiStep,
	double neuStep,
	double aErrCount,
	double aEdgeCount,
	double bSDev,
	double minCount,
	int  connectors,
	long insertStep,
	long deleteStep,
	std::string netFile);
    TSGNG(std::string netFile) : VSupervisedNet(netFile), fXB() {ReadNet("SGNG");};
    TSGNG(const TSGNG& sgng,std::string netFile); // copy constructor
    virtual ~TSGNG();                          //destructor of network  (File will be saved)
       
private:
    TNeuralNetCell* fU;		//! Temp. Cell
    TNeuralNetCellParameters fXB;//  Network parameters
    TNeuralNetCell*  fUbound;	//! Temp. Cell
    TNeuralNetCell*  fUwin1;	//! Temp. Cell
    TNeuralNetCell*  fUwin2;	//! Temp. Cell
    double fMinDistSquare1;
    double fMinDistSquare2;
    
    void AllocNet(void);
    void InitNet(void);
    void WriteText(void);
    void WriteBinary(void);
    void ReadText(void);
    void ReadBinary(void);
    
    int  CondDisconnect(TNeuralNetCell* Up1,TNeuralNetCell* Up2);
    void Connect(TNeuralNetCell* Up1,TNeuralNetCell* Up2);
    void UpdateConnector(TNeuralNetCell* Up1,TNeuralNetCell* Up2);
    
public:
    TNeuralNetCellParameters &GetParameters() { return fXB; }
    int GetNumberOfCells() const { return fXB.fCells; }
    double Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    void   CopyData(const TSGNG& SGNG); // copies data from another sgng network
    
    int Insert(void);  // this function will be called automatically if insert_step>0
    void Prune(void);  // this function will be called automatically if delete_step>0
    
    ClassDef(TSGNG,1)	// Supervised Growing Neural Gas
};

#endif
