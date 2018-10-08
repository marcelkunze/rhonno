#ifndef TSGCS_H
#define TSGCS_H
// TSGCS
//
// Implementation of the SUPERVISED-GROWING-CELL-STRUCTURE (SGCS)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "TNeuralNetCell.h"
#include "VSupervisedNet.h"

class TSGCS : public VSupervisedNet {
public:
    TSGCS() {};
    TSGCS(int   innodes,
          int    cells,
          int    outnodes,
          int    maxCells,
          double winStep,
          double neiStep,
          double neuStep,
          double aErrCount,
          double bSDev,
          int    connectors,
          long   insertStep,
          long   deleteStep,
          std::string  netFile);
    
    TSGCS(std::string netFile) : VSupervisedNet(netFile), fXB() {ReadNet("SGCS");};
    TSGCS(const TSGCS& sgcs,std::string netFile); // copy constructor
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
    int GetNumberOfCells() const { return fXB.fCells; }
    double Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    
    void CopyData(const TSGCS& sgcs); // copies data from another sgcs network
    
    int Insert(void); // this function will be called automatically if insert_step>0
    int Prune(void);  // this function will be called automatically if delete_step>0
    
    ClassDef(TSGCS,1) // Supervised Growing Cell Structure
};


#endif

