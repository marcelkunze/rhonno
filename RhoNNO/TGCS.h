#ifndef TGCS_H
#define TGCS_H
// TGCS
//
// Implementation of the GROWING-CELL-STRUCTURE (GCS)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "RhoNNO/VUnsupervisedNet.h"

class TGCS : public VUnsupervisedNet {
public:
    TGCS() {}
    TGCS(int innodes,int cells,int maxCells,double winStep,double neiStep,
         double aWinCount,int connectors,
         long insertStep,long deleteStep,std::string netFile);
    TGCS(std::string netFile) : VUnsupervisedNet(netFile) {ReadNet("GCS");};
    TGCS(const TGCS& gcs,std::string netFile); // copy constructor
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
    virtual double  Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    virtual double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* =0) { GetWinnerCell(in); return fUwin->fVector; }
    virtual int GetWinnerCell(NNO_INTYPE* in);
    
    void Deviation(void); // calculates s_dev of all cells
    void CopyData(const TGCS& GCS); // copies data from another gcs network
    
    int Insert(void);  // this function will be called automatically if insert_step>0
    //int Insert(NNO_INTYPE* In);
    int Prune(void);  // this function will be called automatically if delete_step>0
    
    ClassDef(TGCS,1)	// Unsupervised Growing Cellstructure
};

#endif
