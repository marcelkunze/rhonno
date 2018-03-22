#ifndef TNEURALNETCELL_H
#define TNEURALNETCELL_H
// Routines for Connected-Cell-Networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "TObject.h"

class TNeuralNetParameters;

class TNeuralNetCellParameters : public TObject{
public:
    TNeuralNetCellParameters();
    virtual ~TNeuralNetCellParameters() {}
    
private:
    Double_t fWinStep;	    //learnstep of winner
    Double_t fNeiStep;	    //learnstep of neighbour
    Double_t fNeuStep;	    //learnstep of neuron
    Int_t    fCells;	    //current used cells
    Int_t    fConnectors;   //maximum number of neighbors
    Double_t fWinCount;	    //decrement of win_count: win_count*=(1-fWinCount)
    Int_t    fMinCells;	    //number of cells required for one simplex
    Int_t    fMaxCells;	    //maximum number of used cells
    Long_t   fInsertStep;   //Learnsteps between two Insertions (if 0: user controlled Insertion)
    Long_t   fDeleteStep;   //Learnsteps between two Removals   (if 0: user controlled Removal)
    Long_t   fInsertCount;
    Long_t   fDeleteCount;
    Double_t fEdgeCount;    //decrement of edge_count: edge_count*=(1-fEdgeCount)
    Double_t fErrCount;	    //decrement of err_count: err_count*=(1-fErrCount)
    Double_t fNeiCount;     //update of s_dev on neighbours: s_dev+=(s_dist-s_dev)*d_s_dev
    Double_t fMinCount;	    //threshold for edge-removal
    Double_t fMainWinCount;
    Double_t fMainErrCount;
    Double_t fMainEdgeCount;
public:
    ClassDef(TNeuralNetCellParameters,1)  // Parameters for unsupervised networks
    
    friend class TNeuralNetCell;
    friend class VUnsupervisedNet;
    friend class TLVQ;
    friend class TGCS;
    friend class TGNG;
    friend class TGNGTracker;
    friend class TSGCS;
    friend class TSGNG;
};

union connector {
    void* fPtr;   // ptr to another cell
    Int_t fID;
};

class TNeuralNetCell : public TObject {
protected:
    Double_t* fVector;	    // !vector
    Int_t         fID;	    // cell ID
    Double_t   fCount;	    // error or winning count
    
    /* help variables */
    Double_t* fDiff ;	    // !difference to input (help variable)
    Int_t    fState ;	    // cell state
    
    connector*   fC ;	    // !connections
    Int_t        fNc;	    // number of connections
    Double_t   fChi2;	    // squared deviation
    
    Double_t* fWeight;	    // !neural weight
    Double_t  fOut;	    // cell output
    Double_t* fAge;	    // !age of connector
    Int_t   fClass;	    // user_definable
    
public:
    TNeuralNetCell();
    virtual ~TNeuralNetCell() {};
    const Double_t *GetVector() const { return fVector; }
    const Int_t     GetID() const { return fID; }
    const Double_t  GetWinnerCount() const { return fCount; }
    const Double_t *GetDifference() const { return fDiff; }
    const connector *GetConnectors() const { return fC; }
    const TNeuralNetCell *GetConnectedCell(UInt_t c) const { if (c<fNc) return (TNeuralNetCell *) (fC[c].fPtr); else return NULL;}
    const Int_t     GetNumberOfConnections() const { return fNc; }
    const Double_t  GetError() const { return fChi2; }
    
    static void Disconnect(TNeuralNetCell* up1,TNeuralNetCell* up2);
    static void Connect   (TNeuralNetCell* up1,TNeuralNetCell* up2,TNeuralNetCellParameters* XB);
    static void ConnectNew(TNeuralNetCell* unew,TNeuralNetCell* Unit1,TNeuralNetCell* Unit2,TNeuralNetCellParameters* XB);
    static void InitVector(TNeuralNetCell* unew,TNeuralNetParameters* b);
    static void InitVector(TNeuralNetCell* unew,TNeuralNetCell* U1,TNeuralNetCell* U2,TNeuralNetParameters* b);
    static void InitWgt   (TNeuralNetCell* unew,TNeuralNetParameters* b);
    static void GetSDev   (TNeuralNetCell* unit,TNeuralNetParameters* b);
    static void InitSDev  (TNeuralNetCell* unew,TNeuralNetParameters* b);
    static void InitCount (TNeuralNetCell* unew);
    static void CheckConnections(TNeuralNetCell* unit);
    static void WriteUnitText  (FILE* f,TNeuralNetCell* unit,TNeuralNetParameters* b);
    static void ReadUnitText   (FILE* f,TNeuralNetCell* unit,TNeuralNetParameters* b);
    static void WriteUnitBinary(FILE* f,TNeuralNetCell* unit,TNeuralNetParameters* b);
    static void ReadUnitBinary (FILE* f,TNeuralNetCell* unit,TNeuralNetParameters* b);
    
    ClassDef(TNeuralNetCell,1)  // Cell for unsupervised networks
    
    friend class TGCS;
    friend class TSGCS;
    friend class TGNG;
    friend class TGNGTracker;
    friend class TSGNG;
    friend class TLVQ;
};

#endif
