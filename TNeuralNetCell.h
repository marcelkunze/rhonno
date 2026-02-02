#ifndef TNEURALNETCELL_H
#define TNEURALNETCELL_H
// Routines for Connected-Cell-Networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TObject.h"

class TNeuralNetParameters;

class TNeuralNetCellParameters : public TObject{
public:
    TNeuralNetCellParameters();
    virtual ~TNeuralNetCellParameters() {}
    
private:
    double fWinStep;	    //learnstep of winner
    double fNeiStep;	    //learnstep of neighbour
    double fNeuStep;	    //learnstep of neuron
    int    fCells;	    //current used cells
    int    fConnectors;   //maximum number of neighbors
    double fWinCount;	    //decrement of win_count: win_count*=(1-fWinCount)
    int    fMinCells;	    //number of cells required for one simplex
    int    fMaxCells;	    //maximum number of used cells
    long   fInsertStep;   //Learnsteps between two Insertions (if 0: user controlled Insertion)
    long   fDeleteStep;   //Learnsteps between two Removals   (if 0: user controlled Removal)
    long   fInsertCount;
    long   fDeleteCount;
    double fEdgeCount;    //decrement of edge_count: edge_count*=(1-fEdgeCount)
    double fErrCount;	    //decrement of err_count: err_count*=(1-fErrCount)
    double fNeiCount;     //update of s_dev on neighbours: s_dev+=(s_dist-s_dev)*d_s_dev
    double fMinCount;	    //threshold for edge-removal
    double fMainWinCount;
    double fMainErrCount;
    double fMainEdgeCount;
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
    int fID;
};

class TNeuralNetCell : public TObject {
protected:
    double* fVector;	    // !vector
    int         fID;	    // cell ID
    double   fCount;	    // error or winning count
    
    /* help variables */
    double* fDiff ;	    // !difference to input (help variable)
    int    fState ;	    // cell state
    
    connector*   fC ;	    // !connections
    int        fNc;	    // number of connections
    double   fChi2;	    // squared deviation
    
    double* fWeight;	    // !neural weight
    double  fOut;	    // cell output
    double* fAge;	    // !age of connector
    int   fClass;	    // user_definable
    
public:
    TNeuralNetCell();
    virtual ~TNeuralNetCell(void) {};
    const double *GetVector(void) const { return fVector; }
    int     GetID(void) const { return fID; }
    double  GetWinnerCount(void) const { return fCount; }
    const double *GetDifference(void) const { return fDiff; }
    const connector *GetConnectors(void) const { return fC; }
    const TNeuralNetCell *GetConnectedCell(int c) const { if (c<fNc) return (TNeuralNetCell *) (fC[c].fPtr); else return NULL;}
    int     GetNumberOfConnections(void) const { return fNc; }
    double  GetError() const { return fChi2; }
    
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
