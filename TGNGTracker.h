#ifndef TGNGTracker_H
#define TGNGTracker_H
// TGNGTracker
//
// Implementation of the GROWING-NEURAL-GAS Trackfinder
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Heidelberg University, 1995

#include "VUnsupervisedNet.h"

class TNtuple;

class TGNGTracker : public VUnsupervisedNet {
public:
    TGNGTracker() {}
    TGNGTracker(int innodes,
                int maxCells,
                double winStep,
                double neiStep,
                double aWinCount,
                double aEdgeCount,
                double minCount,
                int  connectors,
                long insertStep,
                long deleteStep,
                const char* netFile);
    
    TGNGTracker(std::string netFile) : VUnsupervisedNet(netFile) {ReadNet("GNG");};
    TGNGTracker(const TGNGTracker& gng,const char* netFile); // copy constructor
    
    virtual ~TGNGTracker();                          //destructor of network  (File will be saved)
    
protected:
    TNeuralNetCell*  fUbound;	//! Temp. unit
    TNeuralNetCell*  fUwin1;	//! Temp. unit
    TNeuralNetCell*  fUwin2;	//! Temp. unit
    double fMinDistSquare1;
    double fMinDistSquare2;
    
    virtual void AllocNet(void);
    virtual void InitNet(void);
    virtual void WriteText(void);
    virtual void WriteBinary(void);
    virtual void ReadText(void);
    virtual void ReadBinary(void);
    
    int  CondDisconnect(TNeuralNetCell* up1,TNeuralNetCell* up2);
    void Connect(TNeuralNetCell* up1,TNeuralNetCell* up2);
    virtual void UpdateConnector(TNeuralNetCell* up1,TNeuralNetCell* up2);
    
public:
    virtual double  Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    virtual double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* =0) { GetWinnerCell(in); return fUwin1->fVector; }
    virtual int GetWinnerCell(NNO_INTYPE* in);
    
    void Deviation(void); // calculates s_dev of all cells
    void CopyData(const TGNGTracker& gng); // copies data from another gcs network
    
    int Insert(void);  // this function will be called automatically if insert_step>0
    void Prune(void);  // this function will be called automatically if delete_step>0
    
    void SetHits(TNtuple* tuple) { fTuple = tuple; } // sets the input data
    TNtuple* GetHits(void) const { return fTuple; } // gets the input data
    
    int GetCluster(void); // Find a cluster
    
    ClassDef(TGNGTracker,1)	// Unsupervised Growing Neural Gas
};


#endif

