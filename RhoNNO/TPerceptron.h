#ifndef TPERCEPTRON_H
#define TPERCEPTRON_H
// TPerceptron
//
// Implementation of the perceptron (Supervised Learning)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "VSupervisedNet.h"

class PerceptronUnit {
public:   
    double*  fVector;    // synaptic vector
    double*  fDelta;	   // modification vector
    double   fThreshold; // activity threshold
    int      fID;        // ID of this element
};

class PerceptronBase {
public:
};

class TPerceptron : public VSupervisedNet {
public:
    TPerceptron(int fInNodes,       // constructor for first perceptron
	int fOutNodes,
	double fLearnStep,
	TNeuralNetParameters::TRANSFER fTransferId,
	int fPerceptronId);
    
    TPerceptron(TPerceptron* Prev,    // constructor for linked perceptron
	int fOutNodes,
	double learnstep,
	TNeuralNetParameters::TRANSFER fTransferId,
	int fPerceptronId);
    
    TPerceptron(void);                // constructor for first perceptron (data from file)
    
    TPerceptron(TPerceptron* Prev);   // constructor for linked perceptron (data from file)
    
    virtual ~TPerceptron();           // destructor of network
    
    PerceptronUnit* fU;	//! Temp. unit
    void (*Transfer)(double in,double* out,double* deriv);	//! Transfer function
    
    double* fIn;   //!inputvector:  in[fInNodes]
    double* fOut;  //!outputvector: out[fOutNodes];
    double* fDiffSrc; //!error derivation from following net: fDiffSrc[fOutNodes]
    double* fDiffDst; //!error derivation to previous net:    fDiffDst[fInNodes]
    
private:
    TPerceptron* fPrev;    //! previous perceptron
    PerceptronUnit* fUbound;	//! Temp. unit
    void AllocNet(void);
    void InitNet(void);
    
public:
    void SetFile(FILE *f) { fFile = f; }
    virtual void WriteText();
    virtual void WriteBinary();
    virtual void ReadText();
    virtual void ReadBinary();
    
    double Train(NNO_INTYPE* in=0,NNO_OUTTYPE* out=0);  // calling Learnstep, Recallstep must have been performed already
    double* Recall(NNO_INTYPE* in=0,NNO_OUTTYPE* out=0);
    //void CopyData(const TPerceptron& PERC); // copies data from another perceptron
    
    ClassDef(TPerceptron,1)	// Multilayer Perceptron
};

#endif

