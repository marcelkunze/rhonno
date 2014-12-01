#ifndef TPERCEPTRON_H
#define TPERCEPTRON_H
//////////////////////////////////////////////////////////////////////////
//									//
// TPerceptron								//
//									//
// Implementation of the perceptron (Supervised Learning)		//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//	 
//////////////////////////////////////////////////////////////////////////	 

#include "RhoNNO/VSupervisedNet.h"

class PerceptronUnit {
public:   
    Double_t*  fVector;    // synaptic vector
    Double_t*  fDelta;	   // modification vector
    Double_t   fThreshold; // activity threshold
    Int_t      fID;        // ID of this element
};

class PerceptronBase {
public:
};

class TPerceptron : public VSupervisedNet {
public:
    TPerceptron(Int_t fInNodes,       // constructor for first perceptron
	Int_t fOutNodes,
	Double_t fLearnStep,
	TNeuralNetParameters::TRANSFER fTransferId,
	Int_t fPerceptronId);
    
    TPerceptron(TPerceptron* Prev,    // constructor for linked perceptron
	Int_t fOutNodes,
	Double_t learnstep,
	TNeuralNetParameters::TRANSFER fTransferId,
	Int_t fPerceptronId);
    
    TPerceptron(void);                // constructor for first perceptron (data from file)
    
    TPerceptron(TPerceptron* Prev);   // constructor for linked perceptron (data from file)
    
    virtual ~TPerceptron();           // destructor of network
    
    PerceptronUnit* fU;	//! Temp. unit
    void (*Transfer)(Double_t in,Double_t* out,Double_t* deriv);	//! Transfer function
    
    Double_t* fIn;   //!inputvector:  in[fInNodes]
    Double_t* fOut;  //!outputvector: out[fOutNodes];
    Double_t* fDiffSrc; //!error derivation from following net: fDiffSrc[fOutNodes]
    Double_t* fDiffDst; //!error derivation to previous net:    fDiffDst[fInNodes]
    
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
    
    Double_t Train(NNO_INTYPE* in=0,NNO_OUTTYPE* out=0);  // calling Learnstep, Recallstep must have been performed already
    Double_t* Recall(NNO_INTYPE* in=0,NNO_OUTTYPE* out=0);
    //void CopyData(const TPerceptron& PERC); // copies data from another perceptron
    
    ClassDef(TPerceptron,1)	// Multilayer Perceptron
};

#endif

