#ifndef TXMLP_H
#define TXMLP_H
//////////////////////////////////////////////////////////////////////////
//									//
// TXMLP								//
//									//
// Implementation of the Multi-Layer-Perceptron				//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VSupervisedNet.h"
#include "RhoNNO/TPerceptron.h"

class TXMLP : public VSupervisedNet {
private:
    void AllocNet(void);
    void InitNet(void) {};
    void WriteText(void);
    void WriteBinary(void);
    void ReadText(void);
    void ReadBinary(void);
public:
    TXMLP() {};
    TXMLP(Int_t layers,Double_t inputRange,const char* netFile,Int_t innodes,...);
    TXMLP(Int_t layers,Double_t inputRange,const char* netFile,Int_t innodes,Int_t n1,Int_t n2,Int_t n3,Double_t s1,Double_t s2,Double_t s3,
	  TNeuralNetParameters::TRANSFER f1,TNeuralNetParameters::TRANSFER f2,TNeuralNetParameters::TRANSFER f3);
    TXMLP(const char* netFile) : VSupervisedNet(netFile) {ReadNet("XMLP");};
    virtual ~TXMLP();
    Double_t Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    TPerceptron* GetPerceptron(Int_t i) { return fPerc[i]; }
    virtual void SetMomentumTerm(Double_t f);
    
    TPerceptron** fPerc;    //! Temp.unit

    ClassDef(TXMLP,1) // Multi-Layer Perceptron (extended)
};

#endif

