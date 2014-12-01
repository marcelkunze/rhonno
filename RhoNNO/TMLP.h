#ifndef TMLP_H
#define TMLP_H
//////////////////////////////////////////////////////////////////////////
//									//
// TMLP									//
//									//
// Implementation of the Multi-Layer-Perceptron (1 hidden Layer)	//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/TXMLP.h"

class TMLP : public TXMLP {
public:
    TMLP() {}
    TMLP(Double_t hidStep,Double_t outStep,Int_t innodes,Int_t hidnodes,Int_t outnodes,Double_t inputRange,
	const char* netFile,TNeuralNetParameters::TRANSFER f=TNeuralNetParameters::TR_FERMI)
	: TXMLP(2,inputRange,netFile,innodes,hidnodes,outnodes,hidStep,outStep,TNeuralNetParameters::TR_FERMI,f) {}
    TMLP(const char* netFile) : TXMLP(netFile) {ReadNet("XMLP");};
    virtual ~TMLP() {}                          //destructor of network  (File will be saved)
    
    ClassDef(TMLP,1)	    // Multi-Layer Perceptron
};

#endif
