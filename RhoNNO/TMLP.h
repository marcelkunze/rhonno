#ifndef TMLP_H
#define TMLP_H
// TMLP
//
// Implementation of the Multi-Layer-Perceptron (1 hidden Layer)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TXMLP.h"

class TMLP : public TXMLP {
public:
    TMLP() {}
    TMLP(double hidStep,double outStep,int innodes,int hidnodes,int outnodes,double inputRange,
         std::string netFile,TNeuralNetParameters::TRANSFER f=TNeuralNetParameters::TR_FERMI)
    : TXMLP(2,inputRange,netFile,innodes,hidnodes,outnodes,hidStep,outStep,TNeuralNetParameters::TR_FERMI,f) {}
    TMLP(std::string netFile) : TXMLP(netFile) {ReadNet("XMLP");};
    virtual ~TMLP() {}                          //destructor of network  (File will be saved)
    
    ClassDef(TMLP,1)	    // Multi-Layer Perceptron
};

#endif
