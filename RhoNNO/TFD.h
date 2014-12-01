#ifndef TFD_H
#define TFD_H
//////////////////////////////////////////////////////////////////////////
//									//
// TFD									//
//									//
// Implementation of the Fisher discriminant (no hidden Layer)		//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/TXMLP.h"

class TFD : public TXMLP {
public:
    TFD() {}
    TFD(Int_t innodes,Int_t outnodes,TNeuralNetParameters::TRANSFER f,const char* netFile)
	: TXMLP(1,1.0,netFile,innodes,outnodes,0.1,f) {}
    TFD(const char* netFile) : TXMLP(netFile) {ReadNet("XMLP");};
    virtual ~TFD() {}                          //destructor of network  (File will be saved)
    
    ClassDef(TFD,1)	    // Fisher discriminant
};

#endif
