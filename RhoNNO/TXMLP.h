#ifndef TXMLP_H
#define TXMLP_H
// TXMLP
//
// Implementation of the Multi-Layer-Perceptron
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "VSupervisedNet.h"
#include "TPerceptron.h"

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
    TXMLP(int layers,double inputRange,std::string netFile,int innodes,...);
    TXMLP(int layers,double inputRange,std::string netFile,int innodes,int n1,int n2,int n3,double s1,double s2,double s3,
          TNeuralNetParameters::TRANSFER f1,TNeuralNetParameters::TRANSFER f2,TNeuralNetParameters::TRANSFER f3);
    TXMLP(std::string netFile) : VSupervisedNet(netFile) {ReadNet("XMLP");};
    virtual ~TXMLP();
    double Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    TPerceptron* GetPerceptron(int i) { return fPerc[i]; }
    virtual void SetMomentumTerm(double f);
    
    TPerceptron** fPerc;    //! Temp.unit

    ClassDef(TXMLP,1) // Multi-Layer Perceptron (extended)
};

#endif

