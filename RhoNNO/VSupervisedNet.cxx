// VSupervisedNet
//
// Base classes for supervised learning
// Abstract base class of all unsupervised networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "RhoNNO/VSupervisedNet.h"
#include "RhoNNO/VNeuralNetPlotter.h"

#include <iostream>
using namespace std;

ClassImp(VSupervisedNet)

Long_t VSupervisedNet::TrainEpoch(TNtuple *tuple, Bool_t randomize) {
    fTuple = tuple;
    if (fTuple == 0) return 0;
    Long_t nhits = fTuple->GetEntries();
    for (int i=0;i<nhits;i++) {
        Long_t index = i;
        if (randomize) index = rand()%nhits;
        fTuple->GetEvent(index,1);
        Float_t *x=fTuple->GetArgs();
        Learnstep(x, &x[fParm.fInNodes]); // the first fInNodes columns hold input data, the following fOutNodes columns hold the output data
    }
    return nhits;
}
