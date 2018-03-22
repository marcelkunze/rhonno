// VUnsupervisedNet
//
// Base classes for unsupervised learning
// Abstract base class of all unsupervised networks
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "RhoNNO/VUnsupervisedNet.h"
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TNtuple.h>

#include <iostream>
using namespace std;

ClassImp(VUnsupervisedNet)

Long_t  VUnsupervisedNet::TrainEpoch(FILE* file) 
{
    if (file==0) return -1;
    NNO_INTYPE* Buf = new NNO_INTYPE[fParm.fInNodes];
    TestPointer(Buf);
    Long_t records = 0;
    rewind(file);
    while(!feof(file)) {
        if (fread(Buf,sizeof(NNO_INTYPE),fParm.fInNodes,file)==fParm.fInNodes) {
            Train(Buf);
            ++records;
        }
    }
    delete[] Buf;
    return records;
}

Long_t VUnsupervisedNet::TrainEpoch(TNtuple *tuple, Bool_t rand) {
    fTuple = tuple;
    if (fTuple == 0) return 0;
    Long_t nhits = fTuple->GetEntries();
    for (int i=0;i<nhits;i++) {
        Long_t index = i;
        if (rand) index = random()%nhits;
        fTuple->GetEvent(index,1);
        Float_t *x=fTuple->GetArgs();
        Learnstep(x);
    }
    return nhits;
}

void VUnsupervisedNet::Draw (Option_t *option) {
    // Draw the network
    UInt_t numberCells = GetNumberOfCells();
    TPolyMarker3D *cellmarker = new TPolyMarker3D(numberCells);
    
    for (int i=0;i<numberCells;++i) {
        const TNeuralNetCell *c = GetCell(i);
        const Double_t *x1 = c->GetVector();
        // Draw the cell location X1
        int dimension = fParm.fInNodes;
        if (dimension >= 3) {
            cellmarker->SetPoint(i, x1[0], x1[1], x1[2]);
        }
        else {
            cellmarker->SetPoint(i, x1[0], x1[1], 0.0);
        }
        cellmarker->SetMarkerSize(0.5);
        cellmarker->SetMarkerColor(kBlue);
        cellmarker->SetMarkerStyle(kFullDotLarge);
        cellmarker->Draw();

        UInt_t numberConnections = c->GetNumberOfConnections();
        for (int j=0;j<numberConnections;j++) {
            const TNeuralNetCell *u = c->GetConnectedCell(j);
            const Double_t *x2 = u->GetVector();
            // Draw the cell connections (X1->X2)
            TPolyLine3D *connector = new TPolyLine3D(2);
            if (dimension >= 3) {
                connector->SetPoint(0, x1[0], x1[1], x1[2]);
                connector->SetPoint(1, x2[0], x2[1], x2[2]);
            }
            else {
                connector->SetPoint(0, x1[0], x1[1], 0.0);
                connector->SetPoint(1, x2[0], x2[1], 0.0);
            }
            // set attributes
            connector->SetLineWidth(1);
            connector->SetLineColor(kRed);
            connector->Draw();
        }
    }
    
}

void VUnsupervisedNet::Print (Option_t *option) const {
    // Print the network
    cout << endl << "Neural Network Object "<< fParm.fNetId << endl;
    int dimension = fParm.fInNodes;
    UInt_t numberCells = GetNumberOfCells();
    for (int i=0;i<numberCells;++i) {
        const TNeuralNetCell *c = GetCell(i);
        const Double_t *x1 = c->GetVector();
        cout << "Cell " << i << ": X1(";
        for (int d=0;d<dimension-1;d++) cout << x1[d] << ",";
        cout << x1[dimension-1] << ")" << endl;
        UInt_t numberConnections = c->GetNumberOfConnections();
        for (int j=0;j<numberConnections;j++) {
            const TNeuralNetCell *u = c->GetConnectedCell(j);
            const Double_t *x2 = u->GetVector();
            const Int_t id = u->GetID();
            cout << "\t -> " << id << ". X2(";
            for (int d=0;d<dimension-1;d++) cout << x2[d] << ",";
            cout << x2[dimension-1] << ")" << endl;
        }
    }
    
}

