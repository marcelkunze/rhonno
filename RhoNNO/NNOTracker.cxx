// Read HITS spatial data and try a GNG network to segment the tracks...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include "RhoNNO/TGNG.h"

#include <iostream>
using namespace std;

#define MAXEPOCH 500
#define DELETESTEP 500

// The user member function processes one event

TNtuple	    hits("Hits","NNO tracking data","x:y:z");
Long_t      nhits=0, n=0;

int main(int argc, char* argv[]) {
    TString filename("event");
    if (argc > 1) filename = argv[1];
    
    cout << "Reading input file: " << filename << endl;
    FILE* F=fopen(filename,"r");
    if (!F) {
        cerr << "File does not exist!" << endl;
        return EXIT_FAILURE;
    }
    
    UInt_t i,j;
    double x,y,z;
    
    while(!feof(F)) {
        fscanf(F,"%lf %lf %lf",&x,&y,&z);
        hits.Fill(x,y,z);
    }
    fclose(F);
    
    TFile output("NNOTracker.root","RECREATE");
    
    nhits = hits.GetEntries();
    
    // Initialize a 3D canvas and draw the hits
    TCanvas *c1 = new TCanvas("c1","HNNO Tracking: Neural Gas",200,10,700,500);
    // create a pad
    TPad *p1 = new TPad("p1","p1",0.05,0.02,0.95,0.82,46,3,1);
    p1->Draw();
    p1->cd();
    // creating a view
    TView *view = TView::CreateView(1);
    view->SetRange(-200,-200,-200,200,200,200);
    // create a first PolyMarker3D
    TPolyMarker3D *hitmarker = new TPolyMarker3D((UInt_t) nhits);
    for (i=0;i<nhits;i++) {
        hits.GetEvent(i,1);
        Float_t *x=hits.GetArgs();
        hitmarker->SetPoint(i, x[0], x[1], x[2]);
    }
    // set marker size, color & style
    hitmarker->SetMarkerSize(2.0);
    hitmarker->SetMarkerColor(4);
    hitmarker->SetMarkerStyle(2);
    hitmarker->Draw();
    
    cout << endl << "Training Delphi TPC data with a GNG Network";
    cout << endl << " Number of hits:" << nhits << endl;
    
    TGNG net( 3,       // 3 input nodes
             (int) nhits,    // 1 cell per hit
             0.5,       // Learnstep of Winner-Cell
             0.05,      // Learnstep of Neighbours
             0.8,       // a_win_count
             0.01,      // a_edge_count
             0.1,       // min_count
             3,         // connectors
             10,       // insert_step
             DELETESTEP,       // delete_step
             "gng.net"  // Network Filename is "gng.net"
             );
    
    while (n++ < MAXEPOCH) {
        cout << endl << endl << "Epoch: " << n << endl << "Cells:" << net.GetNumberOfCells();
        // Perform training with all hits
        for (i=0;i<nhits;i++) {
            hits.GetEvent(random()%nhits,1);
            Float_t *x=hits.GetArgs();
            net.Learnstep(x);
        }
    }
    
    // Draw the network
    UInt_t numberCells = net.GetNumberOfCells();
    TPolyMarker3D *cellmarker = new TPolyMarker3D(numberCells);
    
    for (i=0;i<numberCells;++i) {
        const TNeuralNetCell *c = net.GetCell(i);
        const Double_t *x1 = c->GetVector();
        printf("\n Cell %d: X1(%f,%f,%f) \n",i,x1[0],x1[1],x1[2]);
        // Draw the cell location X1
        cellmarker->SetPoint(i, x1[0], x1[1], x1[2]);
        UInt_t numberConnections = c->GetNumberOfConnections();
        for (j=0;j<numberConnections;j++) {
            const TNeuralNetCell *u = c->GetConnectedCell(j);
            const Double_t *x2 = u->GetVector();
            const Int_t id = u->GetID();
            printf("\t -> %d. X2(%f,%f,%f) \n",id,x2[0],x2[1],x2[2]);
            // Draw the cell connections (X1->X2)
            TPolyLine3D *connector = new TPolyLine3D(2);
            connector->SetPoint(0, x1[0], x1[1], x1[2]);
            connector->SetPoint(1, x2[0], x2[1], x2[2]);
            // set attributes
            connector->SetLineWidth(2);
            connector->SetLineColor(1);
            connector->Draw();
        }
    }
    // set marker size, color & style
    cellmarker->SetMarkerSize(1.0);
    cellmarker->SetMarkerColor(1);
    cellmarker->SetMarkerStyle(8);
    cellmarker->Draw();
    
    // TBD: Analyze the network
    // Sort out the tracks by following the network connections and fill the corresponding track hits into containers
    
    //TBD: Fit the helix tracks from the hits in the containers
    
    c1->Write();
    c1->Close();
    
    return EXIT_SUCCESS;
}

