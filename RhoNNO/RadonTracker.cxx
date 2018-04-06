// Read HITS spatial data and try a Fuzzy Radon transform...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include "RhoNNO/TRadon.h"

#include <iostream>
using namespace std;

// The user member function processes one event

TNtuple	    hits("Hits","Radon tracking data","x:y:z:sigma");
Long_t      nhits=0, n=0;

int main(int argc, char* argv[]) {
    TFile output("RadonTracker.root","RECREATE");
    
    TRadon radon;
    TString filename("event");
    if (argc > 1) filename = argv[1];
    
    FILE* F=fopen(filename,"r");
    if (F) {
        cout << "Reading input file: " << filename << endl;
        
        double X,Y,Z;
        
        while(!feof(F)) {
            fscanf(F,"%lf %lf %lf",&X,&Y,&Z);
            hits.Fill(X*0.01,Y*0.01,Z*0.01); // transform to meter
        }
        fclose(F);
    }
    else
    {
        // TNtuple *nt, int np, float delta tau, float radius, float phi, float gamma
        radon.GenerateTrack(&hits,50,0.0125,1.0,M_PI/1.0,0.5);
        radon.GenerateTrack(&hits,50,0.0125,1.0,M_PI/2.0,1.0);
        radon.GenerateTrack(&hits,50,0.0125,1.0,M_PI/3.0,1.5);
    }
    
    nhits = hits.GetEntries();
    
    // Initialize a 3D canvas and draw the hits
    TCanvas *c1 = new TCanvas("c1","Fuzzy Radon Tracking",200,10,700,500);
    // create a pad
    TPad *p1 = new TPad("p1","p1",0.05,0.02,0.95,0.82,46,3,1);
    p1->SetFillColor(kBlack);
    p1->Draw();
    p1->cd();
    // creating a view
    TView *view = TView::CreateView(1);
    view->SetRange(-2,-2,-2,2,2,2); // draw in a 2 meter cube
    // Draw axis
    TAxis3D rulers;
    rulers.Draw();
    // create a first PolyMarker3D
    TPolyMarker3D *hitmarker = new TPolyMarker3D((UInt_t) nhits);
    for (int i=0;i<nhits;i++) {
        hits.GetEvent(i,1);
        Float_t *x=hits.GetArgs();
        hitmarker->SetPoint(i, x[0], x[1], x[2]);
    }
    // set marker size, color & style
    hitmarker->SetMarkerSize(1.0);
    hitmarker->SetMarkerColor(kCyan);
    hitmarker->SetMarkerStyle(kStar);
    hitmarker->Draw();
    
    cout << endl << "Fuzzy Radon transform";
    cout << endl << " Number of hits:" << nhits << endl;
    
    // Perform a Radon transformation from hit space to track parameter space
    
    TNtuple *result = radon.Transform(&hits);
    cout << endl << endl << "Number of track candidates:" << result->GetEntries() << endl;
    cout << "Drawing..." << endl;
    radon.Draw();
    
    cout << endl << "Writing..." << endl;
    c1->Write();
    
    cout << endl;
    
    return EXIT_SUCCESS;
}

