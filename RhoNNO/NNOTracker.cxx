// Read HITS spatial data and try a GNG network to segment the tracks...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include "RhoNNO/TGNGTracker.h"
#include "RhoNNO/TRadon.h"

#include <iostream>
#include <fstream>
using namespace std;

#define MAXEPOCH 100
#define INSERTSTEP 100
#define DELETESTEP 500

#define NHITS 25
#define SIGMA 0.002
#define THRESHOLD 100.

// The user member function processes one event

std::vector<Point> hits;

int main(int argc, char* argv[]) {
    
    TFile output("NNOTracker.root","RECREATE");
    
    TRadon radon(SIGMA,THRESHOLD);
    
    TString filename("event");
    TNtuple nt1("Hits","NNO Tracking Data","x:y:z");
    
    if (argc > 1) filename = argv[1];
    
    ifstream infile(filename);
    if (infile) {
        cout << "Reading input file: " << filename << endl;
        double X,Y,Z;
        while (infile >> X >> Y >> Z) {
            X*=0.01; // transform to meter
            Y*=0.01;
            Z*=0.01;
            Point point(X,Y,Z);
            hits.push_back(point);
            //cout << point.x() << "\t" << point.y() << "\t" << point.z() << "\t"<< point.d() << endl;
        }
    }
    else
    {
        // std::vector<Point>, int np, float delta tau, float radius, float phi, float gamma
        radon.GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/1.0,0.5,SIGMA);
        radon.GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/2.0,1.0,SIGMA);
        radon.GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/3.0,1.5,SIGMA);
        radon.GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/3.0,-1.2,SIGMA);
    }
    
    // Sort the hits according to distance from origin
    
    cout << "Sorting hits..." << endl;
    reverse(hits.begin(),hits.end());
    
    // Initialize a 3D canvas and draw the hits
    TCanvas *c1 = new TCanvas("c1","NNO Tracking: Neural Gas",200,10,700,500);
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
    // draw hits as PolyMarker3D
    long nhits = hits.size();
    TPolyMarker3D *hitmarker = new TPolyMarker3D((UInt_t) nhits);
    vector<Point>::iterator it;
    for(it = hits.begin(); it != hits.end(); it++)    {
        static int i = 0;
        Point p=*it;
        hitmarker->SetPoint(i++,p.x(),p.y(),p.z());
        nt1.Fill(p.x(),p.y(),p.z());
    }
    // set marker size, color & style
    hitmarker->SetMarkerSize(1.0);
    hitmarker->SetMarkerColor(kCyan);
    hitmarker->SetMarkerStyle(kStar);
    hitmarker->Draw();
    
    cout << endl << "Training Delphi TPC data with a GNG Network";
    cout << endl << " Number of hits:" << nhits << endl;
    
    TGNGTracker net( 3,            // 3 input nodes
                    (int) nhits,   // 1 cell per hit
                    0.5,           // Learnstep of Winner-Cell
                    0.05,          // Learnstep of Neighbours
                    0.8,           // a_win_count
                    0.01,          // a_edge_count
                    0.1,           // min_count
                    3,             // connectors
                    INSERTSTEP,    // insert_step
                    DELETESTEP,    // delete_step
                    "gng.net"      // Network Filename is "gng.net"
                    );
    
    long n = 0;
    while (n++ < MAXEPOCH) {
        cout << endl << "Epoch: " << n << endl << "Cells:" << net.GetNumberOfCells() << endl;
        // Set the input data
        net.TrainEpoch(&nt1, kFALSE);
    }
    
    // Show the network
    net.Print();
    net.Draw();
    
    cout << endl << "Fuzzy Radon transform";
    cout << endl << " Number of hits:" << nhits << endl;
    
    // Perform a Radon transformation from hit space to track parameter space
    
    std::vector<RADON>&result = radon.Transform(hits);
    cout << endl << endl << "Number of track candidates:" << result.size() << endl;
    cout << "Drawing..." << endl;
    radon.Draw();
    

    // TBD: Analyze the network
    // Sort out the tracks by following the network connections and fill the corresponding track hits into containers
    
    //TBD: Fit the helix tracks from the hits in the containers
    
    cout << endl << "Writing..." << endl;
    c1->Write();
    
    return EXIT_SUCCESS;
}

