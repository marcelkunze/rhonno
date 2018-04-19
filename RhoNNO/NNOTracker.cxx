// Read HITS spatial data and try a GNG network to segment the tracks...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TFile.h>
#include <TVector3.h>
#include "RhoNNO/TGNGTracker.h"

#include <random>
#include <iostream>
#include <fstream>
using namespace std;

#define MAXEPOCH 100
#define INSERTSTEP 100
#define DELETESTEP 500

#define NHITS 50
#define SIGMA 0.001

// The user member function processes one event

std::vector<TVector3> hits;

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

void GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error) {
    default_random_engine generator;
    double tau = 0.025;
    for (int i=0; i<np; i++,tau+=delta)
    {
        Float_t X,Y,Z;
        X = radius * ( sin(phi + (signum(radius)) * tau) - sin(phi));
        Y = radius * (-cos(phi + (signum(radius)) * tau) + cos(phi));
        Z = gamma * tau;
        if (error > 0.0) {
            normal_distribution<float> distribution0(X,error);
            X = distribution0(generator);
            normal_distribution<float> distribution1(Y,error);
            Y = distribution1(generator);
            normal_distribution<float> distribution2(Z,error);
            Z = distribution2(generator);
        }
        points.push_back(TVector3(X,Y,Z));
    }
}

int main(int argc, char* argv[]) {
    
    TFile output("NNOTracker.root","RECREATE");
    
    string filename("event");
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
            TVector3 point(X,Y,Z);
            hits.push_back(point);
            //cout << point.x() << "\t" << point.y() << "\t" << point.z() << "\t"<< point.d() << endl;
        }
    }
    else
    {
        // std::vector<TVector3>, int np, float delta tau, float radius, float phi, float gamma
        GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/1.0,0.5,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/2.0,1.0,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/3.0,1.5,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/3.0,-1.2,SIGMA);
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
    vector<TVector3>::iterator it;
    for(it = hits.begin(); it != hits.end(); it++)    {
        static int i = 0;
        TVector3 p=*it;
        hitmarker->SetPoint(i++,p.x(),p.y(),p.z());
        nt1.Fill(p.x(),p.y(),p.z());
    }
    // set marker size, color & style
    hitmarker->SetMarkerSize(1.0);
    hitmarker->SetMarkerColor(kCyan);
    hitmarker->SetMarkerStyle(kStar);
    hitmarker->Draw();
    
    cout << endl << "Training Track data with a GNG Network";
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

    // TBD: Analyze the network
    // Sort out the tracks by following the network connections and fill the corresponding track hits into containers
    
    //TBD: Fit the helix tracks from the hits in the containers
    
    cout << endl << "Writing..." << endl;
    c1->Write();
    
    return EXIT_SUCCESS;
}


