// Read HITS spatial data and try a Fuzzy Radon transform...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TFile.h>
#include "RhoNNO/TRadon.h"

#include <geneva/Go2.hpp>

#include <iostream>
#include <fstream>
using namespace std;
using namespace Gem::Geneva;

#define NHITS 50
#define SIGMA 0.001
#define THRESHOLD 300000.

// The user member function processes one event

std::vector<TVector3> hits;

int main(int argc, char* argv[]) {
    TFile output("RadonOptimizer.root","RECREATE");
    
    TRadon radon(SIGMA,THRESHOLD);
    
    Go2 go(argc, argv, "./config/Go2.json");
    GFMinIndividualFactory gfi("./config/GFMinIndividual.json");
    GEvolutionaryAlgorithmFactory ea("./config/GEvolutionaryAlgorithm.json");
    std::shared_ptr<GEvolutionaryAlgorithm> ea_ptr = ea.get<GEvolutionaryAlgorithm>();
    go & ea_ptr;
    
    string filename("event");
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
        // std::vector<Point>, int np, float delta tau, float radius, float phi, float gamma
        radon.GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/1.0,0.5,SIGMA);
//        radon.GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/2.0,1.0,SIGMA);
//        radon.GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/3.0,1.5,SIGMA);
//        radon.GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/3.0,-1.2,SIGMA);
    }
    
    // Sort the hits according to distance from origin
    
    cout << "Sorting hits..." << endl;
    reverse(hits.begin(),hits.end());
    
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
    // draw hits as PolyMarker3D
    long nhits = hits.size();
    TPolyMarker3D *hitmarker = new TPolyMarker3D((unsigned int) nhits);
    vector<TVector3>::iterator it;
    for(it = hits.begin(); it != hits.end(); it++)    {
        static int i = 0;
        TVector3 p=*it;
        hitmarker->SetPoint(i++,p.x(),p.y(),p.z());
    }
    // set marker size, color & style
    hitmarker->SetMarkerSize(1.0);
    hitmarker->SetMarkerColor(kCyan);
    hitmarker->SetMarkerStyle(kStar);
    hitmarker->Draw();
    
    cout << endl << "Fuzzy Radon transform";
    cout << endl << " Number of hits:" << nhits << endl;
    
    // Perform a Radon transformation from hit space to track parameter space
    
    double kappa=1.0, phi=M_PI/1.0, gamma=0.5, sigma=SIGMA;
    
    RADON t;
    t.kappa = kappa;
    t.phi   = phi;
    t.gamma = gamma;
    t.sigma = sigma;

    long n = radon.Density(t,hits);
    cout << "Radon hits: " << n << " Density: " << t.density << endl;
    
    std::vector<RADON>&result = radon.Transform(hits);
    cout << endl << endl << "Number of track candidates:" << result.size() << endl;
    cout << "Drawing..." << endl;
    radon.Draw();
    
    cout << endl << "Writing..." << endl;
    c1->Write();
    
    cout << endl;
    
    return EXIT_SUCCESS;
}

