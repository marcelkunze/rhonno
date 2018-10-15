// Read HITS spatial data and try a XMLP network to segment the tracks...
// M.Kunze, Heidelberg university, 2018

#include "TROOT.h"
#include "TCanvas.h"
#include "TView.h"
#include "TPolyMarker3D.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TNtuple.h"

#include "Tracker.h"

#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

#define MAXHITS 150000
#define NHITS 10
#define SIGMA 0.0
#define DRAW true
#define VERBOSE true

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

bool sortFunc( const vector<int>& p1,
              const vector<int>& p2 ) {
    return p1.size() > p2.size();
}


void print(vector<int> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
    cout << endl;
}

inline bool readLine( std::ifstream &in, float f[], int n)
{
    char c;
    for( int i=0; i<n-1; i++ ) in >> f[i] >> c;
    in >> f[n-1];
    return in.good();
}

inline bool readLine( std::ifstream &in, double f[], int n)
{
    char c;
    for( int i=0; i<n-1; i++ ) in >> f[i] >> c;
    in >> f[n-1];
    return in.good();
}

void GenerateTrack(std::vector<Point> &points, int np, double delta, double radius, double phi, double gamma, double error) {
    default_random_engine generator;
    double tau = 0.1;
    static int n = 0;
    for (int i=0; i<np; i++,tau+=delta)
    {
        float X,Y,Z;
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
        Point p;
        p.x = X;
        p.y = Y;
        p.z = Z;
        p.id = n++;
        p.val = -1; // Preset group with no match
        p.r = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
        p.phi = atan2(p.y,p.x);
        p.theta = acos(p.z/p.r);
        p.distance = 0.0;
        points.push_back(p);
    }
}

std::vector<Point> hits;
std::vector<Point> tracks[150000];

int main(int argc, char* argv[]) {
    
    TFile output("MLPTracker.root","RECREATE");
    
    string filename("event");
    TNtuple nt1("Hits","NNO Tracking Data","x:y:z");
    
    if (argc > 1) {
        filename = argv[1];
        
        ifstream in(filename);
        if (in) {
            cout << "Reading input file: " << filename << endl;
        }
        else {
            cout<<"Event "<<filename<<" does not exist!!"<<endl;
            exit(0);
        }
#ifdef TRACKML
        int n = 0;
        char tmpLine[256];
        in.getline(tmpLine,256);
        cout<<tmpLine<<endl;
        while (hits.size()<MAXHITS) {
            double h[7]; //file line: id:x:y:z:volume:layer:module
            if( !readLine(in,h,7) ) break;
            if( h[0]-1 != hits.size() ){
                cout<<"Hit index is wrong: "<<h[0]<<endl;
                exit(0);
            }
            Point p;
            p.x = h[1] * 0.001; // convert mm to m
            p.y = h[2] * 0.001; // convert mm to m
            p.z = h[3] * 0.001; // convert mm to m
            p.id = n++;
            p.val = -1; // Preset group with no match
            p.r = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
            p.phi = atan2(p.y,p.x);
            p.theta = acos(p.z/p.r);
            p.distance = 0.0;
            hits.push_back(p);
            if (hits.size()%1000 == 0) cout << hits.size() << endl;
        }
#else
        int n = 0;
        double X,Y,Z;
        while (in >> X >> Y >> Z) {
            if (n >= MAXHITS) break;
            Point p;
            p.x = X * 0.01; // convert cm to m
            p.y = Y * 0.01; // convert cm to m
            p.z = Z * 0.01; // convert cm to m
            p.id = n++;
            p.val = -1; // Preset group with no match
            p.r = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
            p.phi = atan2(p.y,p.x);
            p.theta = acos(p.z/p.r);
            p.distance = 0.0;
            hits.push_back(p);
        }
#endif
        cout<<" loaded "<<hits.size()<<" hits "<<endl;
        in.close();
    }
    else
    {
        // std::vector<Point>, int np, float delta tau, float radius, float phi, float gamma
        GenerateTrack(hits,NHITS,0.025, 0.5,M_PI/2.0, 1.0,SIGMA); // 00
        GenerateTrack(hits,NHITS,0.025,-1.0,M_PI/4.0, 1.0,SIGMA); // 10
        GenerateTrack(hits,NHITS,0.025, 1.5,M_PI/3.0,-1.0,SIGMA); // 20
        GenerateTrack(hits,NHITS,0.025,-2.0,M_PI/3.0,-1.0,SIGMA); // 30
    }
    
    unsigned long nhits = hits.size();
    cout << endl << "Reconstruct Track data with a XMLP Network";
    cout << endl << " Number of hits:" << nhits << endl;
    
    float x[nhits],y[nhits],z[nhits];
    int labels[nhits];
    int nt;
    for (int i=0; i<nhits; i++) {
        Point hit = hits[i];
        x[i] = hit.x;
        y[i] = hit.y;
        z[i] = hit.z;
    }
    
    nt = Tracker::findTracks((int)nhits,x,y,z,labels);

    // Assemble the tracks from label information
    for(int i=0; i<nt; i++) {
        int track = i+1;
        for (int j=0;j<nhits;j++) {
            if (track != labels[j]) continue;
            hits[j].val = labels[j];
            tracks[i].push_back(hits[j]); // Save the results
        }
    }

#define MAXTRACK 10
    cout << endl << "Number of tracks:" << nt << endl;
    for(int i=0; i<nt; i++) {
        vector<Point> t = tracks[i];
        if (i == MAXTRACK) cout << endl << "..." << endl;
        if (i<MAXTRACK || i>nt-MAXTRACK) {
            cout << "Track " << i+1 << ": ";
            for (vector<Point>::iterator it = t.begin(); it != t.end(); ++it) {
                Point p = *it;
                cout << p.id << " ";
            }
            cout << endl;
        }
    }

#define MAXLABEL 100
    cout << "Labels: ";
    for (int i=0;i<nhits;i++) {
        if (i<MAXLABEL || i>nhits-MAXLABEL) cout << labels[i] << " ";
        if (i == MAXLABEL) cout << endl << "..." << endl;
    }
    
    //TBD: Fit the helix tracks from the hits in the containers
    
    // Initialize a 3D canvas and draw the hits and tracks
    if (DRAW) {
        TCanvas *c1 = new TCanvas("c1","NNO Tracking: XMLP",200,10,700,500);
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
        TPolyMarker3D *hitmarker = new TPolyMarker3D((unsigned int) nhits);
        for (vector<Point>::iterator it = hits.begin(); it != hits.end(); it++)    {
            static int i = 0;
            Point p=*it;
            hitmarker->SetPoint(i++,p.x,p.y,p.z);
            nt1.Fill(p.x,p.y,p.z);
        }
        // set marker size, color & style
        hitmarker->SetMarkerSize(1.0);
        hitmarker->SetMarkerColor(kCyan);
        hitmarker->SetMarkerStyle(kStar);
        hitmarker->Draw();
        
        // Draw the tracks
        for (int i=0;i<nt;i++) {
            //cout << endl << "Drawing track " << i+1 << ": ";
            vector<Point> track = tracks[i];
            int n = 0;
            TPolyLine3D *connector = new TPolyLine3D((int)track.size());
            for (vector<Point>::iterator it = track.begin(); it != track.end(); it++)    {
                Point hit = *it;
                connector->SetPoint(n++, hit.x, hit.y, hit.z);
            }
            connector->SetLineWidth(1);
            connector->SetLineColor(kRed);
            connector->Draw();
        }
        
        cout << endl << "Writing..." << endl;
        c1->Write();
    }
    
    return EXIT_SUCCESS;
}

