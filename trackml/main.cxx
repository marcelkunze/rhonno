// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TView.h"
#include "TPolyMarker3D.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"

#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <stack>
#include <queue>

#define MAXPARTICLES 5
#define MAXHITS 150000
#define DRAW true
#define EVALUATION true

std::string base_path = "/Users/marcel/workspace/train_sample/";

//Which event to run, this may be overwritten by main()'s arguments
int filenum = 21100;

//Not doing much in practice
int debug = 0;
//Not the standard assert!
#define assert(a, m) {if (!(a)) {cout << m << endl; exit(0);}}

void draw(long nhits,float *x,float *y,float *z,std::map<int,std::vector<Point> > tracks);

using namespace std;

int main(int argc, char**argv) {
    //Read which event to run from arguments, default is event # 1000
    //Supports events 0-124 for test, or 1000-1099 for validation (small dataset)

    if (argc >= 2) {
        filenum = atoi(argv[1]);
        cout << "Running on event #" << filenum << endl;
    }
    ios::sync_with_stdio(false);
    cout << fixed;

    if (EVALUATION) {
        Tracker::readBlacklist(base_path,filenum);
        Tracker::readParticles(base_path,filenum);
        Tracker::readTruth(base_path,filenum);
        Tracker::sortTracks();
    }
    Tracker::readHits(base_path,filenum);

    long nParticles = Tracker::truth_tracks.size();
    if (nParticles>MAXPARTICLES) nParticles = MAXPARTICLES;
    cout << "Particles: " << nParticles << endl;
    
    long nhits = Tracker::hits.size();
    float x[nhits],y[nhits],z[nhits];
    int label[nhits],truth[nhits],layer[nhits];

    nhits = 0;
    int n = 0;
    for (auto &track : Tracker::truth_tracks) {
        if (n++ >= MAXPARTICLES) break;
        vector<int> t = track.second;
        for (auto &id : t) {
            auto it = Tracker::track_hits.find(id);
            if (it==Tracker::track_hits.end()) continue;
            point &hit = it->second;
            x[nhits] = hit.x * 0.001; // in m
            y[nhits] = hit.y * 0.001; // in m;
            z[nhits] = hit.z * 0.001; // in m;
            label[nhits] = 0;
            truth[nhits] = id;
            point geo = Tracker::meta[id];
            int vol = geo.x;
            int lay = geo.y;
            layer[nhits] = Tracker::getLayer(vol,lay);
            nhits++;
        }
    }

    if (nhits > MAXHITS) nhits = MAXHITS;
    cout << "Hits: " << nhits << endl;
    
    cout << "Find tracks..." << endl;
    long nt = Tracker::findTracks((int)nhits,x,y,z,layer,label,truth);

    // Show the results
#define MAXLABEL 100
    cout << "Labels: ";
    for (int i=0;i<nhits;i++) {
        if (i<MAXLABEL || i>nhits-MAXLABEL) cout << label[i] << " ";
        if (i == MAXLABEL) cout << endl << "..." << endl;
    }
    cout << endl;

    cout << "Write submission file..." << endl;
    ofstream out("mysubmission.csv");
    if( !out.is_open() ){
        cout<<"Can not open output file"<<endl;
        exit(0);
    }
    for (int ih=0; ih<nhits; ih++ ){
        out<<filenum<<","<<ih+1<<","<<label[ih]<<endl;
    }
    out.close();

    // Assemble tracks
    cout << "Assemble tracks..." << endl;
    map<int,vector<Point> > tracks;
    for(int i=0; i<nt; i++) {
        int track = i+1;
        vector<Point> t;
        for (int j=0;j<nhits;j++) {
            if (track != label[j]) continue;
            Point p(x[j],y[j],z[j],j,track);
            t.push_back(p); // Save the results
        }
        tracks[i] = t;
    }
    
#define MAXTRACK 10
    cout << endl << "Number of tracks: " << nt << endl;
    for (int i=0; i<nt; i++) {
        vector<Point> t = tracks[i];
        if (i == MAXTRACK) cout << endl << "..." << endl;
        if (i<MAXTRACK || i>nt-MAXTRACK) {
            cout << "Track " << i+1 << ": ";
            for (auto it : t) {
                cout << it.id() << " ";
            }
            cout << endl;
        }
    }

    // Show the results in a canvas
    draw(nhits,x,y,z,tracks);
}

void draw(long nhits,float *x,float *y,float *z,map<int,vector<Point> > tracks)
{
    // Initialize a 3D canvas and draw the hits and tracks
    if (DRAW) {
        TString dir = base_path;
        TString filePrefix;
        filePrefix.Form("%sevent%09d",dir.Data(),filenum);
        TString cname = filePrefix+"-canvas.root";
        TFile output(cname,"RECREATE");
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
        for (int i=0;i<nhits;i++) {
            hitmarker->SetPoint(i,x[i],y[i],z[i]);
        }
        // set marker size, color & style
        hitmarker->SetMarkerSize(1.0);
        hitmarker->SetMarkerColor(kCyan);
        hitmarker->SetMarkerStyle(kStar);
        hitmarker->Draw();
        
        // Draw the tracks
        long nt = tracks.size();
        for (int i=0;i<nt;i++) {
            //cout << endl << "Drawing track " << i+1 << ": ";
            vector<Point> track = tracks[i];
            int n = 0;
            TPolyLine3D *connector = new TPolyLine3D((int)track.size());
            for (auto &hit : track)    {
                connector->SetPoint(n++, hit.x(), hit.y(), hit.z());
            }
            connector->SetLineWidth(1);
            connector->SetLineColor(kRed);
            connector->Draw();
        }
        
        cout <<  "Writing " << cname << endl;
        c1->Write();
        output.Close();
    }

}
