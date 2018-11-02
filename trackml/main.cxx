// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include "TNtuple.h"
#include "TRandom.h"
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

#define MAXPARTICLES 20
#define MAXHITS 150000
#define TRAINFILE false
#define DRAW true
#define EVALUATION true

const std::string base_path = "/Users/marcel/workspace/train_sample/";

//Which event to run, this may be overwritten by main()'s arguments
int filenum = 21100;

//Not doing much in practice
int debug = 0;
//Not the standard assert!
#define assert(a, m) {if (!(a)) {cout << m << endl; exit(0);}}

void makeTrain2();
void makeTrain3();
void draw(long nhits,float *x,float *y,float *z,std::map<int,std::vector<Point> > tracks);

using namespace std;

TNtuple *ntuple2,*ntuple3;
TRandom r;

int main(int argc, char**argv) {
    //Read which event to run from arguments, default is event # 1000
    //Supports events 0-124 for test, or 1000-1099 for validation (small dataset)
    
    if (argc >= 2) {
        filenum = atoi(argv[1]);
        cout << "Running on event #" << filenum << endl;
    }
    ios::sync_with_stdio(false);
    cout << fixed;
    
    Tracker::verbose(true);
    
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
    int start[nParticles+1],end[nParticles+1];
    start[0] = 0;
    end[0] = -1;
    for (auto &track : Tracker::truth_tracks) {
        if (n++ >= MAXPARTICLES) break;
        vector<int> t = track.second;
        start[n] = end[n-1]+1;
        end[n] = end[n-1] + (int)t.size();
        for (auto &id : t) {
            auto it = Tracker::track_hits.find(id);
            if (it==Tracker::track_hits.end()) continue;
            point &hit = it->second;
            x[nhits] = hit.x * 0.001; // in m
            y[nhits] = hit.y * 0.001; // in m;
            z[nhits] = hit.z * 0.001; // in m;
            label[nhits] = 0;
            truth[nhits] = n; // true track assignment
            point geo = Tracker::meta[id];
            int vol = geo.x;
            int lay = geo.y;
            layer[nhits] = Tracker::getLayer(vol,lay);
            nhits++;
        }
        if (n<100) cout << "Track " << n << ": " << start[n] << "-" << end[n] << endl;
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
    
    // Generate a training sample for hit pairs
    if (TRAINFILE) {
        TString filePrefix;
        filePrefix.Form("%sevent%09d",base_path.c_str(),filenum);
        TString fname = filePrefix+".root";
        auto f = TFile::Open(fname,"RECREATE");
        cout << "Writing training data to " << fname << endl;
        cout << "Generate training file..." << endl;
        ntuple2 = new TNtuple("tracks","training data","rz1:phi1:z1:rz2:phi2:z2:truth");
        makeTrain2();
        ntuple2->Write();
        delete ntuple2;
        f->Close();
    }
    
    // Show the results in a canvas
    if (DRAW) draw(nhits,x,y,z,tracks);
}

// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain2()
{
    
    const int n=10; // hit pair layer combinations
    pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
    
    long wright=1,wrong=1;
    for (int i = 0; i < n; i++) {
        int tube1 = start_list[i].first;
        for (auto &a : Tracker::tubePoints[tube1]) {
            int tube2 = start_list[i].second;
            for (auto &b : Tracker::tubePoints[tube2]) {
                if (a.truth() == b.truth()) {
                    ntuple2->Fill(a.rz(),a.phi(),a.z(),b.rz(),b.phi(),b.z(),1.0); //wright combination
                    wright++;
                }
                else {
                    if (r.Rndm()<wright/wrong) {
                        ntuple2->Fill(a.rz(),a.phi(),a.z(),b.rz(),b.phi(),b.z(),0.0); //wrong combination
                        wrong++;
                    }
                }
            }
        }
    }
    cout << "Wright: " << wright << " Wrong: " << wrong << endl;
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
