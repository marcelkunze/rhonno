// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include "Graph.h"

#include "TNtuple.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TView.h"
#include "TPolyMarker3D.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#define MAXPARTICLES 10000
#define MAXHITS 150000
#define TRAINFILE true
#define DRAW true
#define EVALUATION true
#define VERBOSE false
#define MAXTRACK 10
#define MAXLABEL 100

const std::string base_path = "/Users/marcel/workspace/train_sample/";

//Which event to run, this may be overwritten by main()'s arguments
int filenum = 21100;

using namespace std;

long checkTracks(map<int,vector<int> >  &tracks);
void trainNetworks(string base_path,int filenum);
void draw(long nhits,float *x,float *y,float *z,map<int,vector<int> > tracks);

vector<pair<int,int> > truepairs;

int main(int argc, char**argv) {
    //Read which event to run from arguments, default is event # 1000
    //Supports events 0-124 for test, or 1000-1099 for validation (small dataset)
    
    if (argc >= 2) {
        filenum = atoi(argv[1]);
        cout << "Running on event #" << filenum << endl;
    }
    ios::sync_with_stdio(false);
    cout << fixed;
    
    Tracker::verbose(VERBOSE);
    
    if (EVALUATION) {
        //Tracker::readGraph("paths.csv",Tracker::paths);
        Tracker::readParticles(base_path,filenum);
        Tracker::readTruth(base_path,filenum);
        Tracker::sortTracks();
    }
    Tracker::readHits(base_path,filenum);
    Tracker::readCells(base_path,filenum);
    
    long nParticles = Tracker::truth_tracks.size();
    if (nParticles>MAXPARTICLES) nParticles = MAXPARTICLES;
    cout << "Particles: " << nParticles << endl;
    
    long nhits = Tracker::hits.size();
    float x[nhits],y[nhits],z[nhits],cx[nhits],cy[nhits],cz[nhits];
    int label[nhits],volume[nhits],layer[nhits],module[nhits],hitid[nhits];
    long long trackid[nhits];
    
    for (int i=0;i<nhits;i++) {
        if (Tracker::hit_dir[nhits][0].x != 0.0)
        cout << i << " " << Tracker::hit_dir[nhits][0].x << " " << Tracker::hit_dir[nhits][0].y << " " << Tracker::hit_dir[nhits][0].z << endl;
    }
    
    // Prepare the trackml data to run the track finder
    // Geberate a graph to represent the track hits in the modules
    
    nhits = 0;
    int n = 0;
    int start[nParticles+1],end[nParticles+1];
    start[0] = 0;
    end[0] = -1;
    for (auto &track : Tracker::particles) {
        vector<int> t = track.hit;
        if (t.size()==0) continue;
        if (n++ >= MAXPARTICLES) break;
        truepairs.push_back(make_pair(nhits,nhits+1));
        point geo = Tracker::meta[t[0]]; // Check the first layer of a hit
        int vol = geo.x;
        int lay = geo.y;
        int first = Tracker::getLayer(vol,lay);
        //if (first!=0 && first!=4 && first!=11) continue; // track does not start at first layers
        start[n] = end[n-1]+1;
        end[n] = end[n-1] + (int)t.size();
        //if (VERBOSE) cout << "Track  " << n << " {";
        int oldl = -1;
        int oldindex = -1;
        for (auto &id : t) {
            auto it = Tracker::track_hits.find(id);
            if (it==Tracker::track_hits.end()) continue;
            point &hit = it->second;
            x[nhits] = hit.x; // in mm
            y[nhits] = hit.y; // in mm
            z[nhits] = hit.z; // in mm
            label[nhits] = 0;
            trackid[nhits] = track.id; // true track assignment
            hitid[nhits] = id;
            point geo = Tracker::meta[id];
            int vol = geo.x;
            int lay = geo.y;
            int mod = geo.z;
            int l = Tracker::getLayer(vol,lay);
            volume[nhits] = vol;
            layer[nhits] = l;
            module[nhits] = mod;
            cx[nhits] = Tracker::hit_dir[nhits][0].x;
            cy[nhits] = Tracker::hit_dir[nhits][0].y;
            cz[nhits] = Tracker::hit_dir[nhits][0].z;
            //cout << nhits << " " << cx[nhits] << " " << cy[nhits] << " " << cz[nhits] << endl;
            int index = MODULES*l + mod;
            // Add the hit pair to the paths graph
            if (oldindex>-1 && oldindex!=index) {
                //Point p1(x[nhits-1],y[nhits-1],z[nhits-1]);
                //Point p2(x[nhits],y[nhits],z[nhits]);
                //double recall = Tracker::recall2(p1, p2)[0];
                //if (recall > THRESHOLD2)
                Tracker::paths.add(oldindex,index,1.0);
            }
            //if (VERBOSE) cout << "{" << index << ","  << l << "," << mod << "},";
            oldl = l;
            oldindex = index;
            nhits++;
        }
        Tracker::paths.add(oldindex,-1);
        if (VERBOSE) {
            //cout << "-1}" << endl;
        }
    }
    
    if (VERBOSE) {
        auto modpath = serialize(Tracker::paths);
        cout << "modpath:" << endl;
        for (auto &it : modpath) Tracker::print(it.second);
        for (int i=1;i<n;i++) cout << "Track " << i << ": " << start[i] << "-" << end[i] << endl;
        
        cout << truepairs.size() << " true pairs" << endl;
        for (auto p : truepairs) cout << "{" << p.first << "," << p.second << "}, ";
        cout << endl;
    }
    
    // Write path data to file
    Tracker::writeGraph("paths.csv",Tracker::paths);
    
    if (nhits > MAXHITS) nhits = MAXHITS;
    cout << "Hits: " << nhits << endl;
    
    cout << endl << "Running Tracker:" << endl;
    long nt = Tracker::findTracks((int)nhits,x,y,z,cx,cy,cz,volume,layer,module,hitid,trackid,label);
    
    // Show the results
    cout << "Labels: ";
    for (int i=0;i<nhits;i++) {
        if (i<MAXLABEL || i>nhits-MAXLABEL) cout << label[i] << " ";
        if (i == MAXLABEL) cout << endl << "..." << endl;
    }
    cout << endl;
    
    cout << "Assig.: ";
    for (int i=0;i<nhits;i++) {
        if (i<MAXLABEL || i>nhits-MAXLABEL) cout << Tracker::assignment[i] << " ";
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
    cout << "Assembling tracks..." << endl;
    map<int,vector<int> > tracks;
    
    int is = VERBOSE ? 0 : 1; // track 0 holds the unassigned points
    for(int track=is; track<=nt; track++) {
        vector<int> t;
        for (int j=0;j<nhits;j++) {
            if (track != label[j]) continue;
            t.push_back(j); // Save the results
        }
        tracks[track] = t;
    }
    
    cout << endl << "Number of tracks: " << nt << endl;
    int i = 0;
    for (auto it : tracks) {
        auto track = it.second;
        if (track.size() == 0) continue;
        if (i<MAXTRACK || i>nt-MAXTRACK) {
            cout << "Track " << it.first << ": ";
            for (auto it : track) cout << it << " ";
            cout << endl;
        }
        if (i++ == MAXTRACK) cout << endl << "..." << endl;
    }
    
    // Check the results
    checkTracks(tracks);
    
    // Generate a training sample for hit pairs and triples
    if (TRAINFILE) trainNetworks(base_path,filenum);
    
    // Show the results in a canvas
    if (DRAW) draw(nhits,x,y,z,tracks);
}

// Check the track hits (evaluation has ascending order)
long checkTracks(map<int,vector<int> >  &tracks) {
    cout << endl << "Checking tracks: " << endl;
    long error = 0;
    int n = 0;
    for (auto it : tracks) {
        if (it.first==0) continue; // track 0 holds the unassigned hits
        if (!VERBOSE && n++>MAXTRACK) break;
        if (it.second.size()==0) continue;
        auto t = it.second;
        long id = t[0];
        long errorid = 0;
        for (auto index : t) {
            if (index != id++) {
                id = index;
                if (errorid == 0) errorid = index;
                error++;
            }
        }
        
        if (errorid != 0) {
            cout << "Track " << it.first << ": ";
            for (auto index : t) {
                if (index == errorid)
                    cout << ">" << index << "< ";
                else
                    cout << index << " ";
            }
            cout << endl;
        }
    }
    return error;
}


void draw(long nhits,float *x,float *y,float *z,map<int,vector<int> > tracks)
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
        view->SetRange(-2000,-2000,-2000,2000,2000,2000); // draw in a 2 meter cube
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
        for (int i=1;i<nt;i++) {
            //cout << endl << "Drawing track " << i+1 << ": ";
            vector<int> track = tracks[i];
            int n = 0;
            TPolyLine3D *connector = new TPolyLine3D((int)track.size());
            for (auto &it : track)    {
                Point hit = Tracker::points[it];
                connector->SetPoint(n++, hit.x(), hit.y(), hit.z());
            }
            connector->SetLineWidth(1);
            connector->SetLineColor(kRed);
            connector->Draw();
        }
        
        cout <<  "Saving canvas file " << cname << endl;
        c1->Write();
        output.Close();
    }
    
}


