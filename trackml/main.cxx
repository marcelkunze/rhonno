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

#define VERBOSE false
#define MAXPAIRS 10
#define MAXPARTICLES 10000
#define MAXHITS 150000
#define TRAINFILE true
#define DRAW true
#define EVALUATION true
#define MAXTRACK 10
#define MAXLABEL 100

const std::string base_path = "/Users/marcel/workspace/train_sample/";

//Which event to run, this may be overwritten by main()'s arguments
int filenum = 21100;
int maxpairs = MAXPAIRS, maxparticles = MAXPARTICLES;
bool verbose = VERBOSE;

using namespace std;

long checkTracks(map<int,vector<int> >  &tracks);
void trainNetworks(string base_path,int filenum);
void draw(long nhits,float *x,float *y,float *z,map<int,vector<int> > tracks);
bool samepart(int,int);

int main(int argc, char**argv) {
    //Read which event to run from arguments, default is event # 1000
    //Supports events 0-124 for test, or 1000-1099 for validation (small dataset)
    
    if (argc >= 2) {
        filenum = atoi(argv[1]);
        cout << "Running on event #" << filenum << endl;
    }
    if (argc >= 3) {
        maxparticles = atoi(argv[2]);
        cout << "Running on event #" << filenum << endl;
    }
    if (argc >= 4) {
        verbose = atoi(argv[3])!=0;
        cout << "Running on event #" << filenum << endl;
    }
    ios::sync_with_stdio(false);
    cout << fixed;
    
    Tracker::verbose(verbose);
    
    if (EVALUATION) {
        //Tracker::readGraph("paths.csv",Tracker::paths);
        Tracker::readParticles(base_path,filenum);
        Tracker::readTruth(base_path,filenum);
        //Tracker::sortTracks();
    }
    Tracker::readHits(base_path,filenum);
    Tracker::readCells(base_path,filenum);
    
    long nParticles = Tracker::truth_tracks.size();
    if (nParticles>maxparticles) nParticles = maxparticles;
    cout << "Particles: " << nParticles << endl;
    
    // Use 1 indexing to be compatible to trackml hit index
    long nhits = Tracker::hits.size();
    float x[nhits+1],y[nhits+1],z[nhits+1],cx[nhits+1],cy[nhits+1],cz[nhits+1];
    int label[nhits+1],volume[nhits+1],layer[nhits+1],module[nhits+1],hitid[nhits+1];
    long long trackid[nhits+1];
    
    for (int i=1;i<=nhits;i++) {
        int hit_id = i+1;
        point &hit = Tracker::hits[hit_id];
        x[i] = hit.x; // in mm
        y[i] = hit.y; // in mm
        z[i] = hit.z; // in mm
        label[i] = 0;
        trackid[i] = Tracker::truth_part[hit_id]; // true track assignment
        hitid[i] = hit_id;
        point geo = Tracker::meta[hit_id];
        int vol = geo.x;
        int lay = geo.y;
        int mod = geo.z;
        int l = Tracker::getLayer(vol,lay);
        volume[i] = vol;
        layer[i] = l;
        module[i] = mod;
        cx[i] = Tracker::hit_dir[hit_id][0].x;
        cy[i] = Tracker::hit_dir[hit_id][0].y;
        cz[i] = Tracker::hit_dir[hit_id][0].z;
    }
    
    // Prepare the trackml data to run the track finder
    // Geberate a graph to represent the track hits in the modules
    
    int tracknumber = 1;
    for (auto &track : Tracker::particles) {
        vector<int> t = track.hit;
        if (t.size()==0) continue;
        Tracker::truepairs.push_back(make_pair(t[0],t[1]));
        point geo = Tracker::meta[t[0]]; // Check the first layer of a hit
        int vol = geo.x;
        int lay = geo.y;
        int first = Tracker::getLayer(vol,lay);
        //if (first!=0 && first!=4 && first!=11) continue;
        
        int oldindex = -1;
        for (auto &hit_id : t) {
            Tracker::truth_assignment[hit_id] = tracknumber;
            point geo = Tracker::meta[hit_id];
            int vol = geo.x;
            int lay = geo.y;
            int mod = geo.z;
            int l = Tracker::getLayer(vol,lay);
            int index = MODULES*l + mod;
            // Add the hit pair to the paths graph
            if (oldindex>-1 && oldindex!=index) {
                Tracker::paths.add(oldindex,index,1.0,true); // incremental mode
            }
            oldindex = index;
        }
        tracknumber++;
        Tracker::paths.add(oldindex,-1);
        if (verbose) {
            //cout << "-1}" << endl;
        }
    }

    // Write path data to file
    //if (verbose) Tracker::paths.print();
    Tracker::writeGraph("paths.csv",Tracker::paths);
    
    if (nhits > MAXHITS) nhits = MAXHITS;
    cout << "Hits: " << nhits << endl;
    cout << "maxparticles: " << maxparticles << " maxpairs: " << maxpairs << endl;
    
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
    
    int is = verbose ? 0 : 1; // track 0 holds the unassigned points
    for(int track=is; track<=nt; track++) {
        vector<int> t;
        for (int j=0;j<nhits;j++) {
            if (track != label[j]) continue;
            t.push_back(j); // Save the results
        }
        tracks[track] = t;
    }
    
    cout << endl << "Number of tracks: " << nt << endl;
    for (auto it : tracks) {
        if (it.first ==0) continue; // Holds unassigned points
        auto track = it.second;
        if (track.size() == 0) continue;
        if (it.first<MAXTRACK || it.first>nt-MAXTRACK) {
            cout << "Track " << it.first << ": ";
            for (auto it : track) cout << it << "(" << Tracker::truth_assignment[it] << ") ";
            cout << endl;
        }
        if (it.first == MAXTRACK) cout << endl << "..." << endl;
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
        if (!verbose && n++>MAXPAIRS) break;
        if (it.second.size()==0) continue;
        auto t = it.second;
        long errorid = 0;
        for (auto index : t) {
            if (!samepart(t[0],index)) {
                if (errorid == 0) errorid = index;
                error++;
            }
        }
        
        if (errorid != 0) {
            cout << "Track " << it.first << ": ";
            for (auto index : t) {
                int tracknumber = Tracker::truth_assignment[index];
                if (index == errorid)
                    cout << ">" << index << "(" << tracknumber << ")< ";
                else
                    cout << index << "(" << tracknumber << ") ";
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


