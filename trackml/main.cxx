// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include "Graph.h"
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

#define MAXPARTICLES 10
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

void makeTrain2();
void makeTrain2pairs();
void makeTrain2PhiTheta();
void makeTrain3();
void makeTrain3triples();
void makeTrain3seed();
void makeTrain3Continuous();
long checkTracks(std::map<int,std::vector<int> >  &tracks);
void draw(long nhits,float *x,float *y,float *z,map<int,vector<int> > tracks);
void transform(Particle &particle, std::vector<Point> &points);

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
    
    Tracker::verbose(VERBOSE);
    
    if (EVALUATION) {
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
    int label[nhits],truth[nhits],layer[nhits],module[nhits];
    
    nhits = 0;
    int n = 0;
    int start[nParticles+1],end[nParticles+1];
    start[0] = 0;
    end[0] = -1;
    for (auto &track : Tracker::truth_tracks) {
        vector<int> t = track.second;
        point geo = Tracker::meta[t[0]]; // Check the first layer of a hit
        int vol = geo.x;
        int lay = geo.y;
        int first = Tracker::getLayer(vol,lay);
        //if (first!=0 && first!=4 && first!=11) continue; // track does not start at first layers
        if (n++ >= MAXPARTICLES) break;
        start[n] = end[n-1]+1;
        end[n] = end[n-1] + (int)t.size();
        if (VERBOSE) cout << "Track  " << n << " {";
        int oldl = -1;
        int oldindex = -1;
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
            int mod = geo.z;
            int l = Tracker::getLayer(vol,lay);
            layer[nhits] = l;
            module[nhits] = mod;
            int index = MODULES*l + mod;
            // Add the hit pair to the paths graph
            if (oldindex>-1 && oldindex!=index) {
                Point p1(x[nhits-1],y[nhits-1],z[nhits-1]);
                Point p2(x[nhits],y[nhits],z[nhits]);
                double recall = Tracker::recall2(p1, p2)[0];
                //if (recall > THRESHOLD2)
                Tracker::paths.add(oldindex,index,1000*recall);
            }
            if (VERBOSE) cout << "{" << index << ","  << l << "," << mod << "},";
            oldl = l;
            oldindex = index;
            nhits++;
        }
        Tracker::paths.add(oldindex,-1);
        if (VERBOSE) {
            cout << "-1}" << endl;
            if (n<100) cout << "Track " << n << ": " << start[n] << "-" << end[n] << endl;
        }
    }
    if (VERBOSE) {
        cout << Tracker::paths << endl;
        auto modpath = serialize(Tracker::paths);
        cout << "modpath:" << endl;
        for (auto &it : modpath) Tracker::print(it.second);
    }
    
    if (nhits > MAXHITS) nhits = MAXHITS;
    cout << "Hits: " << nhits << endl;
    
    cout << endl << "Running Tracker:" << endl;
    long nt = Tracker::findTracks((int)nhits,x,y,z,layer,module,label,truth);
    
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
    
    // Generate a training sample for hit pairs
    if (TRAINFILE) {
        TString filePrefix;
        filePrefix.Form("%sevent%09d",base_path.c_str(),filenum);
        TString fname = filePrefix+".root";
        auto f = TFile::Open(fname,"RECREATE");
        cout << endl << "Generating training data file " << fname << endl;
        ntuple2 = new TNtuple("tracks","training data","rz1:phi1:z1:rz2:phi2:z2:l1:l2:truth");
        makeTrain2pairs();
        ntuple2->Write();
        ntuple3 = new TNtuple("tracks3","training data","rz1:phi1:z1:rz2:phi2:z2:rz3:phi3:z3:l1:l2:l3:truth");
        makeTrain3triples();
        ntuple3->Write();
        delete ntuple2;
        delete ntuple3;
        f->Close();
    }
    
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

void transform(Particle &particle, std::vector<treePoint> &points) {
    
    vector<treePoint> tmpvec;
    long nhits = (long)particle.hit.size();
    static int trackid = 0;
    
    trackid++;
    
    for (int i=0;i<nhits;i++) {
        vector<int> &h = particle.hit;
        int id = h[i];
        point h1 = Tracker::hits[id]*0.001; // in m
        treePoint p(h1.x,h1.y,h1.z,id,trackid,i);
        tmpvec.push_back(p);
    }
    
    sort(tmpvec.begin(),tmpvec.end(),Point::sortRz);
    points.insert(points.end(),tmpvec.begin(),tmpvec.end());
}

// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain2pairs()
{
    long wright=1,wrong=1;
    for (auto pair : Tracker::pairs) {
        treePoint &p1 = Tracker::points[pair.first];
        treePoint &p2 = Tracker::points[pair.second];
        int l1 = p1.layer();
        int l2 = p2.layer();
        if (p1.truth() == p2.truth()) {
            ntuple2->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),l1,l2,1.0); //wright combination
            wright++;
        }
        else {
            if (r.Rndm()<wright/wrong) {
                ntuple2->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),l1,l2,0.0); //wrong combination
                wrong++;
            }
        }
    }
    cout << "makeTrain2pairs Wright: " << wright << " Wrong: " << wrong << endl;
}

// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain2()
{
    
    const int n=50; // hit pair layer combinations
    pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
    
    long wright=1,wrong=1;
    for (int i = 0; i < n; i++) {
        int tube1 = start_list[i].first;
        for (auto start : Tracker::modules[tube1]) { // all modules in first layer
            const auto edgelist = Tracker::paths.edges(start);
            if (edgelist.size() == 0) continue;
            int l1 = start/MODULES;
            //int m1 = start%MODULES;
            for (auto edge : edgelist) {
                int nextindex = edge.first;
                if (nextindex<0) break;
                int l2 = nextindex/MODULES;
                //int m2 = nextindex%MODULES;
                for (auto a : Tracker::module[start]) { // all hits in module
                    treePoint &p1 = Tracker::points[a];
                    for (auto &b : Tracker::module[nextindex]) {
                        treePoint &p2 = Tracker::points[b];
                        if (p1.truth() == p2.truth()) {
                            ntuple2->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),l1,l2,1.0); //wright combination
                            wright++;
                        }
                        else {
                            if (r.Rndm()<wright/wrong) {
                                ntuple2->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),l1,l2,0.0); //wrong combination
                                wrong++;
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "makeTrain2 Wright: " << wright << " Wrong: " << wrong << endl;
}
// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain2PhiTheta()
{
    
    const int n=50; // hit pair layer combinations
    pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
    
    long wright=1,wrong=1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j< PHIDIM; j++) {
            for (int k = 0; k<THEDIM; k++) {
                int tube1 = start_list[i].first;
                for (auto &a : Tracker::tube[tube1][j][k]) {
                    treePoint &p1 = Tracker::points[a];
                    int tube2 = start_list[i].second;
                    for (auto &b : Tracker::tube[tube2][j][k]) {
                        treePoint &p2 = Tracker::points[b];
                        if (p1.truth() == p2.truth()) {
                            ntuple2->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),tube1,tube2,1.0); //wright combination
                            wright++;
                        }
                        else {
                            if (r.Rndm()<wright/wrong) {
                                ntuple2->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),tube1,tube2,0.0); //wrong combination
                                wrong++;
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "makeTrain2PhiTheta Wright: " << wright << " Wrong: " << wrong << endl;
}

// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain3triples()
{
    long wright=1,wrong=1;
    for (auto triple : Tracker::triples) {
        treePoint &p1 = Tracker::points[triple.x];
        treePoint &p2 = Tracker::points[triple.y];
        treePoint &p3 = Tracker::points[triple.z];
        int l1 = p1.layer();
        int l2 = p2.layer();
        int l3 = p3.layer();
        if (p1.truth() == p2.truth() && p1.truth() == p3.truth()) {
            ntuple3->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),p3.rz(),p3.phi(),p3.z(),l1,l2,l3,1.0); //wright combination
            wright++;
        }
        else {
            if (r.Rndm()<wright/wrong) {
                ntuple3->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),p3.rz(),p3.phi(),p3.z(),l1,l2,l3,0.0); //wrong combination
                wrong++;
            }
        }
    }
    cout << "makeTrain3triples Wright: " << wright << " Wrong: " << wrong << endl;
}


// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain3()
{
    long wright=1,wrong=1;
    // Combine 3 hits
    for (auto p : Tracker::particles) {
        long nhits = p.hits;
        if (nhits < 3) continue;
        vector<int> &h = p.hit;
        for (int i=0; i<nhits-3; i++)    {
            int id1 = h[i];
            int id2 = h[i+1];
            point x1 = Tracker::hits[id1]*0.001; // in m
            point x2 = Tracker::hits[id2]*0.001; // in m
            Point p1(x1.x,x1.y,x1.z);
            Point p2(x2.x,x2.y,x2.z);
            point geo = Tracker::meta[id1];
            int vol = geo.x;
            int lay = geo.y;
            int mod = geo.z;
            int l1 = Tracker::getLayer(vol,lay);
            geo = Tracker::meta[id2];
            vol = geo.x;
            lay = geo.y;
            mod = geo.z;
            int l2 = Tracker::getLayer(vol,lay);
            
            // Select 3 continuous points
            int id3 = h[i+2];
            point x3 = Tracker::hits[id3]*0.001; // in m
            Point p3(x3.x,x3.y,x3.z);
            geo = Tracker::meta[id3];
            vol = geo.x;
            lay = geo.y;
            mod = geo.z;
            int l3 = Tracker::getLayer(vol,lay);
            ntuple3->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),p3.rz(),p3.phi(),p3.z(),l1,l2,l3,1.0); //wright combination
            wright++;
            // Select last point randomly in the same module
            int index = MODULES*lay + mod;
            auto hits = Tracker::module[index];
            for (auto idr : hits) {
                if (idr==id3) continue; // Same point
                point x4 = Tracker::hits[idr]*0.001; // in m
                Point p4(x4.x,x4.y,x4.z);
                if (r.Rndm()<wright/wrong) {
                    ntuple3->Fill(p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),p4.rz(),p4.phi(),p4.z(),l1,l2,l3,0.0); //wrong combination
                    wrong++;
                }
            }
        }
    }
    cout << "makeTrain3: " <<  Tracker::particles.size() << " particles. " << "Wright: " << wright << " Wrong: " << wrong << endl;
}

// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain3seed()
{
    long wright=1,wrong=1;
    // Combine 3 hits
    for (auto p : Tracker::particles) {
        vector<treePoint> hits;
        transform(p,hits);
        
        // Sort the hits according to distance from origin
        sort(hits.begin(),hits.end(),Point::sortRad);
        
        // Combine 3 hits
        int nhits = (int)hits.size();
        if (nhits < 3) return;
        
        treePoint &hit1 = hits[0];
        treePoint &hit2 = hits[1];
        
        int istart = 2;
        float d = hit1.distance(hit2); // CHeck for double hits
        if (d<TWINDIST) {
            hit2 = hits[2];
            istart = 3;
        }
        
        double l1(0),l2(0),l3(0);
        point geo = Tracker::meta[hit1.id()];
        int vol = geo.x;
        int lay = geo.y;
        l1 = Tracker::getLayer(vol,lay);
        geo = Tracker::meta[hit2.id()];
        vol = geo.x;
        lay = geo.y;
        l2 = Tracker::getLayer(vol,lay);
        
        for (int i=istart; i<nhits; i++)    {
            treePoint &hit3 = hits[i];
            geo = Tracker::meta[hit3.id()];
            vol = geo.x;
            lay = geo.y;
            l3 = Tracker::getLayer(vol,lay);
            
            ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),hit3.phi(),hit3.z(),l1,l2,l3,hit1.truth()+1); //true combination
            wright++;
            double phi3 = 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r3
            double theta3 = r.Rndm()*M_PI;
            double z3 = hit3.rz() * cos(theta3);
            ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),phi3,z3,l1,l2,l3,0.0); // wrong combination
            wrong++;
        }
        
    }
    cout << "makeTrain3seed: " <<  Tracker::particles.size() << " particles. " << "Wright: " << wright << " Wrong: " << wrong << endl;
}

// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain3Continuous()
{
    long wright=1,wrong=1;
    // Combine 3 hits
    for (auto p : Tracker::particles) {
        vector<treePoint> hits;
        transform(p,hits);
        
        // Sort the hits according to distance from origin
        sort(hits.begin(),hits.end(),Point::sortRad);
        
        // Combine 3 hits
        int nhits = (int)hits.size();
        if (nhits < 3) return;
        for (int i=0; i<nhits-2; i++)    {
            treePoint &hit1 = hits[i];
            treePoint &hit2 = hits[i+1];
            treePoint &hit3 = hits[i+2];
            
            double l1(0),l2(0),l3(0);
            point geo = Tracker::meta[hit1.id()];
            int vol = geo.x;
            int lay = geo.y;
            l1 = Tracker::getLayer(vol,lay);
            geo = Tracker::meta[hit2.id()];
            vol = geo.x;
            lay = geo.y;
            l2 = Tracker::getLayer(vol,lay);
            geo = Tracker::meta[hit3.id()];
            vol = geo.x;
            lay = geo.y;
            l3 = Tracker::getLayer(vol,lay);
            
            ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),hit3.phi(),hit3.z(),l1,l2,l3,hit1.truth()+1); //true combination
            wright++;
            double phi3 = 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r3
            double theta3 = r.Rndm()*M_PI;
            double z3 = hit3.rz() * cos(theta3);
            ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),phi3,z3,hit1.rz(),hit2.rz(),hit3.rz(),0.0); // wrong combination
            wrong++;
        }
    }
    cout << "makeTrain3Continuous: " <<  Tracker::particles.size() << " particles. " << "Wright: " << wright << " Wrong: " << wrong << endl;
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
