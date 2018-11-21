// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#include "TFile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "Tracker.h"

using namespace std;

TRandom r;
TNtuple *ntuple1,*ntuple2,*ntuple3,*ntuple4;

void transform(Particle &particle, std::vector<treePoint> &points) {
    
    vector<treePoint> tmpvec;
    long nhits = (long)particle.hit.size();
    static int trackid = 0;
    
    trackid++;
    
    for (int i=1;i<=nhits;i++) {
        vector<int> &h = particle.hit;
        int id = h[i];
        point h1 = Tracker::hits[id]; // in mm
        treePoint p(h1.x,h1.y,h1.z,id,trackid,i);
        tmpvec.push_back(p);
    }
    
    sort(tmpvec.begin(),tmpvec.end(),Point::sortByRz);
    points.insert(points.end(),tmpvec.begin(),tmpvec.end());
}

// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrainPairs()
{
    long wright=1,wrong=1;
    
    const int n = 100;
    pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
    
    const int features = 6;
    float feature[features];
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j< PHIDIM; j++) {
            for (int k = 0; k<THEDIM; k++) {
                
                for (auto a : Tracker::tube[start_list[i].first][j][k]) {
                    if (Tracker::tube[start_list[i].first][j][k].size()==0) continue;
                    for (auto b : Tracker::tube[start_list[i].second][j][k]) {
                        if (Tracker::tube[start_list[i].second][j][k].size()==0) continue;
                        treePoint &pa = Tracker::points[a];
                        treePoint &pb = Tracker::points[b];
                        double dot = pa.x()*pb.x()+pa.y()*pb.y();
                        double alen = dist2(pa.x(), pa.y());
                        double blen = dist2(pb.x(), pb.y());
                        if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;
                        
                        dot += pa.z()*pb.z();
                        alen += pa.z()*pa.z();
                        blen += pb.z()*pb.z();
                        if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;
                        
                        int g = Tracker::good_pair(pa, pb);
                        
                        if (g==0 && r.Rndm()>0.25*wright/wrong) continue;
                        wright += g!=0;
                        wrong  += g==0;
                        float l1 = pa.layer();
                        float l2 = pb.layer();
                        point v = Tracker::truth_pos[pa.id()];
                        if (Tracker::getFeatures3(pa, pb, feature)) {
                            float x[19]={pa.rz(),pa.phi(),pa.z(),pb.rz(),pb.phi(),pb.z(),feature[0],feature[1],feature[2],feature[3],feature[4],feature[5],l1,l2,(float)v.x,(float)v.y,(float)v.z,(float)g};
                            ntuple1->Fill(x);
                        }
                    }
                }
            }
        }
    }
    cout << "makeTrainPairs Wright: " << wright << " Wrong: " << wrong << endl;
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
                        if (p1.trackid() == p2.trackid()) {
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
                        if (p1.trackid() == p2.trackid()) {
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
void makeTrain3()
{
    long wright=1,wrong=1;
    for (auto triple : Tracker::triples) {
        treePoint &p1 = Tracker::points[triple.x];
        treePoint &p2 = Tracker::points[triple.y];
        treePoint &p3 = Tracker::points[triple.z];
        int l1 = p1.layer();
        int l2 = p2.layer();
        int l3 = p3.layer();
        if (p1.trackid() == p2.trackid() && p1.trackid() == p3.trackid()) {
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
    cout << "makeTrain3 Wright: " << wright << " Wrong: " << wrong << endl;
}

void initDensity3(void);
double scoreTriple(int ai, int bi, int ci);
double scoreTripleLogRadius_and_HitDir(int ai,int bi,int ci,double* L);
double scoreTripleDensity(int a,int b,int c);
int samepart(int a, int b);


// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrainTriples()
{
    long wright=1,wrong=1;

    initDensity3();
    
    for (auto &track : Tracker::particles) {
        vector<int> t = track.hit;
        if (t.size()<=2) continue;

        for (int i=0;i<track.hits-2;i++) {

            int id1 = t[i];
            int id2 = t[i+1];
            int id3 = t[i+2];
            int good = samepart(id1,id2) && samepart(id2,id3);
            if (good==0) continue;
        
            float f[7];
            double L[4];
            L[3] = exp(scoreTripleLogRadius_and_HitDir(id1,id2,id3,L));
            f[0] = L[0];
            f[1] = L[1];
            f[2] = L[2];
            f[3] = L[3];
            f[4] = scoreTriple(id1,id2,id3);
            f[5] = scoreTripleDensity(id1,id2,id3);
            f[6] = scoreTripleDensity(id3,id2,id1);
            //cout << good << ' ' << A << ' ' << B << ' ' << C << ' ' << D << ' ' << f[0] << ' ' << f[1] << ' ' << f[2] << endl;
            for (int j=0;j<7;j++) {
                if (isnan(f[j])) {
                    cout << "makeTrainTriples: NAN " << j << endl;
                }
            }

            point v = Tracker::truth_pos[id1];
            treePoint &p1 = Tracker::points[id1];
            treePoint &p2 = Tracker::points[id2];
            treePoint &p3 = Tracker::points[id3];
            float l1 = p1.layer();
            float l2 = p2.layer();
            float l3 = p3.layer();

            float x[23]={p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),p3.rz(),p3.phi(),p3.z(),f[0],f[1],f[2],f[3],f[4],f[5],f[6],l1,l2,l3,(float)v.x,(float)v.y,(float)v.z,1.0};
            ntuple4->Fill(x); //wright combination
            wright++;

            int index = MODULES*p3.layer() + p3.module();
            for (auto i : Tracker::module[index]) {
                if (r.Rndm()>1.5*wright/wrong) continue;
                treePoint &p3 = Tracker::points[i];
                int idr = p3.id();
                if (idr==id1 || idr==id2 || idr==id3) continue; // Do not take the same hit
                L[3] = exp(scoreTripleLogRadius_and_HitDir(id1,id2,id3,L));
                f[0] = L[0];
                f[1] = L[1];
                f[2] = L[2];
                f[3] = L[3];
                f[4] = scoreTriple(id1,id2,idr);
                f[5] = scoreTripleDensity(id1,id2,idr);
                f[6] = scoreTripleDensity(id1,id2,idr);
                float x[23]={p1.rz(),p1.phi(),p1.z(),p2.rz(),p2.phi(),p2.z(),p3.rz(),p3.phi(),p3.z(),f[0],f[1],f[2],f[3],f[4],f[5],f[6],l1,l2,l3,(float)v.x,(float)v.y,(float)v.z,0.0};
                ntuple4->Fill(x); //wrong combination
                wrong++;
            }
        }
    }
    cout << "makeTrainTriples: " <<  Tracker::particles.size() << " particles. " << "Wright: " << wright << " Wrong: " << wrong << endl;
}


// Look for seeding points by hit pair combinations in the innnermost layers
void makeTrain3PhiTheta()
{
    long wright=1,wrong=1;
    
    // Combine 3 hits
    for (auto p : Tracker::particles) {
        vector<treePoint> hits;
        transform(p,hits);
        
        // Sort the hits according to distance from origin
        sort(hits.begin(),hits.end(),Point::sortByRadius);
        
        // Combine 3 hits
        int nhits = (int)hits.size();
        if (nhits < 3) return;
        for (int i=1; i<=nhits-2; i++)    {
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
            
            ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),hit3.phi(),hit3.z(),l1,l2,l3,(float)hit1.trackid()+1); //true combination
            wright++;
            double phi3 = 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r3
            double theta3 = r.Rndm()*M_PI;
            double z3 = hit3.rz() * cos(theta3);
            ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),phi3,z3,hit1.rz(),hit2.rz(),hit3.rz(),0.0); // wrong combination
            wrong++;
        }
    }
    cout << "makeTrain3PhiTheta: " <<  Tracker::particles.size() << " particles. " << "Wright: " << wright << " Wrong: " << wrong << endl;
}

// Generate a training sample for hit pairs and triples
void trainNetworks(string base_path,int filenum) {
    TString filePrefix;
    filePrefix.Form("%sevent%09d",base_path.c_str(),filenum);
    TString fname = filePrefix+".root";
    auto f = TFile::Open(fname,"RECREATE");
    cout << endl << "Generating training data file " << fname << endl;
    ntuple1 = new TNtuple("pairs","training data","rz1:phi1:z1:rz2:phi2:z2:f0:f1:f2:f3:f4:f5:l1:l2:vx:vy:vz:truth");
    ntuple2 = new TNtuple("tracks","training data","rz1:phi1:z1:rz2:phi2:z2:l1:l2:truth");
    makeTrainPairs();
    ntuple1->Write();
    ntuple2->Write();
    ntuple3 = new TNtuple("tracks3","training data","rz1:phi1:z1:rz2:phi2:z2:rz3:phi3:z3:l1:l2:l3:truth");
    ntuple4 = new TNtuple("triples","training data","rz1:phi1:z1:rz2:phi2:z2:rz3:phi3:z3:f0:f1:f2:f3:f4:f5:f6:l1:l2:l3:vx:vy:vz:truth");
    makeTrainTriples();
    ntuple3->Write();
    ntuple4->Write();
    delete ntuple1;
    delete ntuple2;
    delete ntuple3;
    delete ntuple4;
    f->Close();
}


