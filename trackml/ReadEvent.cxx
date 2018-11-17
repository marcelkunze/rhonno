// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#define DRAW true

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TView.h"
#include "TPolyMarker3D.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"

#include <iostream>
#include <fstream>
#include <random>

#include "DataStructures.h"
#include "Geo.h"
#include "XMLP.h"
#include "Tracker.h"
#include "Point.h"

using namespace std;

std::vector<Hit> mHits;
std::vector<HitMC> mHitsMC;
std::vector<xParticle> mParticles;

int layerNHits[Geo::NLayers];

#define NEVENTS 10
#define MAXHITS 150000
#define MAXPARTICLES 10000
#define MCHITS false
//#define FINDTRACKS
#define NHITS 10
#define SIGMA 0.001
#define VERBOSE false
#define TRACKML
#define DELTAPHI 0.1
#define DELTATHE 0.1

TRandom r;

void transform(xParticle &particle, std::vector<treePoint> &points, bool mc);
void combine(xParticle &p1, xParticle &p2);
void combine2(xParticle &p1, xParticle &p2);
void combine3(xParticle &p1, xParticle &p2);
void combine2(xParticle &p1);
void combine3(xParticle &p1);
void combine3serial(xParticle &p1);
void readEvent( const char *directory, int event, bool loadMC );
void GenerateTrack(std::vector<treePoint> &points, int np, double delta, double radius, double phi, double gamma, double error);

map<unsigned long,int> particlemap;
map<unsigned long,int> recomap;
TNtuple *ntuple2,*ntuple3;

int main()
{
#ifdef TRACKML
    bool analyseTruth = true;
    
    const int nEvents = NEVENTS;
    //const int firstEvent=100021100;
    const int firstEvent=21100;
    //const int firstEvent=1000;
    TString dir = "/Users/marcel/workspace/train_sample/";
    
    ofstream out("mysubmission.csv");
    if( !out.is_open() ){
        cout<<"Can not open output file"<<endl;
        exit(0);
    }
    
    ofstream outtracks("tracks.csv");
    if( !outtracks.is_open() ){
        cout<<"Can not open tracks file"<<endl;
        exit(0);
    }
    
    out<<"event_id,hit_id,track_id"<<endl;
    outtracks<<"event_id,track_id: hits"<<endl;
    
    for (int event = firstEvent; event<firstEvent+nEvents; event++) {
        cout<<"read event "<<event<<endl;
        readEvent( dir.Data(),  event, analyseTruth );
        
        TString filePrefix;
        filePrefix.Form("%sevent%09d",dir.Data(),event);
        TString fname = filePrefix+".root";
        auto f = TFile::Open(fname,"RECREATE");
        cout << "Writing training data to " << fname << endl;
        ntuple2 = new TNtuple("tracks","training data","rz1:phi1:z1:rz2:phi2:z2:l1:l2:truth");
        ntuple3 = new TNtuple("tracks3","training data","rz1:phi1:z1:rz2:phi2:z2:rz3:phi3:z3:l1:l2:l3:truth");
        
        TH1F h0("h0","Track length",100,0,100);
        TH1F h1("h1","x distribution",100,-3,3);
        TH1F h2("h2","y distribution",100,-3,3);
        TH1F h3("h3","z distribution",100,-3,3);
        TH1F h4("h4","rz distribution",100,0,3);
        TH1F h5("h5","r distribution",100,0,3);
        TH1F h6("h6","phi distribution",100,-4,4);
        TH1F h7("h7","theta distribution",100,-4,4);
        
        for (int ih=0; ih<mHitsMC.size(); ih++ ) {
            HitMC &h = mHitsMC[ih];
            if (ih<10) h.Print();
        }
        
        for (int ih=0; ih<mHits.size(); ih++ ) {
            Hit &h = mHits[ih];
            if (ih<10) h.Print();
        }
        
        long start[MAXPARTICLES+1];
        long end[MAXPARTICLES+1];
        std::vector<treePoint> hits;
        hits.reserve(MAXHITS);
        
        size_t nParticles = mParticles.size();
        if (nParticles > MAXPARTICLES) nParticles = MAXPARTICLES;
        cout << "Particles:" << nParticles << endl;
        
        start[0] = 0;
        end[0] = -1;
        map<int,int> particles;
        for (int ip=0; ip<nParticles; ip++ ) {
            xParticle &p1 = mParticles[ip];
            start[ip+1] = end[ip]+1;
            end[ip+1] = end[ip] + p1.hits.size();
        }
        
        unsigned long nh = 0;
        for (int ip=0; ip<nParticles; ip++ ) {
 
            xParticle &p1 = mParticles[ip];
            h0.Fill(p1.hits.size());

            if (VERBOSE && ip<1000) cout << endl << "Track " << ip+1 << "(" << start[ip+1] << "-" << end[ip+1] << "): ";
            
            // Generate a particle hit map to check the reconstruction results
            for (int j=0; j<p1.hits.size(); j++) {
                Hit &h = mHits[p1.hits[j]];
                h.trackID = ip;
                particlemap[nh++] = ip;
                treePoint p(h.x,h.y,h.z);
                h1.Fill(p.x());
                h2.Fill(p.y());
                h3.Fill(p.z());
                h4.Fill(p.rz());
                h5.Fill(p.r());
                h6.Fill(p.phi());
                h7.Fill(p.theta());
                if (VERBOSE) cout << h.hitID << "(" << nh << ") ";
            }
        
            // Push the particle hits into Point vector
            transform(p1,hits,MCHITS);

            // Generate training data
            //int index = ip;
            //while (ip == index) index = rand()%nParticles;
            //Particle &p2 = mParticles[index];
            combine2(p1); // Produce training data (2 hits)
            combine3serial(p1); // Produce training data (3 hits)
            
        }
        
        f->Write();
        delete ntuple2;
        delete ntuple3;
        
#else
        {
            int event = 0;
            std::vector<treePoint> hits;
            hits.reserve(MAXHITS);
            // std::vector<treePoint>, int np, float delta tau, float radius, float phi, float gamma
            GenerateTrack(hits,NHITS,0.025, 1.0,M_PI/2.0, 1.0,SIGMA); // 00
            GenerateTrack(hits,NHITS,0.025,-1.0,M_PI/2.0, 1.0,SIGMA); // 10
            GenerateTrack(hits,NHITS,0.025, 1.5,M_PI/1.0,-1.0,SIGMA); // 20
            GenerateTrack(hits,NHITS,0.025,-2.0,M_PI/3.0,-1.0,SIGMA); // 30
#endif
            size_t nhits = hits.size();
            if (nhits > MAXHITS) nhits = MAXHITS;
            cout << "Hits: " << nhits << endl;
            
            float x[nhits],y[nhits],z[nhits],cx[nhits],cy[nhits],cz[nhits];
            int label[nhits],volume[nhits],layer[nhits],module[nhits],hitid[nhits];
            long long truth[nhits];
            int nt(0);
            
            double minr = FLT_MAX;
            double maxr = 0.0;
            for (int i=0; i<nhits; i++) {
                treePoint &h = hits[i];
                Hit &hit = mHits[h.id()];
                x[i] = h.x();
                y[i] = h.y();
                z[i] = h.z();
                cx[i] = 0.0;
                cy[i] = 0.0;
                cz[i] = 0.0;
                hitid[i] = h.id();
                label[i] = 0;
                volume[i] = hit.volume;
                layer[i] = Tracker::getLayer(hit.volume,hit.layer);
                module[i] = hit.module;
                if (h.r() < minr) minr = h.r();
                if (h.r() > maxr) maxr = h.r();
            }
            cout << "Min. radius: " << minr << endl;
            cout << "Max. radius: " << maxr << endl;
 
#ifdef FINDTRACKS
            cout << "Find tracks..." << endl;
            nt = Tracker::findTracks((int)nhits,x,y,z,cx,cy,cz,volume,layer,module,hitid,truth,label);
#endif
            // Assemble the tracks from label information
            nh = 0;
            map<int,vector<treePoint> > tracks;
            for(int i=0; i<nt; i++) {
                int track = i+1;
                vector<treePoint> t;
                for (int j=0;j<nhits;j++) {
                    if (track != label[j]) continue;
                    hits[j].setlabel(label[j]);
                    t.push_back(hits[j]); // Save the results
                    recomap[nh++] = i;
                }
                tracks[i] = t;
            }
            
            // Compare the results with truth
            cout << "Check the results..." << endl;
            
            nh = 0;
            unsigned long np = 0;
            int delta = 0;
            for(int ip=0; ip<nt; ip++) {
                vector<treePoint> t = tracks[ip];
                for (int j=0;j<t.size();j++) {
                    auto itp = particlemap.find(np++);
                    int trackp = itp->second + delta;
                    auto ith = recomap.find(nh++);
                    int trackh = ith->second;
                    if (trackp == trackh) {
                        if (VERBOSE) cout << nh << ":" << trackh << "/" << trackp << endl;
                        particlemap.erase(itp++);
                    }
                    else {
                        if (VERBOSE) cout << nh << ":" << trackh << "/" << trackp << " NOK" << endl;
                        while (trackp<trackh) { trackp++; delta++; }
                        while (trackp>trackh) { trackp--; delta--; }
                    }
                }
            }
            
            if (particlemap.size()>0) cout << "Wrongly assigned:" << particlemap.size() << endl;
            
#define MAXTRACK 10
            cout << endl << "Number of tracks:" << nt << endl;
            for (int i=0; i<nt; i++) {
                vector<treePoint> t = tracks[i];
                if (i == MAXTRACK) cout << endl << "..." << endl;
                if (i<MAXTRACK || i>nt-MAXTRACK) {
                    cout << "Track " << i+1 << ": ";
                    for (vector<treePoint>::iterator it = t.begin(); it != t.end(); ++it) {
                        treePoint p = *it;
                        cout << p.id() << " ";
                    }
                    cout << endl;
                }
            }
            
#define MAXLABEL 100
            cout << "Labels: ";
            for (int i=0;i<nhits;i++) {
                if (i<MAXLABEL || i>nhits-MAXLABEL) cout << label[i] << " ";
                if (i == MAXLABEL) cout << endl << "..." << endl;
            }
            cout << endl;
            
#ifdef TRACKML
            cout << "Write tracks file..." << endl;
            for (int ip=0; ip<mParticles.size(); ip++ ) {
                xParticle &p1 = mParticles[ip];
                outtracks<<event<<","<<ip<<": ";
                for (int j=0;j<p1.hits.size();j++) outtracks<<p1.hits[j]<<" ";
                outtracks << endl;
            }
            outtracks.close();
            
            cout << "Write submission file..." << endl;
            for (int ih=0; ih<nhits; ih++ ){
                out<<event<<","<<ih+1<<","<<label[ih]<<endl;
            }
            out.close();
#endif
            // Initialize a 3D canvas and draw the hits and tracks
            if (DRAW) {
                TString dir = "/Users/marcel/workspace/train_sample/";
                TString filePrefix;
                filePrefix.Form("%sevent%09d",dir.Data(),event);
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
                for (auto &p : hits)    {
                    static int i = 0;
                    hitmarker->SetPoint(i++,p.x(),p.y(),p.z());
                }
                // set marker size, color & style
                hitmarker->SetMarkerSize(1.0);
                hitmarker->SetMarkerColor(kCyan);
                hitmarker->SetMarkerStyle(kStar);
                hitmarker->Draw();
                
                // Draw the tracks
                for (int i=0;i<nt;i++) {
                    //cout << endl << "Drawing track " << i+1 << ": ";
                    vector<treePoint> track = tracks[i];
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
        
        return 0;
}

void transform(xParticle &particle, std::vector<treePoint> &points, bool mc=false) {
    
    vector<treePoint> tmpvec;
    int nhits = (int)particle.hits.size();
    static int trackid = 0;
    
    trackid++;
    
    if (!mc) {
        for (int i=0;i<nhits;i++) {
            Hit &h1 = mHits[particle.hits[i]];
            treePoint p(h1.x,h1.y,h1.z,0.,0.,0.,h1.hitID,trackid,h1.id);
            tmpvec.push_back(p);
        }
    }
    else {
        for (int i=0;i<nhits;i++) {
            HitMC &h1 = mHitsMC[particle.hits[i]];
            //Hit &h2 = mHits[particle.hits[i]];
            treePoint p(h1.x,h1.y,h1.z,0.,0.,0.,h1.hitID+1,trackid,h1.id);
            tmpvec.push_back(p);
        }
    }
    
    sort(tmpvec.begin(),tmpvec.end(),Point::sortRz);
    points.insert(points.end(),tmpvec.begin(),tmpvec.end());
}

void combine(xParticle &p1, xParticle &p2)
{
    int nhits1 = (int)p1.hits.size();
    int nhits2 = (int)p2.hits.size();
    
    vector<treePoint> hits1, hits2;
    transform(p1,hits1,MCHITS);
    transform(p2,hits2,MCHITS);
    
    for (int i=0; i<nhits1; i++)    {
        Hit &h1 = mHits[p1.hits[i]];
        treePoint &hit1 = hits1[i];
        for (int j=0; j<nhits1; j++)    {
            Hit &h2 = mHits[p1.hits[j]];
            treePoint &hit2 = hits1[j];
            if (i!=j) ntuple2->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),h1.values,h2.values,hit1.truth()+1); //true combination
        }
        for (int j=0; j<nhits2; j++)    {
            Hit &h2 = mHits[p2.hits[j]];
            treePoint &hit2 = hits2[j];
            if (i!=j) ntuple2->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),h1.values,h2.values,0.0); // wrong combination
        }
    }
}

void combine2(xParticle &p)
{
    vector<treePoint> hits;
    transform(p,hits,MCHITS);
    
    // Sort the hits according to distance from origin
    sort(hits.begin(),hits.end(),Point::sortRad);

    int nhits = (int)hits.size();
    if (nhits < 2) return;
    for (int i=0; i<nhits-1; i++)    {
        treePoint &hit1 = hits[i];
        treePoint &hit2 = hits[i+1];
        double l1(0),l2(0);
        ntuple2->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),l1,l2,1.0); //true combination
        double phi2 = 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r2
        double theta2 = r.Rndm()*M_PI;
        double z2 = hit2.rz() * cos(theta2);
        ntuple2->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),phi2,z2,l1,l2,0.0); // wrong combination
    }
}

void combine3(xParticle &p)
{
    vector<treePoint> hits;
    transform(p,hits,MCHITS);
    
    // Sort the hits according to distance from origin
    sort(hits.begin(),hits.end(),Point::sortRad);
    
    // Combine 3 hits
    int nhits = (int)hits.size();
    if (nhits < 3) return;
    
    treePoint &hit1 = hits[0];
    treePoint &hit2 = hits[1];
    double l1(0),l2(0),l3(0);
    int istart = 2;
    float d = hit1.distance(hit2); // Check for double hits
    if (d<TWINDIST) {
        hit2 = hits[2];
        istart = 3;
    }
    for (int i=istart; i<nhits; i++)    {
       treePoint &hit3 = hits[i];
        ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),hit3.phi(),hit3.z(),l1,l2,l3,1.0); //true combination
        double phi3 = 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r3
        double theta3 = r.Rndm()*M_PI;
        double z3 = hit3.rz() * cos(theta3);
        ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),phi3,z3,l1,l2,l3,0.0); // wrong combination
    }
    
}
    
void combine3serial(xParticle &p)
{
    vector<treePoint> hits;
    transform(p,hits,MCHITS);

    // Sort the hits according to distance from origin
    sort(hits.begin(),hits.end(),Point::sortRad);

    // Combine 3 hits
    int nhits = (int)hits.size();
    if (nhits < 3) return;
    double l1(0),l2(0),l3(0);
    for (int i=0; i<nhits-2; i++)    {
        treePoint &hit1 = hits[i];
        treePoint &hit2 = hits[i+1];
        treePoint &hit3 = hits[i+2];
        ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),hit3.phi(),hit3.z(),l1,l2,l3,1.0); //true combination
        double phi3 = hit3.phi() - DELTAPHI * 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r3
        double theta3 = hit3.theta() - DELTATHE * 2.*(0.5-r.Rndm())*M_PI;
        double z3 = hit3.rz() * cos(theta3);
        ntuple3->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),hit3.rz(),phi3,z3,hit1.rz(),hit2.rz(),hit3.rz(),0.0); // wrong combination
    }
}

void combine2(xParticle &p1, xParticle &p2)
{
    vector<treePoint> hits1, hits2;
    transform(p1,hits1,MCHITS);
    transform(p2,hits2,MCHITS);
    
    // Sort the hits according to distance from origin
    sort(hits1.begin(),hits1.end(),Point::sortRad);
    sort(hits2.begin(),hits2.end(),Point::sortRad);
    
    // Combine 2 hits of the same track
    if (hits1.size() < 2) return;
    double l1(0),l2(0);
    for (auto it=hits1.begin(); it!=hits1.end()-1; it++)    {
        treePoint hit1 = *it;
        treePoint hit2 = *(it+1);
        ntuple2->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),l1,l2,1.0); //true combination
       it->settruth(1);
    }
    
    if (hits2.size() < 2) return;
    for (auto it=hits2.begin(); it!=hits2.end()-1; it++)    {
        treePoint hit1 = *it;
        treePoint hit2 = *(it+1);
        ntuple2->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),l1,l2,1.0); //true combination
        it->settruth(2);
    }
    
    hits1.insert(hits1.end(),hits2.begin(),hits2.end());
    sort(hits1.begin(),hits1.end(),Point::sortRad);
    
    // Combine 2 hits of mixed tracks
    for (auto it1=hits1.begin(); it1!=hits1.end()-1; it1++)    {
        treePoint hit1 = *it1;
        for (auto it2=it1; it2!=hits1.end(); it2++)    {
            treePoint hit2 = *it2;
            if (hit1.truth()==hit2.truth()) continue;
            double d = hit1.distance(hit2);
            if (d > DISTANCE*hit1.r()) continue; // Consider hits within a DISTANCE*r
            ntuple2->Fill(hit1.rz(),hit1.phi(),hit1.z(),hit2.rz(),hit2.phi(),hit2.z(),l1,l2,0.0); // wrong combination
        }
    }
}

void combine3(xParticle &p1, xParticle &p2)
{
    vector<treePoint> hits1, hits2;
    vector<treePoint>::iterator itb1,ite1,itb2,ite2;
    transform(p1,hits1,MCHITS);
    transform(p2,hits2,MCHITS);

    // Sort the hits according to distance from origin
    sort(hits1.begin(),hits1.end(),Point::sortRad);
    sort(hits2.begin(),hits2.end(),Point::sortRad);

    // Combine 3 hits
    if (hits1.size() < 3) return;
    if (hits2.size() < 3) return;
    
    itb1=hits1.begin();
    ite1=hits1.end();
    itb2=hits2.begin();
    ite2=hits2.end();
    
    double l1(0),l2(0),l3(0);
    for (auto it1=itb1, it2=itb2; it1!=ite1-4&&it2!=ite2-4; it1++, it2++)    {
        treePoint &hit11 = *it1; // Three consecutive hits of track1
        treePoint &hit12 = *(it1+1);
        treePoint &hit13 = *(it1+2);
        treePoint &hit21 = *it2; // Three consecutive hits of track2
        treePoint &hit22 = *(it2+1);
        treePoint &hit23 = *(it2+2);
        
        ntuple3->Fill(hit11.rz(),hit11.phi(),hit11.z(),hit12.rz(),hit12.phi(),hit12.z(),hit13.rz(),hit13.phi(),hit13.z(),l1,l2,l3,1.0); //true combination
        ntuple3->Fill(hit21.rz(),hit21.phi(),hit21.z(),hit22.rz(),hit22.phi(),hit22.z(),hit23.rz(),hit23.phi(),hit23.z(),l1,l2,l3,1.0); //true combination

        // Consider hits within a DISTANCE*r
        if (hit11.distance(hit21)<DISTANCE*hit11.r()) ntuple3->Fill(hit11.rz(),hit11.phi(),hit11.z(),hit12.rz(),hit12.phi(),hit12.z(),hit21.rz(),hit21.phi(),hit21.z(),l1,l2,l3,0.0); // wrong combination
        if (hit12.distance(hit21)<DISTANCE*hit12.r()) ntuple3->Fill(hit11.rz(),hit11.phi(),hit11.z(),hit12.rz(),hit12.phi(),hit12.z(),hit21.rz(),hit21.phi(),hit21.z(),l1,l2,l3,0.0); // wrong combination
        if (hit11.distance(hit22)<DISTANCE*hit11.r()) ntuple3->Fill(hit11.rz(),hit11.phi(),hit11.z(),hit12.rz(),hit12.phi(),hit12.z(),hit22.rz(),hit22.phi(),hit22.z(),l1,l2,l3,0.0); // wrong combination
        if (hit12.distance(hit22)<DISTANCE*hit12.r()) ntuple3->Fill(hit11.rz(),hit11.phi(),hit11.z(),hit12.rz(),hit12.phi(),hit12.z(),hit22.rz(),hit22.phi(),hit22.z(),l1,l2,l3,0.0); // wrong combination
        if (hit11.distance(hit23)<DISTANCE*hit11.r()) ntuple3->Fill(hit11.rz(),hit11.phi(),hit11.z(),hit12.rz(),hit12.phi(),hit12.z(),hit23.rz(),hit23.phi(),hit23.z(),l1,l2,l3,0.0); // wrong combination
        if (hit12.distance(hit23)<DISTANCE*hit12.r()) ntuple3->Fill(hit11.rz(),hit11.phi(),hit11.z(),hit12.rz(),hit12.phi(),hit12.z(),hit23.rz(),hit23.phi(),hit23.z(),l1,l2,l3,0.0); // wrong combination
        if (hit21.distance(hit11)<DISTANCE*hit21.r()) ntuple3->Fill(hit21.rz(),hit21.phi(),hit21.z(),hit22.rz(),hit22.phi(),hit22.z(),hit11.rz(),hit11.phi(),hit11.z(),l1,l2,l3,0.0); // wrong combination
        if (hit22.distance(hit11)<DISTANCE*hit22.r()) ntuple3->Fill(hit21.rz(),hit21.phi(),hit21.z(),hit22.rz(),hit22.phi(),hit22.z(),hit11.rz(),hit11.phi(),hit11.z(),l1,l2,l3,0.0); // wrong combination
        if (hit21.distance(hit12)<DISTANCE*hit21.r()) ntuple3->Fill(hit21.rz(),hit21.phi(),hit21.z(),hit22.rz(),hit22.phi(),hit22.z(),hit12.rz(),hit12.phi(),hit12.z(),l1,l2,l3,0.0); // wrong combination
        if (hit22.distance(hit12)<DISTANCE*hit22.r()) ntuple3->Fill(hit21.rz(),hit21.phi(),hit21.z(),hit22.rz(),hit22.phi(),hit22.z(),hit12.rz(),hit12.phi(),hit12.z(),l1,l2,l3,0.0); // wrong combination
        if (hit21.distance(hit13)<DISTANCE*hit21.r()) ntuple3->Fill(hit21.rz(),hit21.phi(),hit21.z(),hit22.rz(),hit22.phi(),hit22.z(),hit13.rz(),hit13.phi(),hit13.z(),l1,l2,l3,0.0); // wrong combination
        if (hit22.distance(hit13)<DISTANCE*hit22.r()) ntuple3->Fill(hit21.rz(),hit21.phi(),hit21.z(),hit22.rz(),hit22.phi(),hit22.z(),hit13.rz(),hit13.phi(),hit13.z(),l1,l2,l3,0.0); // wrong combination
    }
    
}

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

void GenerateTrack(std::vector<treePoint> &points, int np, double delta, double radius, double phi, double gamma, double error) {
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
        treePoint p(X,Y,Z,n++,-1);
        points.push_back(p);
    }
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

bool checkLayers( const int *baseLayers, int nBaseLayers )
{
    int icross=0;
    for( int il=0; il<nBaseLayers; il++ ){
        if( layerNHits[baseLayers[il]]>0 ) icross++;
    }
    return ( icross==nBaseLayers );
}


void readEvent( const char *directory, int event, bool loadMC )
{
    TString filePrefix;
    filePrefix.Form("%sevent%09d",directory,event);
    bool mcFlag;
    
    Geo::init();
    
    mcFlag = loadMC;
    
    mHits.clear();
    //mHitRecInfo.clear();
    //mTracks.clear();
    mHitsMC.clear();
    mParticles.clear();
    mHits.reserve(250000);
    
    // ===== load hits
    {
        TString fname = filePrefix+"-hits.csv";
        ifstream in(fname.Data());
        if( !in.is_open() ){
            cout<<"Event "<<event<<" does not exist!!"<<endl;
            exit(0);
        }
        char tmpLine[256];
        in.getline(tmpLine,256);
        cout<<tmpLine<<endl;
        while (1) {
            double h[7]; //file line: id:x:y:z:volume:layer:module
            if( !readLine(in,h,7) ) break;
            if( h[0]-1 != mHits.size() ){
                cout<<"Hit index is wrong: "<<h[0]<<endl;
                exit(0);
            }
            Hit hit;
            hit.hitID = h[0];
            hit.x = h[1]; // [mm]
            hit.y = h[2]; // [mm]
            hit.z = h[3]; // [mm]
            hit.r = sqrt( hit.x*hit.x + hit.y*hit.y );
            hit.phi = atan2( hit.y, hit.x );
            hit.module = h[6];
            hit.layer = h[5];
            hit.volume = h[4];
            int vol;
            int lay = ( (int) h[5])/2 -1;
            switch( (int) h[4] )
            {
                case  8: vol = 0; break;
                case  7: vol = 1; hit.layer = 6-hit.layer; break;
                case  9: vol = 2; break;
                case 13: vol = 3; break;
                case 12: vol = 4; hit.layer = 5-hit.layer; break;
                case 14: vol = 5; break;
                case 17: vol = 6; break;
                case 16: vol = 7; hit.layer = 5-hit.layer; break;
                case 18: vol = 8; break;
                default:
                    cout<<"Unknown detector volume: "<< (int) h[4] << endl;
                    exit(0);
            };
            hit.layer = h[5];
            
            if( lay<0 || lay>= Geo::volumes[vol].nLayers ){
                cout<<"Unknown detector layer: "<<hit.layer<<endl;
                exit(0);
            }
            hit.layerID = Geo::volumes[vol].layerIDs[lay];
            Layer1 &layer = Geo::layers[hit.layerID];
            if( Geo::volumes[vol].type==0 ){
                hit.t = hit.z / hit.r * layer.r;
            } else {
                hit.t = hit.r / hit.z * layer.z ;
            }
            hit.isUsed=0;
            hit.trackID=0;
            hit.BzFwd = Geo::OriginBzkG;
            hit.BzMid = Geo::OriginBzkG;
            hit.BzBck = Geo::OriginBzkG;
            hit.cells = 0;
            hit.values = 0.0;
            mHits.push_back(hit);
        }
        cout<<" loaded "<<mHits.size()<<" hits "<<endl;
        in.close();
        
    } // load hits
    
    // ===== load cells
    {
        TString fname = filePrefix+"-cells.csv";
        ifstream in(fname.Data());
        if( !in.is_open() ){
            cout<<"Event "<<event<<" does not exist!!"<<endl;
            exit(0);
        }
        char tmpLine[256];
        in.getline(tmpLine,256);
        cout<<tmpLine<<endl;
        int n = 0;
        while (1) {
            double h[4]; //file line: hit_id:ch0:ch1:value
            if( !readLine(in,h,4) ) break;
            int ih = h[0];
            mHits[ih].cells++;
            mHits[ih].values += h[3];
            n++;
        }
        cout<<" loaded "<<n<<" cells "<<endl;
        in.close();
        
    } // load cells
    
    // create particle ID->index map
    std::map<unsigned long long int,int> partIDmap;
    
    // ========= load particles with reindexing
    {
        mParticles.clear();
        
        TString fname = filePrefix+"-particles.csv";
        ifstream in(fname.Data());
        if( !in.is_open() ){
            cout<<"Particle file for event "<<event<<" does not exist!!"<<endl;
            exit(0);
        }
        char tmpLine[512];
        in.getline(tmpLine,512);
        cout<<tmpLine<<endl;
        while (1) {
            double f[10]; //  particle_id,vx,vy,vz,px,py,pz,q,nhits
            if( !readLine(in,f,10) ) break;
            //cout << f[0] << endl;
            int nhits = (int) f[9];
            if( nhits==0 ) continue; // no hits belong to the particle
            xParticle p( nhits );
            p.id = f[0];
            p.x = f[2]; // [mm]
            p.y = f[3]; // [mm]
            p.z = f[4]; // [mm]
            p.r = sqrt(p.x*p.x+p.y*p.y);
            p.px = f[5];
            p.py = f[6];
            p.pz = f[7];
            p.q  = f[8];
            p.xl = 0;
            p.yl = 0;
            p.zl = 0;
            p.rl = 0;
            p.pt = sqrt(f[5]*f[5] + f[6]*f[6]);
            p.p  = sqrt(f[5]*f[5] + f[6]*f[6] + f[7]*f[7]);
            p.prim = fabs(p.z)<1.2 && p.r<0.05;
            p.w = 0;
            partIDmap[ (unsigned long long int) f[0] ] = (int)mParticles.size();
            mParticles.push_back(p);
        }
        cout <<" loaded "<<mParticles.size() <<" particles in event "<<event<<endl;
        in.close();
    } // particles
    
    
    { // ============  read  mc truth
        
        mHitsMC.clear();
        mHitsMC.reserve(150000);
        
        TString fname = filePrefix+"-truth.csv";
        ifstream in(fname.Data());
        if( !in.is_open() ){
            cout<<"Truth file for event "<<event<<" does not exist!!"<<endl;
            exit(0);
        }
        
        char tmpLine[512];
        in.getline(tmpLine,512);
        cout<<tmpLine<<endl;
        
        while (1) {
            double mc[9]; //mc: hit_id,particle_id,tx,ty,tz,tpx,tpy,tpz,weight
            if( !readLine(in,mc,9) ) break;
            if( mc[0]-1 != mHitsMC.size() ){
                cout<<"MC hit index is wrong: "<<mc[0]<< " " << mHitsMC.size() << endl;
                exit(0);
            }
            HitMC hitmc;
            Hit &hit = mHits[mHitsMC.size()];
            
            hitmc.hitID = (int)mHitsMC.size();
            hitmc.id = mc[1];
            hit.id = mc[1];
            hitmc.x = mc[2];
            hitmc.y = mc[3];
            hitmc.z = mc[4];
            hitmc.partID = -1;
            hitmc.w = mc[8]; // weight
            hitmc.px = mc[5];
            hitmc.py = mc[6];
            hitmc.pz = mc[7];
            hitmc.pt = sqrt(hitmc.px*hitmc.px + hitmc.py*hitmc.py); // pt
            hitmc.dOrigin2 = -1.; // squared distance to particle origin
            hitmc.p = sqrt( hitmc.px*hitmc.px + hitmc.py*hitmc.py + hitmc.pz*hitmc.pz ); // momentum  for sorting
            hitmc.q = -1000;
            hitmc.r = sqrt( hitmc.x*hitmc.x + hitmc.y*hitmc.y );
            hitmc.phi = atan2( hitmc.y, hitmc.x );
            if( Geo::volumes[hit.volume].type==0 ) hitmc.t = hitmc.z;
            else hitmc.t = hitmc.r;
            
            
            // find mapped particle id
            long unsigned int id = mc[1];
            int newID = 0;
            if( id==0 ){ // hit is not associated to any particle
                mHitsMC.push_back(hitmc);
                continue;
            }
            std::map<unsigned long long int, int>::iterator it = partIDmap.find(id);
            if( it==partIDmap.end() ){
                cout<<"Particle ID not found in map!!!"<<endl;
                cout<<"ID= "<<id<<" hit "<<id<<" iterator at ID "<<it->first<<endl;
                exit(0);
            }
            newID = it->second;
            if( newID < 0 || newID>= mParticles.size() ){
                cout<<"Mapped particle ID is wrong!!!"<<endl;
                cout<<"ID= "<<id<<" new ID "<<newID<<endl;
                exit(0);
            }
            hitmc.partID = newID;
            xParticle &p = mParticles[newID];
            
            double dx = hitmc.x - p.x;
            double dy = hitmc.y - p.y;
            double dz = hitmc.z - p.z;
            hitmc.dOrigin2 = dx*dx + dy*dy + dz*dz;
            
            p.w += hitmc.w;
            hitmc.q = p.q;
            
            p.hits.push_back( hitmc.hitID );
            mHitsMC.push_back(hitmc);
        }
        cout << " read "<<mHitsMC.size() << " mc hits for event "<<event<<endl;
        in.close();
        if( mHitsMC.size() != mHits.size() ){
            cout<<"number of MC hits is not equal to the number of hits"<<endl;
            exit(0);
        }
    } // read mc hits
    
    
    // sort particle hits
    
    for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){
        xParticle &p = mParticles[ipart];
        if( p.hits.size()<1 ) continue;
        std::vector<HitMC> v;
        for( int i=0; i<p.hits.size(); i++ ){
            v.push_back( mHitsMC[p.hits[i]] );
        }
        std::sort(v.begin(), v.end() );
        for( int i=0; i<p.hits.size(); i++ ){
            p.hits[i] = v[i].hitID;
        }
        p.hitClusterIds.resize(p.hits.size());
        int il=-1;
        int lastLid=-1;
        for( int i=0; i<p.hits.size(); i++ ){
            if( mHits[p.hits[i]].layerID != lastLid ){
                il++;
                lastLid = mHits[p.hits[i]].layerID;
            }
            p.hitClusterIds[i] = il;
        }
        
        {
            int ihlast = p.hits[p.hits.size()-1];
            HitMC &mc = mHitsMC[ihlast];
            p.xl = mc.x;
            p.yl = mc.y;
            p.zl = mc.z;
            p.rl = sqrt(mc.x*mc.x+mc.y*mc.y);
        }
    }
    
    for( unsigned int ipart=0; ipart<mParticles.size(); ipart++ ){
        xParticle &p = mParticles[ipart];
        for( int i=0; i<Geo::NLayers; i++ ) layerNHits[i]=0;
        for( int i=0; i<p.hits.size(); i++ ){
            layerNHits[mHits[p.hits[i]].layerID]++;
        }
        p.nLayers = 0;
        for( int i=0; i<Geo::NLayers; i++ ){
            if( layerNHits[i]>0 ) p.nLayers ++;
        }
        
        p.baseV0 = 0;
        p.baseV1 = 0;
        p.baseV2 = 0;
        p.baseV  = 0;
        if( p.prim ){
            int baseV0[2] ={ Geo::getLayerID(0,0), Geo::getLayerID(0,1) };
            int baseV1[2] ={ Geo::getLayerID(1,0), Geo::getLayerID(1,1) };
            int baseV2[2] ={ Geo::getLayerID(2,0), Geo::getLayerID(2,1) };
            p.baseV0 = checkLayers(baseV0,2);
            p.baseV1 = checkLayers(baseV1,2);
            p.baseV2 = checkLayers(baseV2,2);
            p.baseV  =  ( p.baseV0 || p.baseV1 || p.baseV2 );
        }
        
        int baseA0[3] ={ Geo::getLayerID(0,0), Geo::getLayerID(0,1), Geo::getLayerID(0,2) };
        p.baseA0 = checkLayers(baseA0,3);
        
        int baseA1[3] ={ Geo::getLayerID(1,0), Geo::getLayerID(1,1), Geo::getLayerID(1,2) };
        p.baseA1 = checkLayers(baseA1,3);
        
        int baseA2[3] ={ Geo::getLayerID(2,2), Geo::getLayerID(2,3), Geo::getLayerID(2,4) };
        p.baseA2 = checkLayers(baseA2,3);
        
        int baseA3[3] ={ Geo::getLayerID(3,0), Geo::getLayerID(3,1), Geo::getLayerID(3,2) };
        p.baseA3 = checkLayers(baseA3,3);
        
    } // particles
    
}
