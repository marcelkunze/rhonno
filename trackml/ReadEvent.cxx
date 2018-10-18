// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#define TRACKML
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

#include "DataStructures.h"
#include "Geo.h"
#include "XMLP.h"
#include "Tracker.h"

using namespace std;

std::vector<Hit> mHits;
std::vector<HitMC> mHitsMC;
std::vector<Particle> mParticles;

int layerNHits[Geo::NLayers];

#define NEVENTS 1
#define MAXHITS 150000
#define MAXPARTICLES 100
#define MCHITS true

TRandom r;
TNtuple *ntuple,*ntuple3;

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
            hit.x = 0.001 * h[1]; // [m]
            hit.y = 0.001 * h[2]; // [m]
            hit.z = 0.001 * h[3]; // [m]
            hit.r = sqrt( hit.x*hit.x + hit.y*hit.y );
            hit.phi = atan2( hit.y, hit.x );
            hit.module = h[6];
            hit.layer = ( (int) h[5])/2 -1;
            switch( (int) h[4] )
            {
                case  8: hit.volume = 0; break;
                case  7: hit.volume = 1; hit.layer = 6-hit.layer; break;
                case  9: hit.volume = 2; break;
                case 13: hit.volume = 3; break;
                case 12: hit.volume = 4; hit.layer = 5-hit.layer; break;
                case 14: hit.volume = 5; break;
                case 17: hit.volume = 6; break;
                case 16: hit.volume = 7; hit.layer = 5-hit.layer; break;
                case 18: hit.volume = 8; break;
                default:
                    cout<<"Unknown detector volume: "<< (int) h[4] << endl;
                    exit(0);
            };
            
            if( hit.layer<0 || hit.layer>= Geo::volumes[hit.volume].nLayers ){
                cout<<"Unknown detector layer: "<<hit.layer<<endl;
                exit(0);
            }
            hit.layerID = Geo::volumes[hit.volume].layerIDs[hit.layer];
            Layer &layer = Geo::layers[hit.layerID];
            if( Geo::volumes[hit.volume].type==0 ){
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
            Particle p( nhits );
            p.x = 0.001 * f[2]; // [m]
            p.y = 0.001 * f[3]; // [m]
            p.z = 0.001 * f[4]; // [m]
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
            hitmc.x = 0.001 * mc[2]; // convert to [m]
            hitmc.y = 0.001 * mc[3]; // convert to [m]
            hitmc.z = 0.001 * mc[4]; // convert to [m]
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
            Particle &p = mParticles[newID];
            
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
        Particle &p = mParticles[ipart];
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
        Particle &p = mParticles[ipart];
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

void transform(Particle &particle, std::vector<Point> &points, bool mc=false) {
    int nhits = (int)particle.hits.size();
    if (!mc) {
        for (int i=0;i<nhits;i++) {
            Hit &h1 = mHits[particle.hits[i]];
            Point p(h1.x,h1.y,h1.z,h1.hitID,h1.trackID);
            points.push_back(p);
        }
    }
    else {
        for (int i=0;i<nhits;i++) {
            HitMC &h1 = mHitsMC[particle.hits[i]];
            Point p(h1.x,h1.y,h1.z,h1.hitID+1,-1);
            points.push_back(p);
        }
    }
}

void combine(Particle &p1, Particle &p2)
{
    int nhits1 = (int)p1.hits.size();
    int nhits2 = (int)p2.hits.size();
    
    vector<Point> hits1, hits2;
    transform(p1,hits1,MCHITS);
    transform(p2,hits2,MCHITS);
    
    for (int i=0; i<nhits1; i++)    {
        Hit &h1 = mHits[p1.hits[i]];
        Point &hit1 = hits1[i];
        for (int j=0; j<nhits1; j++)    {
            Hit &h2 = mHits[p1.hits[j]];
            Point &hit2 = hits1[j];
            if (i!=j) ntuple->Fill(hit1.r(),hit1.phi(),hit1.theta(),hit2.r(),hit2.phi(),hit2.theta(),h1.values,h2.values,1.0); //true combination
        }
        for (int j=0; j<nhits2; j++)    {
            Hit &h2 = mHits[p2.hits[j]];
            Point &hit2 = hits2[j];
            if (i!=j) ntuple->Fill(hit1.r(),hit1.phi(),hit1.theta(),hit2.r(),hit2.phi(),hit2.theta(),h1.values,h2.values,0.0); // wrong combination
        }
    }
}

void combine2(Particle &p)
{
    vector<Point> hits;
    transform(p,hits,MCHITS);
    
    // Sort the hits according to distance from origin
    sort(hits.begin(),hits.end(),Point::sortDist);

    int nhits = (int)hits.size();
    for (int i=0; i<nhits-1; i++)    {
        Hit &h1 = mHits[p.hits[i]];
        Point &hit1 = hits[i];
        Hit &h2 = mHits[p.hits[i+1]];
        Point &hit2 = hits[i+1];
        ntuple->Fill(hit1.r(),hit1.phi(),hit1.theta(),hit2.r(),hit2.phi(),hit2.theta(),h1.values,h2.values,1.0); //true combination
        float phi2 = 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r2
        float theta2 = r.Rndm()*M_PI;
        ntuple->Fill(hit1.r(),hit1.phi(),hit1.theta(),hit2.r(),phi2,theta2,h1.values,h2.values,0.0); // wrong combination
    }
}

void combine3(Particle &p)
{
    vector<Point> hits;
    transform(p,hits,MCHITS);

    // Sort the hits according to distance from origin
    sort(hits.begin(),hits.end(),Point::sortDist);

    // Combine 3 hits
    int nhits = (int)hits.size();
    for (int i=0; i<nhits-2; i++)    {
        Hit &h1 = mHits[p.hits[i]];
        Hit &h2 = mHits[p.hits[i+1]];
        Hit &h3 = mHits[p.hits[i+2]];
        Point &hit1 = hits[i];
        Point &hit2 = hits[i+1];
        Point &hit3 = hits[i+2];
        ntuple3->Fill(hit1.r(),hit1.phi(),hit1.theta(),hit2.r(),hit2.phi(),hit2.theta(),hit3.r(),hit3.phi(),hit3.theta(),h1.values,h2.values,h3.values,1.0); //true combination
        float phi3 = 2.*(0.5-r.Rndm())*M_PI; // Generate a random point on sphere with r3
        float theta3 = r.Rndm()*M_PI;
        ntuple3->Fill(hit1.r(),hit1.phi(),hit1.theta(),hit2.r(),hit2.phi(),hit2.theta(),hit3.r(),phi3,theta3,h1.values,h2.values,h3.values,0.0); // wrong combination
    }
    
}

int main()
{
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
        ntuple = new TNtuple("tracks","training data","r1:phi1:theta1:r2:phi2:theta2:v1:v2:truth");
        ntuple3 = new TNtuple("tracks3","training data","r1:phi1:theta1:r2:phi2:theta2:r3:phi3:theta3:v1:v2:v3:truth");
        
        filePrefix.Form("%sevent%09d",dir.Data(),event+100000000);
        TString hname = filePrefix+"-hits.csv";
        ofstream outhits(hname);
        cout << "Writing hits data to " << hname << endl;
        
        for (int ih=0; ih<mHitsMC.size(); ih++ ) {
            HitMC &h = mHitsMC[ih];
            if (ih<10) h.Print();
        }
        
        for (int ih=0; ih<mHits.size(); ih++ ) {
            Hit &h = mHits[ih];
            if (ih<10) h.Print();
        }
        
        std::vector<Point> hits;
        hits.reserve(MAXHITS);
        
        
        size_t nParticles = mParticles.size();
        if (nParticles > MAXPARTICLES) nParticles = MAXPARTICLES;
        cout << "Particles:" << nParticles << endl;
        
        int n = 0;
        for (int ip=0; ip<nParticles; ip++ ) {
            int jp = r.Rndm() * nParticles;
            while (ip == jp) jp = r.Rndm() * nParticles; // Do not combine the particle with itself
            Particle &p1 = mParticles[ip];
            Particle &p2 = mParticles[jp];
            combine2(p1); // Produce training data (2 hits)
            combine3(p1); // Produce training data (3 hits)

            // Print the hit ids
            transform(p1,hits,MCHITS);
            if (ip == 10) cout << endl << "..." << endl;
            if (ip<10 || ip>nParticles-10) {
                cout << "Track " << ip+1 << ": ";
                for (int j=0; j<p1.hits.size(); j++) {
                    Hit &h = mHits[p1.hits[j]];
                    cout << h.hitID << "(" << n++ << ") ";
                }
                cout << endl;
            }
        }
        
        f->Write();
        delete ntuple;
        delete ntuple3;
        
        size_t nhits = hits.size();
        if (nhits > MAXHITS) nhits = MAXHITS;
        cout << "Hits: " << nhits << endl;
        
        float x[nhits],y[nhits],z[nhits];
        int labels[nhits];
        int nt;
        
        for (int i=0; i<nhits; i++) {
            labels[i] = 0;
            Point &h = hits[i];
            x[i] = h.x();
            y[i] = h.y();
            z[i] = h.z();
        }
        
        cout << "Find tracks..." << endl;
        nt = Tracker::findTracks((int)nhits,x,y,z,labels);
        
        // Assemble the tracks from label information
        std::vector<Point> tracks[nt];
        for(int i=0; i<nt; i++) {
            int track = i+1;
            for (int j=0;j<nhits;j++) {
                if (track != labels[j]) continue;
                hits[j].setval(labels[j]);
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
                    cout << p.id() << " ";
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
        
        // Wite a hits file with hits sorted by tracks
        cout << endl << "Write hits file..." << endl;
        //outhits<<"hit_id,x,y,z,volume_id,layer_id,module_id"<<endl;
        int nh = 0;
        for (int ip=0; ip<mParticles.size(); ip++ ) {
            Particle &p1 = mParticles[ip];
            for (int j=0;j<p1.hits.size();j++) {
                if (nh++>MAXHITS) break;
                int index = p1.hits[j];
                Hit h = mHits[index];
                outhits<<100.*h.x<<"\t"<<100.*h.y<<"\t"<<100.*h.z<<endl; // Coordinates in cm
                //outhits<<h.volume<<","<<h.layer<<","<<h.module<<endl;
            }
        }
        
        cout << "Write tracks file..." << endl;
        for (int ip=0; ip<mParticles.size(); ip++ ) {
            Particle &p1 = mParticles[ip];
            outtracks<<event<<","<<ip<<": ";
            for (int j=0;j<p1.hits.size();j++) outtracks<<p1.hits[j]<<" ";
            outtracks << endl;
        }
        outtracks.close();
        
        cout << "Write submission file..." << endl;
        for( int ih=0; ih<nhits; ih++ ){
            out<<event<<","<<ih+1<<","<<labels[ih]<<endl;
        }
        out.close();
        
        // Initialize a 3D canvas and draw the hits and tracks
        if (DRAW) {
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
            for (vector<Point>::iterator it = hits.begin(); it != hits.end(); it++)    {
                static int i = 0;
                Point p=*it;
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
                vector<Point> track = tracks[i];
                int n = 0;
                TPolyLine3D *connector = new TPolyLine3D((int)track.size());
                for (vector<Point>::iterator it = track.begin(); it != track.end(); it++)    {
                    Point hit = *it;
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
    
}
