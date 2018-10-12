#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom.h"

#include <iostream>
#include <fstream>

#include "DataStructures.h"
#include "Geo.h"
#include "XMLP.h"

using namespace std;

std::vector<Hit> mHits;
std::vector<HitMC> mHitsMC;
std::vector<Particle> mParticles;
std::vector<int> tracks[200000];

int layerNHits[Geo::NLayers];

#define NEVENTS 1
#define MAXPARTICLES 200000

#define MAXHITS 1000
#define NETFILE "/Users/marcel/workspace/rhonno/trackml/NNO0200-6-25-15-1.TXMLP"
//#define NETFILE "/Users/marcel/workspace/rhonno/trackml/NNO0128-6-20-10-1.TXMLP"
//#define NETFILE "/Users/marcel/workspace/rhonno/Networks/NNO0100.TXMLP"
#define TRACKLET 3
#define THRESHOLD 90
#define DISTANCE 1.0
#define DELTAR   0.1
#define DELTAPHI 0.1
#define DELTATHETA 0.1

#define VERBOSE false

Double_t* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist);
int findTracks(int nhits, float *x, float *y, float *z, int* labels);

bool sortFunc( const vector<int>& p1,
              const vector<int>& p2 ) {
    return p1.size() > p2.size();
}

struct Point
{
    int id;             // Hit id of point
    int val;            // Group of point
    double x, y, z;     // Cartesian coordinate of point
    double r,phi,theta; // Spherical coordinates of point
    double distance;    // Distance from test point
};

// Calculate the circle center through 3 points
Point circleCenter(const Point &p1,const Point &p2,const Point &p3)
{
    Point center;
    center.x = 0.0;
    center.y = 0.0;
    center.z = 0.0;
    
    static double TOL = 0.0000001;
    double offset = pow(p2.x,2) + pow(p2.y,2);
    double bc =   ( pow(p1.x,2) + pow(p1.y,2) - offset )/2.0;
    double cd =   (offset - pow(p3.x, 2) - pow(p3.y, 2))/2.0;
    double det =  (p1.x - p2.x) * (p2.y - p3.y) - (p2.x - p3.x)* (p1.y - p2.y);
    
    if (abs(det) < TOL) { return center; }
    
    double idet = 1/det;
    
    center.x =  (bc * (p2.y - p3.y) - cd * (p1.y - p2.y)) * idet;
    center.y =  (cd * (p1.x - p2.x) - bc * (p2.x - p3.x)) * idet;
    
    return center;
}

// Calculate the circle radius through 3 points
double circleRadius(const Point &p1,const Point &p2,const Point &p3)
{
    Point center = circleCenter(p1,p2,p3);
    double radius = sqrt( pow(p2.x - center.x,2) + pow(p2.y-center.y,2));
    return radius;
}

// Used to sort an array of points by increasing
// order of distance from origin
bool sortDist(const Point &a,const Point &b)
{
    return (a.r < b.r);
}

// Used to sort an array of points by increasing
// order of distance from origin
bool sortId(const Point &a,const Point &b)
{
    return (a.id < b.id);
}

// Used to sort an array of points by increasing
// order of distance
bool comparison(const Point &a,const Point &b)
{
    return (a.distance < b.distance);
}

double distance(const Point &a,const Point &b)
{
    double d =  sqrt((a.x - b.x) * (a.x - b.x) +
                     (a.y - b.y) * (a.y - b.y) +
                     (a.z - b.z) * (a.z - b.z));
    return d;
}

// This function finds classification of point p using
// k nearest neighbour algorithm. It assumes only two
// groups and returns 0 if p belongs to group 0, else
// 1 (belongs to group 1).
int classifyAPoint(Point arr[], int n, int k, Point p)
{
    // Fill distances of all points from p
    for (int i = 0; i < n; i++) arr[i].distance = distance(arr[i],p);
    
    // Sort the Points by distance from p
    sort(arr, arr+n, comparison);
    
    // Now consider the first k elements and only
    // two groups
    int freq1 = 0;     // Frequency of group 0
    int freq2 = 0;     // Frequency of group 1
    for (int i = 0; i < k; i++)
    {
        if (arr[i].val == 0)
            freq1++;
        else if (arr[i].val == 1)
            freq2++;
    }
    
    return (freq1 > freq2 ? 0 : 1);
}


TRandom r;
TNtuple *ntuple;

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
            if( 1|| ( hit.volume==4 || hit.volume==5 || hit.volume==7 || hit.volume==8 ) ){
                hit.BzFwd = Geo::layers[hit.layerID].getFieldFwd(hit.phi, hit.t);
                hit.BzMid = Geo::layers[hit.layerID].getFieldMid(hit.phi, hit.t);
                hit.BzBck = Geo::layers[hit.layerID].getFieldBck(hit.phi, hit.t);
            }
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
            partIDmap[ (unsigned long long int) f[0] ] = mParticles.size();
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
            
            hitmc.hitID = mHitsMC.size();
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
        if( 0 && ipart==467 ){
            for( int i=0; i<p.hits.size(); i++ ){
                Hit &h = mHits[p.hits[i]];
                h.Print();
            }
            for( int i=0; i<p.hits.size(); i++ ){
                HitMC &mc = mHitsMC[p.hits[i]];
                mc.Print();
            }
        }
    }
    //exit(0);
    
    
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


void combine(Particle &p1, Particle &p2, Double_t truth)
{
    int nhits1 = p1.hits.size();
    int nhits2 = p2.hits.size();
    
    //if (nhits1 > 3) nhits1 = 3;
    //if (nhits2 > 3) nhits2 = 3;
    
    for(int i=0; i<nhits1; i++)    {
        Hit &hit1 = mHits[p1.hits[i]];
        float r1 = sqrt(hit1.x*hit1.x+hit1.y*hit1.y+hit1.z*hit1.z); // convert to m
        float phi1 = atan2(hit1.y,hit1.x);
        float theta1 = acos(hit1.z/r1);
        for(int j=0; j<nhits2; j++)    {
            if (nhits1 == nhits2 && i == j) continue; // Do not combine the hit with itself
            Hit &hit2 = mHits[p2.hits[j]];
            float r2 = sqrt(hit2.x*hit2.x+hit2.y*hit2.y+hit2.z*hit2.z); // convert to m
            float phi2 = atan2(hit2.y,hit2.x);
            float theta2 = acos(hit2.z/r2);
            //printf("%i %i: x1=%8f, y1=%8f, z1=%8f x2=%8f, y2=%8f, z2=%8f, v1=%8f, v2=%8f, v=%8f, l=%8f, m=%8f, t=%8f \n",i,j,hit1.x,hit1.y,hit1.z,hit2.x,hit2.y,hit2.z,hit1.values,hit2.values,(double)hit1.volume,(double)hit1.layer,(double)hit1.module,truth);
            //double dist = sqrt((hit1.x-hit2.x)*(hit1.x-hit2.x) + (hit1.y-hit2.y)*(hit1.y-hit2.y) + (hit1.z-hit2.z)*(hit1.z-hit2.z));
            if (i!=j) ntuple->Fill(hit1.x,hit1.y,hit1.z,hit2.x,hit2.y,hit2.z,r1,phi1,theta1,r2,phi2,theta2,hit1.values,hit2.values,truth);
        }
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
    
    long nParticles = 0;
    
    for (int event = firstEvent; event<firstEvent+nEvents; event++){
        cout<<"read event "<<event<<endl;
        readEvent( dir.Data(),  event, analyseTruth );
        
        TString filePrefix;
        filePrefix.Form("%sevent%09d",dir.Data(),event);
        TString fname = filePrefix+".root";
        auto f = TFile::Open(fname,"RECREATE");
        cout << "Writing training data to " << fname << endl;
        ntuple = new TNtuple("tracks","training data","x1:y1:z1:x2:y2:z2:r1:phi1:theta1:r2:phi2:theta2:v1:v2:truth");
        
        filePrefix.Form("%sevent%09d",dir.Data(),event+100000000);
        TString hname = filePrefix+"-hits.csv";
        ofstream outhits(hname);
        cout << "Writing hits data to " << hname << endl;
        
        for( int ih=0; ih<mHitsMC.size(); ih++ ){
            HitMC &h = mHitsMC[ih];
            if (ih<10) h.Print();
        }
        
        for( int ih=0; ih<mHits.size(); ih++ ){
            Hit &h = mHits[ih];
            if (ih<10) h.Print();
        }
        
        cout << "Particles:" << mParticles.size() << endl;
        for (int ip=0; ip<mParticles.size(); ip++ ) {
            if (nParticles++ > MAXPARTICLES) break;
            int jp = r.Rndm() * mParticles.size();
            if (ip == jp) continue; // Do not combine the particle with itself
            if (VERBOSE) cout << "Combine " << ip << " " << jp << endl;
            Particle &p1 = mParticles[ip];
            Particle &p2 = mParticles[jp];
            combine(p1,p1,1.0); // wright pairs
            combine(p1,p2,0.0); // wrong pairs
        }
        f->Write();
        delete ntuple;
        
        size_t nhits = mHits.size(); //mHits.size();
        if (nhits > MAXHITS) nhits = MAXHITS;
        cout << "Hits: " << nhits << endl;
        
        float x[MAXHITS],y[MAXHITS],z[MAXHITS];
        Point p[MAXHITS];
        int labels[MAXHITS];
        int nt;
        for (int i=0; i<nhits; i++) {
            x[i] = mHits[i].x;
            y[i] = mHits[i].y;
            z[i] = mHits[i].z;
            p[i].x = mHits[i].x;
            p[i].y = mHits[i].y;
            p[i].z = mHits[i].z;
        }
        
        cout << "Find tracks..." << endl;
        nt = findTracks(nhits,x,y,z,labels);
        
#define MAXTRACK 25
        cout << "Number of tracks:" << nt << endl;
        for(int i=0; i<nt; i++) {
            if (i<MAXTRACK || i>nt-MAXTRACK) {
                cout << "Track " << i+1 << ": ";
                print(tracks[i]);
            }
            if (i == MAXTRACK) cout << endl << "..." << endl;
        }
        
#define MAXLABEL 250
        cout << "Labels: ";
        for (int i=0;i<nhits;i++) {
            if (i<MAXLABEL || i>nhits-MAXLABEL) cout << labels[i] << " ";
            if (i == MAXLABEL) cout << endl << "..." << endl;
        }

        // Wite a hits file with hits sorted by tracks
        cout << "Write hits file..." << endl;
        //outhits<<"hit_id,x,y,z,volume_id,layer_id,module_id"<<endl;
        int nh = 0;
        for( int ip=0; ip<mParticles.size(); ip++ ){
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
        for( int ip=0; ip<mParticles.size(); ip++ ){
            Particle &p1 = mParticles[ip];
            outtracks<<event<<","<<ip<<": ";
            for (int j=0;j<p1.hits.size();j++) outtracks<<p1.hits[j]<<" ";
            outtracks << endl;
        }
        
        cout << "Write submission file..." << endl;
        for( int ih=0; ih<nhits; ih++ ){
            out<<event<<","<<ih+1<<","<<labels[ih]<<endl;
        }
    }
    
    out.close();
    outtracks.close();
}

// Assign track labels to hits (x,y,z)
// The hit pair quality is assessed by the neural network
// The quality is noted in the hit pair matrix m[nhits][nhits]
int findTracks(int nhits, float *x, float *y, float *z, int* labels)
{
    std::clock_t c_start = std::clock();
    
    Point *p = new Point[nhits];
    vector<Point> points;
    
    // Set up a cache for the point coordinates
    //cout << "Set up points cache..." << endl;
    for (int i=0;i<nhits;i++) {
        labels[i] = 0;
        p[i].id = i;
        p[i].val = -1; // Preset group with no match
        p[i].x = x[i]; // Cache the point coordinates
        p[i].y = y[i];
        p[i].z = z[i];
        p[i].r = sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
        p[i].phi = atan2(p[i].y,p[i].x);
        p[i].theta = acos(p[i].z/p[i].r);
        p[i].distance = 0.0;
        points.push_back(p[i]);
    }
    
    // Sort the hits according to distance from origin
    cout << "Sorting hits..." << endl;
    sort(points.begin(),points.end(),sortDist);
    
    // Search neighbouring hits, the neural network recall identifies the hit belonging to a tracklet
    cout << "Find tracklets..." << endl;
    int nd(0), nr(0), np(0), nt(0);
    Point vertex;
    vertex.x = 0;
    vertex.y = 0;
    vertex.z = 0;
    vector<vector<int>> tracklet;
    for (vector<Point>::iterator it1 = points.begin(); it1 != points.end(); ++it1) {
        Point p1 = *it1; // Seeding point
        vector<Point> pvec;
        // Conformal mapping of circle to straight line
        pvec.push_back(p1); //Note the seeding point in the first place
        double r1 = 0.0;
        if (VERBOSE) cout << endl << p1.id << "(0) ";
        for (vector<Point>::iterator it2 = it1+1; it2 != points.end(); ++it2) { //
            Point p2 = *it2;
            double dist = distance(p1,p2); // Only consider points in the neighborhood
            double r2 = 0.0;
            if (pvec.size()==2) r1 = circleRadius(pvec[0],pvec[1],p2);
            if (pvec.size()>=2)  r2 = circleRadius(pvec[0],pvec[1],p2);
            double deltar = abs(r1-r2);
            double deltaphi = abs(abs(p1.phi)-abs(p2.phi));
            double deltatheta = abs(p1.theta-p2.theta);
            //cout << endl << p1.id << " " << p2.id << " " << dist << " " << deltar << " " << deltaphi << " " << deltatheta << endl;
            if (dist > DISTANCE) {
                //if (VERBOSE) cout << p1.id << " " << p2.id << " dist:" << dist << endl;
                nd++;
                continue;
            }
            if (deltar > DELTAR) {
                //if (VERBOSE && pvec.size()>=2) cout << pvec[0].id << " " << pvec[1].id << " " << p2.id << " deltar:" << deltar << endl;
                nr++;
                continue;
            }
            if (deltaphi > DELTAPHI) {
                //if (VERBOSE) cout << p1.id << " " << p2.id << " deltaphi:" << deltaphi << endl;
                np++;
                continue;
            }
            if (deltatheta > DELTATHETA) {
                //if (VERBOSE) cout << p1.id << " " << p2.id << " deltatheta:" << deltatheta << endl;
                nt++;
                continue;
            }
            int recall1 = (int) 100. * Recall(p1.r,p1.phi,p1.theta,p2.r,p2.phi,p2.theta,dist)[0]; // Recall the hit pair matching quality
            int recall2 = (int) 100. * Recall(p2.r,p2.phi,p2.theta,p1.r,p1.phi,p1.theta,dist)[0]; // Recall the hit pair matching quality
            if (recall1 < THRESHOLD) recall1 = 0; // Apply a cut on the quality
            if (recall2 < THRESHOLD) recall2 = 0;
            int recall  = (recall1>recall2) ? recall1:recall2;
            if (recall>THRESHOLD) {
                pvec.push_back(p2); // Note the columns with a good combination
                if (VERBOSE) cout << p2.id << "(" << recall << ") ";
                points.erase(it2);  // Remove the corresponding point from the set
                *it2--;
                p1 = p2; // Note the assigned hit
                continue;
            }
        }
        points.erase(it1);  // Remove the corresponding point from the set
        *it1--;
        
        if (pvec.size() < TRACKLET) continue; // Perform a cut on tracklet size
        sort(pvec.begin(), pvec.end(), sortId); // Sort the hits acording to the Id
        vector<int> tmpvec;
        for (int ip=0;ip<pvec.size();ip++) tmpvec.push_back(pvec[ip].id); // Note the hit indices
        tracklet.push_back(tmpvec);
    }
    
    cout << endl << "Number of tracklets: " << tracklet.size() << endl;
    if (tracklet.size() == 0) exit(0);
    
    // Sort the tracklet vector according to the tracklet length
    
    sort(tracklet.begin(), tracklet.end(), sortFunc);
    
    if (VERBOSE) {
        cout << "Tracks:" << endl;
        for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    }
    
    int n = 0;
    for (int i=0;i<tracklet.size();i++) {
        tracks[i] = tracklet[i]; // Save the results
        for (int j=0;j<tracklet[i].size();j++) {
            int hit = tracklet[i][j];
            labels[hit] = i+1;
            p[hit].val = i+1;
            n++;
        }
    }
    
    cout << "Number of assigned points: " << n << endl;
    cout << "Distance <" << DISTANCE << ": " << nd <<endl;
    cout << "Radius   <" << DELTAR << ": " << nr <<endl;
    cout << "Phi      <" << DELTAPHI << ": " << np <<endl;
    cout << "Theta    <" << DELTATHETA << ": " << nt <<endl;
    
    delete [] p;
    
    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
    
    return (int) tracklet.size();
}

// Recall function on normalised network input
double* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist)
{
    static XMLP net(NETFILE);
    float x[7];
    x[0]     = x1;    // r1
    x[1]     = y1;    // phi1
    x[2]     = z1;    // theta1
    x[3]     = x2;    // r2
    x[4]     = y2;    // phi2
    x[5]     = z2;    // theta2
    x[6]     = dist;
    return net.Recallstep(x);
}
