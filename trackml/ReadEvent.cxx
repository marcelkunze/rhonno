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
#include "RhoNNO/TXMLP.h"

using namespace std;

std::vector<Hit> mHits;
std::vector<HitMC> mHitsMC;
std::vector<Particle> mParticles;
std::vector<int> tracks[200000];

int layerNHits[Geo::NLayers];

#define NETFILE "/Users/marcel/workspace/rhonno/trackml/NNO0036.TXMLP"
//#define NETFILE "/Users/marcel/workspace/rhonno/Networks/NNO0100.TXMLP"
#define THRESHOLD 50
#define MAXHITS 1000
#define MAXPARTICLES 5000

double* Recall(double *invec);
int bestMatchingHit(size_t nhits, int **m, int row);
int bestHitPair(size_t nhits, int **m, int &row, int &col);
int findTracks(int nhits, float *x, float *y, float *z, float *v, int* labels);

TRandom r;
TNtuple *ntuple;

void print(vector<int> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
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
            hit.x = h[1]*0.1; // convert to [cm]
            hit.y = h[2]*0.1; // convert to [cm]
            hit.z = h[3]*0.1; // convert to [cm]
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
            p.x = f[2]*0.1; // [cm]
            p.y = f[3]*0.1; // [cm]
            p.z = f[4]*0.1; // [cm]
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
            hitmc.x = mc[2]*0.1; // convert to [cm]
            hitmc.y = mc[3]*0.1; // convert to [cm]
            hitmc.z = mc[4]*0.1; // convert to [cm]
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
        for(int j=0; j<nhits2; j++)    {
            Hit &hit2 = mHits[p2.hits[j]];
            //printf("%i %i: x1=%8f, y1=%8f, z1=%8f x2=%8f, y2=%8f, z2=%8f, v=%8f, l=%8f, m=%8f, t=%8f \n",i,j,hit1.x,hit1.y,hit1.z,hit2.x,hit2.y,hit2.z,(double)hit1.volume,(double)hit1.layer,(double)hit1.module,truth);
            Double_t cos = sqrt(hit1.x*hit2.x+hit1.y*hit2.y+hit1.z+hit2.z);
            cos /= sqrt(hit1.x*hit1.x+hit1.y*hit1.y+hit1.z+hit1.z);
            cos /= sqrt(hit2.x*hit2.x+hit2.y*hit2.y+hit2.z+hit2.z);
            if (i!=j) ntuple->Fill(hit1.x,hit1.y,hit1.z,hit2.x,hit2.y,hit2.z,cos,hit1.values,hit2.values,truth,p1.p);
        }
    }
}

int main()
{
    bool analyseTruth = true;
    
    const int nEvents = 1;
    const int firstEvent=21130;
    //const int firstEvent=1000;
    TString dir = "/Users/marcel/workspace/train_sample/";
    
    const long nParticles = MAXPARTICLES; // Number of particles to extract to ROOT file
    auto f = TFile::Open("tracks.root","RECREATE");
    ntuple = new TNtuple("tracks","training data","x1:y1:z1:x2:y2:z2:cos:v1:v2:truth:p");
    
    long int currentID=1;
    
    for( int event = firstEvent; event<firstEvent+nEvents; event++){
        cout<<"read event "<<event<<endl;
        readEvent( dir.Data(),  event, analyseTruth );
        
         for( int ih=0; ih<mHitsMC.size(); ih++ ){
            HitMC &h = mHitsMC[ih];
            if (ih<10) h.Print();
        }
        
        for( int ih=0; ih<mHits.size(); ih++ ){
            Hit &h = mHits[ih];
            if (ih<10) h.Print();
        }
        
        cout << "Particles:" << mParticles.size() << endl;
        for( int ip=0; ip<nParticles; ip++ ){
            int n1 = r.Rndm() * mParticles.size();
            Particle &p1 = mParticles[n1];
            int n2 = r.Rndm() * mParticles.size();
            cout << "Combine " << n1 << " " << n2 << endl;
            Particle &p2 = mParticles[n2];
            combine(p1,p1,1.0); // wright pairs
            combine(p2,p2,1.0); // wright pairs
            combine(p1,p2,0.0); // wrong pairs
            combine(p2,p1,0.0); // wrong pairs
        }
        
    }
    f->Write();
    
    size_t nhits = MAXHITS; //mHits.size();
    float x[nhits],y[nhits],z[nhits], v[nhits];
    int labels[nhits];
    int nt;
    for (int i=0; i<nhits; i++) {
        x[i] = mHits[i].x;
        y[i] = mHits[i].y;
        z[i] = mHits[i].z;
        v[i] = mHits[i].values;
    }

    cout << "Find tracks:" << endl;
    nt = findTracks(nhits,x,y,z,v,labels);
    
    cout << "Number of tracks:" << nt << endl;
    for(int i=0; i<nt; i++) {
        cout << "Track " << i << ":";
        print(tracks[i]);
        cout << endl;
    }
    
    cout << "Labels: ";
    for (int i=0;i<nhits;i++) cout << labels[i] << " ";
    cout << endl;
    
}

// Assign track labels to hits (x,y,z)
// The hit pair quality is assessed by the neural network
// The quality is noted in the hit pair matrix m[nhits][nhits]
int findTracks(int nhits, float *x, float *y, float *z, float *v, int* labels)
{
    for (int i=0;i<nhits;i++) labels[i] = -1; // Preset with no match
    
    double in1[7], in2[7], *out;
    
    // Allocate a nhits*nhits hit pair matrix as one continuous memory block
    int** m = new int*[nhits];
    if (nhits)
    {
        m[0] = new int[nhits * nhits];
        for (int i = 1; i < nhits; ++i)
            m[i] = m[0] + i * nhits;
    }
    
    for(int i=0; i<nhits-1; i++)    {
        in1[0] = x[i];
        in1[1] = y[i];
        in1[2] = z[i];
        in2[3] = x[i];
        in2[4] = y[i];
        in2[5] = z[i];
        for(int j=i+1; j<nhits; j++)    {
            in1[3] = x[j];
            in1[4] = y[j];
            in1[5] = z[j];
            in2[0] = x[j];
            in2[1] = y[j];
            in2[2] = z[j];
            double cos = sqrt(in1[0]*in2[0] + in1[1]*in2[1] + in1[2]*in2[2]);
            cos /= sqrt(in1[0]*in1[0] + in1[1]*in1[1] + in1[2]*in1[2]);
            cos /= sqrt(in2[0]*in2[0] + in2[1]*in2[1] + in2[2]*in2[2]);
            in1[6] = cos;
            in2[6] = cos;
            m[i][j] = (int) 100. * Recall(in1)[0]; // Recall the hit pair matching quality
            m[j][i] = (int) 100. * Recall(in2)[0];
            if (m[i][j]< THRESHOLD) m[i][j] = 0; // Apply a cut on the quality
            if (m[j][i]< THRESHOLD) m[j][i] = 0;
        }
    }
    
    // Analyze the hit pair matrix
    // Sort out the tracks by following the network connections and fill the corresponding track hits into containers
    int ntracks = 0;
    
    int row, col;
    int seed = bestHitPair(nhits, m, row, col); // Look for seed
    while (seed > -1) {
        cout << "Best hit pair: (" << row << "," << col << ")" << endl;
        while (seed>-1) {
            labels[seed] = ntracks;
            tracks[ntracks].push_back(seed);
            seed = bestMatchingHit(nhits, m, seed);
        }
        ntracks++;
        seed = bestHitPair(nhits, m, row, col); // Look for new seed
    }
    
    if (nhits) delete [] m[0]; // Clean up the memory
    delete [] m;
    
    return ntracks;
}

// Examine the hit pair matrix to find the best matching hit
// The corresponding rows and columns are crossed out (-1)
int bestMatchingHit(size_t nhits, int **m, int row)
{
    int imax = -1;
    static int imaxold = -1;
    
    cout << "Seed hit: " << row << endl;
    Double_t max = 0.0;
    for(size_t i=0; i<nhits; i++)    {
        if (i!=row && m[row][i] > max) {
            max = m[row][i]; // Best matching hit
            imax = (int) i;
        }
    }
    cout << "  Best matching hit: " << imax << endl;
    for(size_t i=0; i<nhits; i++) m[row][i] = -1.0;
    for(size_t i=0; i<nhits; i++) m[i][row] = -1.0;
    if (imax == imaxold) {
        return -1;
    }
    imaxold = imax;
    return imax;
}

// Find the best matching hit pair in the matrix
int bestHitPair(size_t nhits, int **m, int &row, int &col)
{
    col = -1;
    row = -1;
    
    Double_t max = -1.0;
    for(int i=0; i<nhits; i++)    {
        for(int j=0; j<nhits; j++)    {
            if (m[i][j] > max) {
                max = m[i][j]; // Best hit pair
                row = i;
                col = j;
            }
        }
    }
    return row;
}

// Recall function on normalised network input
double* Recall(double *invec)
{
    static TXMLP net(NETFILE);
    float x[7],y[1];
    x[0]     = 7.75563     *    invec[0];    // x1
    x[1]     = 7.16209     *    invec[1];    // y1
    x[2]     = 1.00536     *    invec[2];    // z1
    x[3]     = 7.75563     *    invec[3];    // x2
    x[4]     = 7.16209     *    invec[4];    // y2
    x[5]     = 1.00536     *    invec[5];    // z2
    x[6]     = 0.000245409 *    invec[6];    // dot
    return net.Recallstep(x,y);
}
