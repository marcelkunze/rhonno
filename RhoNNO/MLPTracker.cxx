// Read HITS spatial data and try a XMLP network to segment the tracks...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TFile.h>
#include <TVector3.h>
#include "RhoNNO/TXMLP.h"

#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

#define TRACKML

#define NHITS 5
#define SIGMA 0.001

#ifdef TRACKML
#define MAXHITS 10000
#define NETFILE "/Users/marcel/workspace/rhonno/trackml/NNO0100.TXMLP"
#define TRACKLET 3
#define DISTANCE 5000
#define THRESHOLD 65
#else
#define MAXHITS 500
#define NETFILE "/Users/marcel/workspace/rhonno/RhoNNO/NNO0069.TXMLP"
#define TRACKLET 3
#define DISTANCE 100
#define THRESHOLD 65
#endif

#define DRAW true
#define VERBOSE false
//#define LOOKUP

// The user member function processes one event

std::vector<TVector3> hits;
std::vector<int> tracks[20000];

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

int findTracks(int nhits, float *x, float *y, float *z, int* labels);
Double_t* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist);

bool sortFunc( const vector<int>& p1,
              const vector<int>& p2 ) {
    return p1.size() > p2.size();
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

void GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error) {
    default_random_engine generator;
    double tau = 0.025;
    for (int i=0; i<np; i++,tau+=delta)
    {
        Float_t X,Y,Z;
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
        points.push_back(TVector3(X,Y,Z));
    }
}


int main(int argc, char* argv[]) {
    
    TFile output("MLPTracker.root","RECREATE");
    
    string filename("event");
    TNtuple nt1("Hits","NNO Tracking Data","x:y:z");
    
    if (argc > 1) {
        filename = argv[1];
        
        ifstream in(filename);
        if (in) {
            cout << "Reading input file: " << filename << endl;
        }
        else {
            cout<<"Event "<<filename<<" does not exist!!"<<endl;
            exit(0);
        }
#ifdef TRACKML
        char tmpLine[256];
        in.getline(tmpLine,256);
        cout<<tmpLine<<endl;
        while (hits.size()<MAXHITS) {
            double h[7]; //file line: id:x:y:z:volume:layer:module
            if( !readLine(in,h,7) ) break;
            if( h[0]-1 != hits.size() ){
                cout<<"Hit index is wrong: "<<h[0]<<endl;
                exit(0);
            }
            TVector3 hit(h[1],h[2],h[3]);
            //hit *= 0.001; // convert mm to m
            hits.push_back(hit);
            if (hits.size()%1000 == 0) cout << hits.size() << endl;
        }
#else
        double X,Y,Z;
        while (in >> X >> Y >> Z) {
            X*=0.01; // transform to meter
            Y*=0.01;
            Z*=0.01;
            TVector3 point(X,Y,Z);
            hits.push_back(point);
            //cout << point.x() << "\t" << point.y() << "\t" << point.z() << "\t"<< point.d() << endl;
        }
#endif
        cout<<" loaded "<<hits.size()<<" hits "<<endl;
        in.close();
        
        /*        double X,Y,Z;
         while (infile >> X >> Y >> Z) {
         X*=0.01; // transform to meter
         Y*=0.01;
         Z*=0.01;
         TVector3 point(X,Y,Z);
         hits.push_back(point);
         cout << point.x() << "\t" << point.y() << "\t" << point.z() << "\t"<< point.Mag() << endl;
         }
         */
        
    }
    else
    {
        // std::vector<TVector3>, int np, float delta tau, float radius, float phi, float gamma
        GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/1.0,0.5,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/3.0,1.5,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/4.0,-1.2,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/2.0,1.0,SIGMA);
    }
    
    // Sort the hits according to distance from origin
    
    //cout << "Sorting hits..." << endl;
    //reverse(hits.begin(),hits.end());
    
    unsigned long nhits = hits.size();
    cout << endl << "Reconstruct Track data with a XMLP Network";
    cout << endl << " Number of hits:" << nhits << endl;
    
    float x[nhits],y[nhits],z[nhits];
    int labels[nhits];
    int nt;
    for (int i=0; i<nhits; i++) {
        TVector3 hit = hits[i];
        x[i] = hit.X();
        y[i] = hit.Y();
        z[i] = hit.Z();
    }
    nt = findTracks((int)nhits,x,y,z,labels);
    cout << "Number of tracks:" << nt << endl;
    for(int i=0; i<nt; i++) {
        cout << "Track " << i << ":";
        print(tracks[i]);
    }
    
    
    cout << "Labels: ";
    for (int i=0;i<nhits;i++) cout << labels[i] << " ";
    
    //TBD: Fit the helix tracks from the hits in the containers
    
    // Initialize a 3D canvas and draw the hits and tracks
    if (DRAW) {
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
        TPolyMarker3D *hitmarker = new TPolyMarker3D((UInt_t) nhits);
        vector<TVector3>::iterator it;
        for(it = hits.begin(); it != hits.end(); it++)    {
            static int i = 0;
            TVector3 p=*it;
            hitmarker->SetPoint(i++,p.x(),p.y(),p.z());
            nt1.Fill(p.x(),p.y(),p.z());
        }
        // set marker size, color & style
        hitmarker->SetMarkerSize(1.0);
        hitmarker->SetMarkerColor(kCyan);
        hitmarker->SetMarkerStyle(kStar);
        hitmarker->Draw();
        
        // Draw the tracks
        TPolyMarker3D *trackmarker = new TPolyMarker3D((UInt_t) hits.size());
        static int n = 0;
        for (int i=0;i<nt;i++) {
            vector<int> h = tracks[i];
            for (int j=0;j<h.size();j++) {
                TVector3 hit = hits[h[j]];
                trackmarker->SetPoint(n++,hit.x(),hit.y(),hit.z());
            }
        }
        // set marker size, color & style
        trackmarker->SetMarkerSize(0.25);
        trackmarker->SetMarkerColor(kRed);
        trackmarker->SetMarkerStyle(kFullDotLarge);
        trackmarker->Draw();
        
        cout << endl << "Writing..." << endl;
        c1->Write();
    }
    
    return EXIT_SUCCESS;
}

// Assign track labels to hits (x,y,z)
// The hit pair quality is assessed by the neural network
// The quality is noted in the hit pair matrix m[nhits][nhits]
int findTracks(int nhits, float *x, float *y, float *z, int* labels)
{
    std::clock_t c_start = std::clock();

    for (int i=0;i<nhits;i++) labels[i] = -1; // Preset with no match

#ifdef LOOKUP
    Double_t in1[7], in2[7];
    
    // Allocate a nhits*nhits hit pair matrix as one continuous memory block
    int **m = new int*[nhits];
    if (nhits)
    {
        m[0] = new int[nhits * nhits];
        for (int i = 1; i < nhits; ++i)
            m[i] = m[0] + i * nhits;
    }
    
    // Allocate a nhits*nhits hit pair distance matrix as one continuous memory block
    int **d = new int*[nhits];
    if (nhits)
    {
        d[0] = new int[nhits * nhits];
        for (int i = 1; i < nhits; ++i)
            d[i] = d[0] + i * nhits;
    }
    
    for(int i=0; i<nhits-1; i++)    {
        in1[0] = x[i];
        in1[1] = y[i];
        in1[2] = z[i];
        in2[3] = x[i];
        in2[4] = y[i];
        in2[5] = z[i];
        d[i][i] = 0;
        m[i][i] = 0;
        // Polar coordinates
        //double phi = atan2(y[i],x[i]);
        //double r = sqrt(x[i]*x[i]+y[i]*y[i]);
        //double zz = z[i];
        for(int j=i+1; j<nhits; j++)    {
            in1[3] = x[j];
            in1[4] = y[j];
            in1[5] = z[j];
            in2[0] = x[j];
            in2[1] = y[j];
            in2[2] = z[j];
            double dist = sqrt((in1[0]-in2[0])*(in1[0]-in2[0]) + (in1[1]-in2[1])*(in1[1]-in2[1]) + (in1[2]-in2[2])*(in1[2]-in2[2]));
            in1[6] = dist;
            in2[6] = dist;
            d[i][j] = 1000.*dist; // Cache the distance between the points
            d[j][i] = 1000.*dist;
            d[j][j] = 0;
            m[j][j] = 0;
            m[i][j] = (int) 100. * Recall(in1[0],in1[1],in1[2],in1[3],in1[4],in1[5],in1[6])[0]; // Recall the hit pair matching quality
            m[j][i] = (int) 100. * Recall(in2[0],in2[1],in2[2],in2[3],in2[4],in2[5],in2[6])[0];
            if (m[i][j]< THRESHOLD) m[i][j] = 0; // Apply a cut on the quality
            if (m[j][i]< THRESHOLD) m[j][i] = 0;
        }
    }
#endif
    
    // Search neighbouring hits, the neural network recall identifies the hit belonging to a tracklet
    vector<vector<int>> tracklet;
    for(int i=0; i<nhits-1; i++)    {
        vector<int> tmpvec;
        tmpvec.push_back(i); //Note the row in the first place
        for(int j=i+1; j<nhits; j++)    {
#ifdef LOOKUP
            int dist = d[i][j];
            if (dist > DISTANCE) continue;
            int recall  = (m[i][j]>m[j][i]) ? m[i][j]:m[j][i];
#else
            double d = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));
            int dist = 1000.*d;
            if (dist > DISTANCE) continue;
            int recall1 = (int) 100. * Recall(x[i],y[i],z[i],x[j],y[j],z[j],d)[0]; // Recall the hit pair matching quality
            int recall2 = (int) 100. * Recall(x[j],y[j],z[j],x[i],y[i],z[i],d)[0]; // Recall the hit pair matching quality
            if (recall1 < THRESHOLD) recall1 = 0; // Apply a cut on the quality
            if (recall2 < THRESHOLD) recall2 = 0;
            int recall  = (recall1>recall2) ? recall1:recall2;
#endif
            if (recall>THRESHOLD) {
                tmpvec.push_back(j); // Note the columns with a good combination
            }
        }
        
        tracklet.push_back(tmpvec);
    }
    
    cout << "Number of tracklets: " << tracklet.size() << endl;
    
    // Sort the tracklet vector according to the tracklet length
    
    sort(tracklet.begin(), tracklet.end(), sortFunc);
    
    // Print out the sorted vector
    if (VERBOSE) {
        cout << "Sorted tracklets:" << endl;
        for( int i=0; i<tracklet.size(); i++ ) {
            for( int j=0; j<tracklet[i].size(); j++ ) {
                int row = tracklet[i][0];
                int col = tracklet[i][j];
#ifdef LOOKUP
                cout << tracklet[i][j] << "(" << m[row][col] << ") ";
#else
                double dist = sqrt((x[row]-x[col])*(x[row]-x[col]) + (y[row]-y[col])*(y[row]-y[col]) + (z[row]-z[col])*(z[row]-z[col]));
                int recall1 = (int) 100. * Recall(x[row],y[row],z[row],x[col],y[col],z[col],dist)[0]; // Recall the hit pair matching quality
                int recall2 = (int) 100. * Recall(x[col],y[col],z[col],x[row],y[row],z[row],dist)[0]; // Recall the hit pair matching quality
                if (recall1 < THRESHOLD) recall1 = 0; // Apply a cut on the quality
                if (recall2 < THRESHOLD) recall2 = 0;
                int recall  = (recall1>recall2) ? recall1:recall2;
                if (j==0) recall = 0.;
                cout << tracklet[i][j] << "(" << recall << ") ";
#endif
            }
            cout << endl;
        }
    }

    cout << "Seed: " << tracklet[0][0] << " length: " << tracklet[0].size() << endl;
    
    // Prune the tracklets by removing short tracks
    for (vector<vector<int>>::iterator it = tracklet.begin(); it != tracklet.end(); ++it) {
        vector<int> row = *it;
        if (row.size() < TRACKLET) { // Remove short tracklets
            tracklet.erase(it);
            *it--;
            continue;
        }
    }
    
    // Print out the pruned vector
    if (VERBOSE) {
        cout << "Pruned tracklets (remove short paths):" << endl;
        for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    }

    // Prune the tracklets by removing rows with identical entries
    // Assemble tracks from the corresponding tracklets
    vector<vector<int>> track;
    for (vector<vector<int>>::iterator it = tracklet.begin(); it != tracklet.end(); ++it) {
        vector<int> row = *it;
        vector<int> tmpvec = row; // Vector to assemble the track
        for (vector<int>::iterator it2 = row.begin(); it2 != row.end(); ++it2) { // Search next seeding hits in remainder tracklet list
            int prune = *it2;
            for (vector<vector<int>>::iterator it3 = it+1; it3 != tracklet.end(); ++it3) {
                vector<int> nextrow = *it3;
                bool hitExists = find(nextrow.begin(),nextrow.end(),prune) != nextrow.end();
                if (hitExists) {
                    for (int j=0;j<nextrow.size();j++) tmpvec.push_back(nextrow[j]); // Append the hits to track before erasing the row
                    tracklet.erase(it3);
                    *it3--;
                    continue;
                }
            }
        }
        set<int> s( tmpvec.begin(), tmpvec.end() ); // Remove duplicates
        tmpvec.assign( s.begin(), s.end() );
        track.push_back(tmpvec);
    }
    
    // Print out the pruned vector
    if (VERBOSE) {
        cout << "Pruned tracklets (Removed duplicates):" << endl;
        for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    }

    // Print out the tracks vector
    cout << "Tracks:" << endl;
    for( int i=0; i<track.size(); i++ ) print(track[i]);
    
#ifdef LOOKUP
    if (nhits) delete [] m[0]; // Clean up the memory
    delete [] m;
    
    if (nhits) delete [] d[0]; // Clean up the memory
    delete [] d;
#endif

    for (int i=0;i<track.size();i++) {
        tracks[i] = track[i]; // Save the results
        for (int j=0;j<track[i].size();j++) {
            int hit = track[i][j];
            labels[hit] = i;
        }
    }
    
    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

    return (int) track.size();
}

#ifdef TRACKML
// Recall function on normalised network input
double* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist)
{
    static TXMLP net(NETFILE);
    float x[7];
    x[0]     = 0.775527    *    x1;    // x1
    x[1]     = 0.71624     *    y1;    // y1
    x[2]     = 0.100536    *    z1;    // z1
    x[3]     = 0.775527    *    x2;    // x2
    x[4]     = 0.71624     *    y2;    // y2
    x[5]     = 0.100536    *    z2;    // z2
    x[6]     = 0.0000245   *    dist;
    return net.Recallstep(x);
}
#else
// Recall function on normalised network input
double* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist)
{
    static TXMLP net(NETFILE);
    Float_t x[7];
    x[0]     = 2631.15    *    x1;    // x1
    x[1]     = 383.788    *    y1;    // y1
    x[2]     = 490.839    *    z1;    // z1
    x[3]     = 2631.15    *    x2;    // x2
    x[4]     = 383.788    *    y2;    // y2
    x[5]     = 490.839    *    z2;    // z2
    x[6]     = 126.061    *    dist;
    return net.Recallstep(x);
}
#endif
