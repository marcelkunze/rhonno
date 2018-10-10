// Read HITS spatial data and try a XMLP network to segment the tracks...

#include "TROOT.h"
#include "TCanvas.h"
#include "TView.h"
#include "TPolyMarker3D.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TNtuple.h"

#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>

using namespace std;

//#define TRACKML

#define NHITS 5
#define SIGMA 0.001

#ifdef TRACKML
#define MAXHITS 150000
#define NETFILE "/Users/marcel/workspace/rhonno/trackml/NNO0200-6-25-15-1.TXMLP"
#define TRACKLET 3
#define THRESHOLD 90
#define DISTANCE 100
#define DELTAR   0.1
#define DELTAPHI 0.1
#define DELTATHETA 0.1

#else
#define MAXHITS 500
#define NETFILE "/Users/marcel/workspace/rhonno/RhoNNO/NNO0100.TXMLP"
#define TRACKLET 2
#define THRESHOLD 50
#define DISTANCE 100
#define DELTAR   0.1
#define DELTAPHI 0.1
#define DELTATHETA 0.1
#endif

#define DRAW true
#define VERBOSE false

// The user member function processes one event

std::vector<TVector3> hits;
std::vector<int> tracks[20000];

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

int findTracks(int nhits, float *x, float *y, float *z, int* labels);
double* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist);

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

// Used to sort an array of points by increasing
// order of distance
bool comparison(const Point &a,const Point &b)
{
    return (a.distance < b.distance);
}

// This function finds classification of point p using
// k nearest neighbour algorithm. It assumes only two
// groups and returns 0 if p belongs to group 0, else
// 1 (belongs to group 1).
int classifyAPoint(Point arr[], int n, int k, Point p)
{
    // Fill distances of all points from p
    for (int i = 0; i < n; i++)
        arr[i].distance =
        sqrt((arr[i].x - p.x) * (arr[i].x - p.x) +
             (arr[i].y - p.y) * (arr[i].y - p.y) +
             (arr[i].z - p.z) * (arr[i].z - p.z));
    
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
            hit *= 0.001; // convert mm to m
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
    
#define MAXTRACK 25
    cout << "Number of tracks:" << nt << endl;
    for(int i=0; i<nt; i++) {
        if (i<MAXTRACK || i>nt-MAXTRACK) {
            cout << "Track " << i+1 << ":";
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
        TPolyMarker3D *hitmarker = new TPolyMarker3D((unsigned int) nhits);
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
        TPolyMarker3D *trackmarker = new TPolyMarker3D((unsigned int) hits.size());
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

    // Search neighbouring hits, the neural network recall identifies the hit belonging to a tracklet
    cout << "Find tracklets..." << endl;
    vector<vector<int>> tracklet;
    for (vector<Point>::iterator it1 = points.begin(); it1 != points.end(); ++it1) {
        Point point1 = *it1;
        int i = point1.id;
        vector<int> tmpvec;
        tmpvec.push_back(i); //Note the index of point in the first place
        for (vector<Point>::iterator it2 = it1+1; it2 != points.end(); ++it2) { //
            Point point2 = *it2;
            int j = point2.id; // Index of second point
            double d = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));
            int dist = 1000.*d;
            if (dist > DISTANCE) continue;
            float deltar = abs(p[i].r-p[j].r);
            if (deltar > DELTAR) continue;
            float deltaphi = abs(abs(p[i].phi)-abs(p[j].phi));
            if (deltaphi > DELTAPHI) continue;
            float deltatheta = abs(p[i].theta-p[j].theta);
            if (deltatheta > DELTATHETA) continue;
            int recall1 = (int) 100. * Recall(p[i].r,p[i].phi,p[i].theta,p[j].r,p[j].phi,p[j].theta,d)[0]; // Recall the hit pair matching quality
            int recall2 = (int) 100. * Recall(p[j].r,p[j].phi,p[j].theta,p[i].r,p[i].phi,p[i].theta,d)[0]; // Recall the hit pair matching quality
            if (recall1 < THRESHOLD) recall1 = 0; // Apply a cut on the quality
            if (recall2 < THRESHOLD) recall2 = 0;
            int recall  = (recall1>recall2) ? recall1:recall2;
            if (recall>THRESHOLD) {
                tmpvec.push_back(j); // Note the columns with a good combination
                points.erase(it2);  // Remove the corresponding point from the set
                *it2--;
                continue;
            }
        }
        points.erase(it1);  // Remove the corresponding point from the set
        *it1--;

        if (tmpvec.size() < TRACKLET) continue; // Perform a cut on tracklet size
        tracklet.push_back(tmpvec);
    }
    
    cout << "Number of unassigned points:" << points.size() << endl;
    cout << "Number of tracklets: " << tracklet.size() << endl;
    if (tracklet.size() == 0) exit(0);
    
    // Sort the tracklet vector according to the tracklet length
    
    sort(tracklet.begin(), tracklet.end(), sortFunc);
    
    // Print out the sorted vector
    if (VERBOSE) {
        cout << "Sorted tracklets:" << endl;
        for (vector<Point>::iterator it1 = points.begin(); it1 != points.end(); ++it1) {
            Point point1 = *it1;
            int i = point1.id;
            vector<int> tmpvec;
            tmpvec.push_back(i); //Note the index of point in the first place
            for (vector<Point>::iterator it2 = it1+1; it2 != points.end(); ++it2) { //
                Point point2 = *it2;
                int j = point2.id; // Index of second point
                int recall1 = (int) 100. * Recall(p[i].r,p[i].phi,p[i].theta,p[j].r,p[j].phi,p[j].theta,0)[0]; // Recall the hit pair matching quality
                int recall2 = (int) 100. * Recall(p[j].r,p[j].phi,p[j].theta,p[i].r,p[i].phi,p[i].theta,0)[0]; // Recall the hit pair matching quality
                int recall  = (recall1>recall2) ? recall1:recall2;
                cout << tracklet[i][j] << "(" << recall << ") ";
            }
            cout << endl;
        }
    }
/*
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
*/
    // Print out the tracks vector
    if (VERBOSE) {
        cout << "Tracks:" << endl;
        for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    }
    
    for (int i=0;i<tracklet.size();i++) {
        tracks[i] = tracklet[i]; // Save the results
        for (int j=0;j<tracklet[i].size();j++) {
            int hit = tracklet[i][j];
            labels[hit] = i+1;
            p[hit].val = i+1;
        }
    }
    
    delete [] p;
    
    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

    return (int) tracklet.size();
}

#include "TXMLP.h"

#ifdef TRACKML
// Recall function on normalised network input
double* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist)
{
    static TXMLP net(NETFILE);
    float x[7];
    x[0]     = x1;    // x1
    x[1]     = y1;    // y1
    x[2]     = z1;    // z1
    x[3]     = x2;    // x2
    x[4]     = y2;    // y2
    x[5]     = z2;    // z2
    x[6]     = dist;
    return net.Recallstep(x);
}
#else
// Recall function on normalised network input
double* Recall(float x1, float y1, float z1, float x2, float y2, float z2, float dist)
{
    static TXMLP net(NETFILE);
    float x[7];
    x[0]     = 1.98932    *    x1;    // r1
    x[1]     = 1.48747    *    y1;    // phi1
    x[2]     = 0.0427916  *    z1;    // theta1
    x[3]     = 1.98932    *    x2;    // r2
    x[4]     = 1.48747    *    y2;    // phi2
    x[5]     = 0.0427916  *    z2;    // theta2
    x[6]     = 126.061    *    dist;
    return net.Recallstep(x);
}
#endif
