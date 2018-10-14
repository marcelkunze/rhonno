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

#define NHITS 10
#define SIGMA 0.0

#ifdef TRACKML
#define MAXHITS 50000
#define NETFILE2 "/Users/marcel/workspace/rhonno/trackml/NNO0200-6-25-15-1.TXMLP"
#define TRACKLET 3
#define THRESHOLD 98
#define DISTANCE 1.0
#define DELTAR   0.1
#define DELTAPHI 0.01
#define DELTATHETA 0.05

#else
#define MAXHITS 150000
#define NETFILE2 "/Users/marcel/workspace/rhonno/RhoNNO/NNO0100.TXMLP"
#define NETFILE3 "/Users/marcel/workspace/rhonno/RhoNNO/NNO0099.TXMLP"
#define TRACKLET 3
#define THRESHOLD 90
#define DISTANCE 1.0
#define DELTAR   0.1
#define DELTAPHI 0.043
#define DELTATHETA 0.08
#endif

#define DRAW true
#define VERBOSE true

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

int findTracks(int nhits, float *x, float *y, float *z, int* labels);
double* Recall2(float,float,float,float,float,float);
double* Recall3(float,float,float,float,float,float,float,float,float);

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

void GenerateTrack(std::vector<Point> &points, int np, double delta, double radius, double phi, double gamma, double error) {
    default_random_engine generator;
    double tau = 0.1;
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
        Point point;
        point.x = X;
        point.y = Y;
        point.z = Z;
        points.push_back(point);
    }
}

std::vector<Point> hits;
std::vector<int> tracks[20000];

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
            Point point;
            point.x = h[1] * 0.001; // convert mm to m
            point.y = h[2] * 0.001; // convert mm to m
            point.z = h[3] * 0.001; // convert mm to m
            hits.push_back(point);
            if (hits.size()%1000 == 0) cout << hits.size() << endl;
        }
#else
        int n = 0;
        double X,Y,Z;
        while (in >> X >> Y >> Z) {
            if (n++ >= MAXHITS) break;
            Point point;
            point.x = X * 0.01; // convert cm to m
            point.y = Y * 0.01; // convert cm to m
            point.z = Z * 0.01; // convert cm to m
            hits.push_back(point);
            //cout << point.x() << "\t" << point.y() << "\t" << point.z() << "\t"<< point.d() << endl;
        }
#endif
        cout<<" loaded "<<hits.size()<<" hits "<<endl;
        in.close();
    }
    else
    {
        // std::vector<Point>, int np, float delta tau, float radius, float phi, float gamma
        GenerateTrack(hits,NHITS,0.025, 0.5,M_PI/2.0, 1.0,SIGMA); // 00
        GenerateTrack(hits,NHITS,0.025,-1.0,M_PI/4.0, 1.0,SIGMA); // 10
        GenerateTrack(hits,NHITS,0.025, 1.5,M_PI/3.0,-1.0,SIGMA); // 20
        GenerateTrack(hits,NHITS,0.025,-2.0,M_PI/3.0,-1.0,SIGMA); // 30
    }
    
    unsigned long nhits = hits.size();
    cout << endl << "Reconstruct Track data with a XMLP Network";
    cout << endl << " Number of hits:" << nhits << endl;
    
    float x[nhits],y[nhits],z[nhits];
    int labels[nhits];
    int nt;
    for (int i=0; i<nhits; i++) {
        Point hit = hits[i];
        x[i] = hit.x;
        y[i] = hit.y;
        z[i] = hit.z;
    }
    
    nt = findTracks((int)nhits,x,y,z,labels);
    
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
        for (vector<Point>::iterator it = hits.begin(); it != hits.end(); it++)    {
            static int i = 0;
            Point p=*it;
            hitmarker->SetPoint(i++,p.x,p.y,p.z);
            nt1.Fill(p.x,p.y,p.z);
        }
        // set marker size, color & style
        hitmarker->SetMarkerSize(1.0);
        hitmarker->SetMarkerColor(kCyan);
        hitmarker->SetMarkerStyle(kStar);
        hitmarker->Draw();
        
        // Draw the tracks
        for (int i=0;i<nt;i++) {
            //cout << endl << "Drawing track " << i+1 << ": ";
            vector<int> h = tracks[i];
            int n = 0;
            TPolyLine3D *connector = new TPolyLine3D((int)h.size());
            for (vector<int>::iterator it = h.begin(); it != h.end(); it++)    {
                int index = *it;
                //cout << index << " ";
                Point hit = hits[index];
                connector->SetPoint(n++, hit.x, hit.y, hit.z);
            }
            connector->SetLineWidth(1);
            connector->SetLineColor(kRed);
            connector->Draw();
        }
        
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

    // Sort the hits according to distance from origin
    cout << "Sorting hits..." << endl;
    sort(points.begin(),points.end(),sortDist);

    // Search neighbouring hits, the neural network recall identifies the hit belonging to a tracklet
    cout << "Find tracklets..." << endl;
    int nd(0), nr(0), np(0), nt(0), nn(0);
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
            int recall1 = (int) 100. * Recall2(p1.r,p1.phi,p1.theta,p2.r,p2.phi,p2.theta)[0]; // Recall the hit pair matching quality
            int recall2 = (int) 100. * Recall2(p2.r,p2.phi,p2.theta,p1.r,p1.phi,p1.theta)[0]; // Recall the hit pair matching quality
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
            else
                nn++;

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
    cout << "Threshold <" << THRESHOLD << ": " << nn <<endl;
    cout << "Distance <" << DISTANCE << ": " << nd <<endl;
    cout << "Radius   <" << DELTAR << ": " << nr <<endl;
    cout << "Phi      <" << DELTAPHI << ": " << np <<endl;
    cout << "Theta    <" << DELTATHETA << ": " << nt <<endl;

    delete [] p;
    
    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    return (int) tracklet.size();
}

#include "TXMLP.h"

// Recall function for 2 points
double* Recall2(float x1,float x2,float x3,float x4,float x5,float x6)
{
    static TXMLP net(NETFILE2);
    float x[6];
    x[0]     = x1;    // r1
    x[1]     = x2;    // phi1
    x[2]     = x3;    // theta1
    x[3]     = x4;    // r2
    x[4]     = x5;    // phi2
    x[5]     = x6;    // theta2
    return net.Recallstep(x);
}

// Recall function for 3 points
double* Recall3(float x1,float x2,float x3,float x4,float x5,float x6,float x7, float x8, float x9)
{
    static TXMLP net(NETFILE3);
    float x[9];
    x[0]     = x1;    // r1
    x[1]     = x2;    // phi1
    x[2]     = x3;    // theta1
    x[3]     = x4;    // r2
    x[4]     = x5;    // phi2
    x[5]     = x6;    // theta2
    x[6]     = x7;    // r3
    x[7]     = x8;    // phi3
    x[8]     = x9;    // theta3
    return net.Recallstep(x);
}

