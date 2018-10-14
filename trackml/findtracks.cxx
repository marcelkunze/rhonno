#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

#include "XMLP.h"

using namespace std;

#define NETFILE "~/workspace/rhonno/trackml/NNO0200-6-25-15-1.TXMLP"
#define TRACKLET 3
#define THRESHOLD 98
#define DISTANCE 1.0
#define DELTAR   0.1
#define DELTAPHI 0.043
#define DELTATHETA 0.08

#define VERBOSE false

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

double distanceBetween(const Point &a,const Point &b)
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
    for (int i = 0; i < n; i++) arr[i].distance = distanceBetween(arr[i],p);
    
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

// Recall function on network input
double* Recall(float x1, float x2, float x3, float x4, float x5, float x6, const char *netfile=NETFILE)
{
    //cout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << endl;
    static XMLP net(netfile);
    float x[6];
    x[0]     = x1;    // r1
    x[1]     = x2;    // phi1
    x[2]     = x3;    // theta1
    x[3]     = x4;    // r2
    x[4]     = x5;    // phi2
    x[5]     = x6;    // theta2
    return net.Recallstep(x);
}

// Assign track labels to hits (x,y,z)
// The hit pair quality is assessed by the neural network
// The quality is noted in the hit pair matrix m[nhits][nhits]
int findTracks(int nhits,float *x,float *y,float *z,int* labels,float threshold=THRESHOLD,float distance=DISTANCE,float deltar=DELTAR,float deltaphi=DELTAPHI,float deltatheta=DELTATHETA,const char *netfile=NETFILE)
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
    cout << "Finding tracklets..." << endl;
    int nd(0), nr(0), np(0), nt(0), nn(0);
    Point vertex;
    vertex.x = 0;
    vertex.y = 0;
    vertex.z = 0;
    vector<vector<int> > tracklet;
    for (vector<Point>::iterator it1 = points.begin(); it1 != points.end(); ++it1) {
        Point p1 = *it1; // Seeding point
        vector<Point> pvec;
        // Conformal mapping of circle to straight line
        pvec.push_back(p1); //Note the seeding point in the first place
        double r1 = 0.0;
        if (VERBOSE) cout << endl << p1.id << "(0) ";
        for (vector<Point>::iterator it2 = it1+1; it2 != points.end(); ++it2) { //
            Point p2 = *it2;
            double dist = distanceBetween(p1,p2); // Only consider points in the neighborhood
            double r2 = 0.0;
            if (pvec.size()==2) r1 = circleRadius(pvec[0],pvec[1],p2);
            if (pvec.size()>=2)  r2 = circleRadius(pvec[0],pvec[1],p2);
            double dr = abs(r1-r2);
            double dphi = abs(abs(p1.phi)-abs(p2.phi));
            double dtheta = abs(p1.theta-p2.theta);
            //cout << endl << p1.id << " " << p2.id << " " << dist << " " << deltar << " " << deltaphi << " " << deltatheta << endl;
            if (dist > distance) {
                //if (VERBOSE) cout << p1.id << " " << p2.id << " dist:" << dist << endl;
                nd++;
                continue;
            }
            if (dr > deltar) {
                //if (VERBOSE && pvec.size()>=2) cout << pvec[0].id << " " << pvec[1].id << " " << p2.id << " deltar:" << dr << endl;
                nr++;
                continue;
            }
            if (dphi > deltaphi) {
                //if (VERBOSE) cout << p1.id << " " << p2.id << " deltaphi:" << dphi << endl;
                np++;
                continue;
            }
            if (dtheta > deltatheta) {
                //if (VERBOSE) cout << p1.id << " " << p2.id << " deltatheta:" << dtheta << endl;
                nt++;
                continue;
            }
            int recall1 = (int) 100. * Recall(p1.r,p1.phi,p1.theta,p2.r,p2.phi,p2.theta,netfile)[0]; // Recall the hit pair matching quality
            int recall2 = (int) 100. * Recall(p2.r,p2.phi,p2.theta,p1.r,p1.phi,p1.theta,netfile)[0]; // Recall the hit pair matching quality
            if (recall1 < threshold) recall1 = 0; // Apply a cut on the quality
            if (recall2 < threshold) recall2 = 0;
            int recall  = (recall1>recall2) ? recall1:recall2;
            if (recall>threshold) {
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

    // Print out the tracks vector
    if (VERBOSE) {
        cout << "Tracks:" << endl;
        for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    }
    
    int n = 0;
    for (int i=0;i<tracklet.size();i++) {
        for (int j=0;j<tracklet[i].size();j++) {
            int hit = tracklet[i][j];
            labels[hit] = i+1;
            p[hit].val = i+1;
            n++;
        }
    }
    
    cout << "Number of assigned points: " << n << endl;
    cout << "Threshold <" << threshold << ": " << nn <<endl;
    cout << "Distance  <" << distance << ": " << nd <<endl;
    cout << "Radius    <" << deltar << ": " << nr <<endl;
    cout << "Phi       <" << deltaphi << ": " << np <<endl;
    cout << "Theta     <" << deltatheta << ": " << nt <<endl;
    
    delete [] p;
    
    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
    
    return (int) tracklet.size();
}

extern "C" {
    int processHits(int nhits,float *x,float *y,float *z,int* labels,float threshold=THRESHOLD,float distance=DISTANCE,float deltar=DELTAR,float deltaphi=DELTAPHI,float deltatheta=DELTATHETA)
    {
        const char *netfile(NETFILE);
        return findTracks(nhits,x,y,z,labels,threshold,distance,deltar,deltaphi,deltatheta,netfile);
    }
    
    float processRecall(float x1, float x2, float x3, float x4, float x5, float x6, const char *netfile=NETFILE)
    {
        return (float) Recall(x1,x2,x3,x4,x5,x6,netfile)[0];
    }

}

