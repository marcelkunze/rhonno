// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"

// Assign track labels to hits (x,y,z)
// Hits are sorted according to their distance from the origin
// The hit pair quality is assessed by neural network function recall2
// Starting from the innermost hit a seeding hit pair is identified and extended to a tracklet
// The tracklet quality is assessed by neural network function recall3
// The tracklets are extended until no further suiting hit is found
// Corresponding labels are assigned to hits

#include "TXMLP.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

Point::Point(double x, double y, double z, int id, int val)
{
    _id = id;
    _val = val;
    _x = x;
    _y = y;
    _z = z;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
}

Point::Point(float x, float y, float z, int id=-1, int val=-1)
{
    _id = id;
    _val = val;
    _x = x;
    _y = y;
    _z = z;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
}

Point::Point(const Point &p)
{
    _id = p._id;
    _val = p._val;
    _x = p._x;
    _y = p._y;
    _z = p._z;
    _r = p._r;
    _phi = p._phi;
    _theta = p._theta;
    _distance = p._distance;
}

// Print an integer hit id track vector
void Tracker::print(vector<int> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
    cout << endl;
}

// Sort an hit id track vector by length
bool Tracker::sortFunc( const vector<int>& p1,
              const vector<int>& p2 ) {
    return p1.size() > p2.size();
}

// Recall function for 2 points
double* Tracker::Recall2(Point &p1, Point &p2)
{
    static TXMLP net(NETFILE2);
    static float x[6]={0.,0.,0.,0.,0.,0.};
    static double bad[1]={-1.0};
    static double null[1]={0.0};

    double angle = acos(Point::dot(p1,p2)); // Check angle between the vectors (origin,p1) and (origin,p2)
    if (angle > DELTAPHI) return null;

    double deltar = abs(p1.r()-p2.r()); // Check the radial distance
    if (deltar > DELTAR) return bad;
    
    x[0]     = p1.r();     // r1
    x[1]     = p1.phi();   // phi1
    x[2]     = p1.theta(); // theta1
    x[3]     = p2.r();     // r2
    x[4]     = p2.phi();   // phi2
    x[5]     = p2.theta(); // theta2
    
    return net.Recallstep(x);
}

// Recall function for 3 points
double* Tracker::Recall3(Point &p1, Point &p2, Point &p3)
{
    static TXMLP net(NETFILE3);
    static float x[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    static double bad[1]={-1.0};
    static double null[1]={0.0};

    Point v1 = p2 - p1;
    Point v2 = p3 - p1;
    double angle = acos(Point::dot(v1,v2)); // Check angle between the last points
    if (angle > DELTAPHI) return null;
    
    x[0]     = p1.r();     // r1
    x[1]     = p1.phi();   // phi1
    x[2]     = p1.theta(); // theta1
    x[3]     = p2.r();     // r2
    x[4]     = p2.phi();   // phi2
    x[5]     = p2.theta(); // theta2
    x[6]     = p3.r();     // r3
    x[7]     = p3.phi();   // phi3
    x[8]     = p3.theta(); // theta3
    
    return net.Recallstep(x);
}

int Tracker::findTracks(int nhits, float *x, float *y, float *z, int* labels)
{
    std::clock_t c_start = std::clock();
    
    Point *p = new Point[nhits];
    vector<Point> points;
    points.reserve(nhits);
    
    unsigned long n2(0), n3(0);

    // Set up a cache for the point coordinates
    //cout << "Set up points cache..." << endl;
    for (int i=0;i<nhits;i++) {
        labels[i] = 0;
        p[i] = Point(x[i],y[i],z[i],i);
        points.push_back(p[i]);
    }
    
    // Sort the hits according to distance from origin
    cout << "Sorting hits..." << endl;
    sort(points.begin(),points.end(),Point::sortDist);
    
    // Search neighbouring hits, the neural network recall identifies the hit belonging to a tracklet
    cout << "Find tracklets..." << endl;
    int nd(0), nr(0), nn(0), np(0);
    vector<vector<int> > tracklet;
    for (vector<Point>::iterator it1 = points.begin(); it1 != points.end(); ++it1) {
        Point p0 = *it1; // Seeding point
        Point p1 = *it1;
        vector<Point> pvec;
        pvec.push_back(p1); //Note the seeding point in the first place
        if (VERBOSE) cout << endl << p1.id() << "(0) ";
        int ntry = 0;
        for (vector<Point>::iterator it2 = it1+1; it2 != points.end(); ++it2) { //
            Point p2 = *it2;
            double dist = p2.distance(p1); // Only consider points in the neighborhood
            if (dist > DISTANCE) {
                nd++;
                continue;
            }
            
            // Seeding a second point: Get network track quality of 2 points
            double recall = 0.0;
            if (pvec.size() <2) {
                recall = Recall2(p1,p2)[0]; // try the initial pair
                n2++;
                if (recall == 0 && ntry++<5) { np++; continue; } // try the next 5 points
                if (recall < 0.0) { nr++; break; } // Out of bounds; finish tracklet
            }
            else
            {
                recall = Recall3(p0,p1,p2)[0]; // Get network track quality of 3 consecutive points
                n3++;
                if (recall <= 0.0) { np++; continue; } // No straight conection between the three points
            }

            if (recall>THRESHOLD) {
                pvec.push_back(p2); // Note the columns with a good combination
                if (VERBOSE) cout << p2.id() << "(" << (int) 100*recall << ") ";
                points.erase(it2);  // Remove the corresponding point from the set
                *it2--;
                p0 = p1;// Note the assigned hits
                p1 = p2;
                continue;
            }
            else
                nn++;
            
        }
        points.erase(it1);  // Remove the corresponding point from the set
        *it1--;
        
        if (pvec.size() < TRACKLET) continue; // Perform a cut on tracklet size
        sort(pvec.begin(), pvec.end(), Point::sortId); // Sort the hits acording to the Id
        vector<int> tmpvec;
        for (int ip=0;ip<pvec.size();ip++) tmpvec.push_back(pvec[ip].id()); // Note the hit indices
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
            p[hit].setval(i+1);
            n++;
        }
    }
    
    cout << "Number of assigned points: " << n << endl;
    cout << "Threshold <" << THRESHOLD << ": " << nn <<endl;
    cout << "Distance <" << DISTANCE << ": " << nd <<endl;
    cout << "Radius   <" << DELTAR << ": " << nr <<endl;
    cout << "Phi      <" << DELTAPHI << ": " << np <<endl;
    cout << "Recalls2 : " << n2 << endl;
    cout << "Recalls3 : " << n3 << endl;

    delete [] p;
    
    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    return (int) tracklet.size();
}


