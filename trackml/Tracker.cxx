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

#include "XMLP.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>

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
    for (int i=0; i<NEIGHBOURS;i++) _neighbour[i] = -1;
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
    for (int i=0; i<NEIGHBOURS;i++) _neighbour[i] = -1;
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
    for (int i=0; i<NEIGHBOURS;i++) _neighbour[i] = p._neighbour[i];
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
    static XMLP net(NETFILE2);
    static float x[6]={0.,0.,0.,0.,0.,0.};

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
    static XMLP net(NETFILE3);
    static float x[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
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

// Recall function for 4 points
double* Tracker::Recall4(Point &p1, Point &p2, Point &p3, Point &p4)
{
    static XMLP net(NETFILE4);
    static float x[12]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    //static double bad[1]={-1.0};
    static double null[1]={0.0};
    
    Point v1 = p2 - p1;
    Point v2 = p3 - p1;
    double angle = acos(Point::dot(v1,v2)); // Check angle between the last points
    if (angle > DELTAPHI) return null;
    
    Point v3 = p3 - p2;
    Point v4 = p4 - p2;
    angle = acos(Point::dot(v3,v4)); // Check angle between the last points
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
    x[9]     = p4.r();     // r4
    x[10]    = p4.phi();   // phi4
    x[11]    = p4.theta(); // theta4
    
    return net.Recallstep(x);
}

int Tracker::findTracks(int nhits, float *x, float *y, float *z, int* labels)
{
    std::clock_t c_start = std::clock();
    
    Point *p = new Point[nhits];
    vector<Point> points;
    points.reserve(nhits);
    map<int,Point> hitmap;

    unsigned long n1(0), n2(0), n3(0);

    // Set up a cache for the point coordinates
    //cout << "Set up points cache..." << endl;
    for (int i=0;i<nhits;i++) {
        labels[i] = 0;
        p[i] = Point(x[i],y[i],z[i],i);
        points.push_back(p[i]);
    }
    
    // Sort the hits according to distance from origin
    for (int i=0;i<nhits;i++) points[i] = p[i];
    cout << "Sorting " << nhits << " hits..." << endl;
    sort(points.begin(),points.end(),Point::sortRad);
    
    // Search neighbouring hits
    cout << "Searching neighbours..." << endl;
    
    for (int i=0;i<nhits;i++) {
        int index = i;
        if (i>nhits-NEIGHBOURS) index = nhits-NEIGHBOURS;
        int n = 0;
        vector<Point> neighbours;
        Point p0 = points[i];
        
        while (n<NEIGHBOURS && index++<nhits-1) {
            
            if (index == i) continue;
            Point &p1 = points[index];
            
            double angle = acos(Point::dot(p0,p1)); // Check angle between the vectors (origin,p1) and (origin,p2)
            if (angle > DELTAPHI) continue;
            
            double deltar = abs(p1.r()-p0.r()); // Check the radial distance
            if (deltar > DELTAR) continue;
            
            double d = p1.distance(points[i]);
            if (d < DISTANCE) {
                if (Recall2(p0,p1)[0]>THRESHOLD) {
                    neighbours.push_back(p1);
                    n++;
                }
            }
        }
        
        sort(neighbours.begin(),neighbours.end(),Point::sortDist);
        
        for (int j=0;j<neighbours.size();j++) {
            points[i].setneighbour(neighbours[j].id(),j);
            p[points[i].id()].setneighbour(neighbours[j].id(),j);
        }
        
    }
        
    if (VERBOSE) {
        for (int i=0;i<nhits;i++) {
            cout << points[i].id() << "(" ;
            for (int j=0;j<NEIGHBOURS;j++) cout << points[i].neighbour(j) << " ";
            cout << ") ";
        }
        cout << endl;
    }
    
    // fill the hitmap
    for (int i=0;i<nhits;i++) hitmap[points[i].id()] = points[i];

    // Swimmer
    cout << "Find tracklets..." << endl;
    vector<vector<int> > tracklet;
    vector<vector<int> > shortpath;
    int nd(0), nr(0), nn(0), np(0);
    
    while (hitmap.size()>0) {
        if (VERBOSE && n1%1000==0) cout << endl << n1 << " size " << hitmap.size() << endl;
        n1++;
        vector<Point> pvec;
        auto it = hitmap.begin();
        if (VERBOSE) cout << endl << it->first << "(0) ";
        Point p0 = it->second;
        pvec.push_back(p0);
        hitmap.erase(it++);
        while (it != hitmap.end()) {
            int neighbour = p0.neighbour();
            if (neighbour == -1) break;
            Point p1 = p[neighbour];
            
            double recall = 0.0;
            
            if (pvec.size()<2 || p1.neighbour() == -1) {
                recall = Recall2(p0,p1)[0]; // try the initial or last pair
                if (VERBOSE) cout << "/" << p0.id() << " " << p1.id() << "/ ";
                n2++;
            }
            else
            {
                Point p2 = p[p1.neighbour()];
                recall = Recall3(p0,p1,p2)[0]; // Get network track quality of 3 consecutive points
                n3++;
                if (recall == 0.0) { // No straight conection between the three points
                    p2 = p[p1.neighbour(1)]; // Try next neighbour
                    recall = Recall3(p0,p1,p2)[0]; // Get network track quality of 3 consecutive points
                    n3++;
                    if (recall == 0.0) { // No straight conection between the three points
                        if (VERBOSE) cout << endl << "break " << p0.id() << " " << p1.id() << " " << p2.id() << "(" << recall << ")";
                        np++;
                        break;
                    }
                }
            }
  
            if (recall>THRESHOLD) {
                if (VERBOSE) cout << "*" << p1.id() << "(" << (int) 100*recall << ") ";
                pvec.push_back(p1);
                it = hitmap.find(neighbour);
                if (it == hitmap.end()) break;
                hitmap.erase(it++);
                p0 = p1;
            }
            else {
                if (VERBOSE) cout << "*" << p1.id() << "{" << (int) 100*recall << "} ";
                nn++;
                break;
            }
        }
        
        sort(pvec.begin(), pvec.end(), Point::sortId); // Sort the hits acording to the Id
        vector<int> tmpvec;
        for (int ip=0;ip<pvec.size();ip++) tmpvec.push_back(pvec[ip].id()); // Note the hit indices
        if (pvec.size() >= TRACKLET)
            tracklet.push_back(tmpvec);
        else
            shortpath.push_back(tmpvec);
    }

/*
    for (auto it1 = hitmap.begin(); it1 != hitmap.end(); n1++) {
        Point p0 = (*it1).second;
        int neighbour = p0.neighbour();
        vector<Point> pvec;
        while (neighbour != -1) {
            pvec.push_back(p0);
            auto it = hitmap.find(p0.id());
            if (it == hitmap.end()) { it1++; break; }
            it1 = hitmap.erase(it);
            p0 = points[neighbour];
            neighbour = p0.neighbour();
        }
        sort(pvec.begin(), pvec.end(), Point::sortId); // Sort the hits acording to the Id
        vector<int> tmpvec;
        for (int ip=0;ip<pvec.size();ip++) tmpvec.push_back(pvec[ip].id()); // Note the hit indices
        if (pvec.size() >= TRACKLET)
            tracklet.push_back(tmpvec);
        else
            shortpath.push_back(tmpvec);
    }
 */
/*
    // Search neighbouring hits, the neural network recall identifies the hit belonging to a tracklet
    cout << "Find tracklets..." << endl;
    int nd(0), nr(0), nn(0), np(0);
    for (auto it1 = hitmap.begin(); it1 != hitmap.end(); it1++, n1++) {
        if (VERBOSE && n1%1000==0) cout << n1 << " size " << points.size() << endl;
        Point p0 = (*it1).second;// it1->second; // Seeding point
        Point p1 = (*it1).second;// it1->second;
        vector<Point> pvec;
        pvec.push_back(p1); //Note the seeding point in the first place
        if (VERBOSE) cout << endl << p1.id() << "(0) ";
        int ntry = 0;
        auto tmp = it1;
        for (auto it2 = ++tmp; it2 != hitmap.end(); it2++) { //
            if (it1==it2) continue;
            Point p2 = (*it2).second;//it2->second;
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
                hitmap.erase(it2);  // Remove the corresponding point from the set
                *it2--;
                p0 = p1;// Note the assigned hits
                p1 = p2;
                continue;
            }
            else {
                if (VERBOSE) cout << p2.id() << "{" << (int) 100*recall << "} ";
                nn++;
            }
            
        }
        hitmap.erase(it1);  // Remove the corresponding point from the set
        *it1--;

        sort(pvec.begin(), pvec.end(), Point::sortId); // Sort the hits acording to the Id
        vector<int> tmpvec;
        for (int ip=0;ip<pvec.size();ip++) tmpvec.push_back(pvec[ip].id()); // Note the hit indices
        if (pvec.size() >= TRACKLET)
            tracklet.push_back(tmpvec);
        else
            shortpath.push_back(tmpvec);
    }
*/
    cout << endl << "Number of tracklets   : " << tracklet.size() << endl;
    cout << "Number of short tracks: " << shortpath.size() << endl;
    if (tracklet.size() == 0) exit(0);
    
    // Sort the tracklet vector according to the tracklet length
    
    sort(tracklet.begin(), tracklet.end(), sortFunc);
    
    // Print out the tracks vector
    if (VERBOSE) {
        cout << "Tracklets:" << endl;
        for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    }
    
    // Print out the shot tracks vector
    if (VERBOSE) {
        cout << "Short Tracklets:" << endl;
        for( int i=0; i<shortpath.size(); i++ ) print(shortpath[i]);
    }
    
    cout << endl << "Re-assign n/a hits to tracklets..." << endl;
    
    // Re-assign the short paths to tracklets
    int napoints = 0;
    int rapoints = 0;
    for (vector<vector<int> >::iterator it = shortpath.begin(); it != shortpath.end(); ++it) napoints += (*it).size();

#ifdef REASSIGN
    for (vector<vector<int> >::iterator it1 = tracklet.begin(); it1 != tracklet.end(); ++it1) {
        vector<int> &t1 = *it1;
        int seed1 = t1[t1.size()-1]; // last hit
        int seed2 = t1[t1.size()-2]; // 2nd last hit
        for (vector<vector<int> >::iterator it2 = shortpath.begin(); it2 != shortpath.end(); ++it2) {
            vector<int> &t2 = *it2;
            int seed3 = t2[0];
            double recall = Recall3(points[seed1],points[seed2],points[seed3])[0];
            if (recall>THRESHOLD) {
                for (int i=0;i<t2.size();i++) {
                    t1.push_back(t2[i]);
                    rapoints++;
                }
                shortpath.erase(it2);
                *it2--;
                break;
            }
        }
    }

    // Print out the tracks vector
    if (VERBOSE) {
        cout << "Tracklets after re-assignment:" << endl;
        for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    }

#endif
    
    // Assign labels
    int n = 0;
    int nc = 0;
    for (int i=0;i<tracklet.size();i++) {
        int id = 0;
        for (int j=0;j<tracklet[i].size();j++) {
            int hit = tracklet[i][j];
            if (j==0) id = hit;
            labels[hit] = i+1;
            p[hit].setval(i+1);
            if (hit == id++) nc ++;
            n++;
        }
    }
    
    /*
     // Gather the not assigned points
     for (int i=0;i<nhits;i++) {
     if (labels[i] == 0) {
     points.push_back(p[i]);
     }
     }
     */
    
    cout << "Number of hits                      : " << nhits << endl;
    cout << "Number of assigned points           : " << n << endl;
    cout << "Number of correctly assigned points : " << nc << endl;
    cout << "Number of unassigned points         : " << napoints << endl;
    cout << "Number of reassigned points : " << rapoints << endl;
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


