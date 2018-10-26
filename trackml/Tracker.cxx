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
#include <set>

using namespace std;

Point::Point(double x, double y, double z, int id, int label, int truth)
{
    _id = id;
    _label = label;
    _truth = truth;
    _x = x;
    _y = y;
    _z = z;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
    for (int i=0; i<NEIGHBOURS;i++) _neighbour[i] = _recall[i] = -1;
}

Point::Point(float x, float y, float z, int id, int label, int truth)
{
    _id = id;
    _label = label;
    _truth = truth;
    _x = x;
    _y = y;
    _z = z;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
    for (int i=0; i<NEIGHBOURS;i++) _neighbour[i] = _recall[i] = -1;
}

Point::Point(const Point &p)
{
    _id = p._id;
    _label = p._label;
    _truth = p._truth;
    _x = p._x;
    _y = p._y;
    _z = p._z;
    _r = p._r;
    _rz = p._rz;
    _phi = p._phi;
    _theta = p._theta;
    _distance = p._distance;
    for (int i=0; i<NEIGHBOURS;i++) {
        _neighbour[i] = p._neighbour[i];
        _recall[i] = p._recall[i];
    }
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

    x[0]     = p1.rz();     // r1
    x[1]     = p1.phi();    // phi1
    x[2]     = p1.z();      // z1
    x[3]     = p2.rz();     // r2
    x[4]     = p2.phi();    // phi2
    x[5]     = p2.z();      // z2
    
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
    if (angle > DELTANN) { nx++; return null; }
    
    x[0]     = p1.rz();     // r1
    x[1]     = p1.phi();    // phi1
    x[2]     = p1.z();      // z1
    x[3]     = p2.rz();     // r2
    x[4]     = p2.phi();    // phi2
    x[5]     = p2.z();      // z2
    x[6]     = p3.rz();     // r3
    x[7]     = p3.phi();    // phi3
    x[8]     = p3.z();      // z3
    
    return net.Recallstep(x);
}

//id==49
#define TBD id==20
#define REF true

void Tracker::kNearestNeighbour(std::vector<Point> &points)
{
    nr = nd = np = nt = nx = n1 = n2 = n3 = n4 =0;
    
    if (points.size()<2) return;
    
    for (auto it=points.begin(); it!=points.end();it++) {
        
        //if (i%10000==0) cout << i << endl;
        
        auto it2 = points.begin(); //it+1;
        
        Point &p0 = *it; // Seeding hit
        int id = p0.id();
        if (TBD&&REF) cout << "---------------------------> id: " << id << endl;
        
        int knn = 0;
        vector<Point> neighbours;
        while (it2++ != points.end()) {
            if (it==it2) continue;
            if (knn++ > MAXKNN) break; // Max. number of neighbouring points reached
            
            Point &p1 = *it2;
            if (TBD&&REF) cout << p1.id() << endl;
            
            double deltar = abs(p1.r()-p0.r()); // Check the radial distance
            if (deltar > DELTAR) { if (TBD&&REF) cout << "R " << deltar << endl; nr++; continue; }
            
            double deltat = abs(p1.theta()-p0.theta()); // Check the elongation
            if (deltat > DELTATHE) { if (TBD&&REF) cout << "T " << deltat << endl; nt++; continue; }
            
            double d = p1.distance(p0); // Check the euclidean distance
            if (d > DISTANCE*p0.r()) { if (TBD&&REF) cout << "D " << d << endl; nd++; continue; }
            
            //double angle = acos(Point::dot(p0,p1)); // Check angle between (p0,p1)
            //if (angle > DELTAPHI) { if (TBD&&REF) cout << "P " << angle << endl; np++; continue; }
            
            neighbours.push_back(p1);
            
            if (TBD&&REF) {
                cout << p1.id() << ": R " << deltar << " T " << deltat << " D " << d;
            }
            
        }
        if (TBD&&REF) cout << "<--------------------------- id: " << id << endl;
        
        //sort(neighbours.begin(),neighbours.end(),Point::sortDist);
        
        // Generate seeding points
        vector<Point> seed;
        for (auto it=neighbours.begin(); it != neighbours.end(); it++)
        {
            Point &p1 = *it;
            bool ok = checkTracklet(p0,p1);
            if (ok) {
                if (TBD) cout << p1.id() << ": R2 OK " << p1.recall() << endl;
                seed.push_back(p1);
            }
        }
        
        long size = seed.size();
        if (size==0) continue;
        if (size>NEIGHBOURS) size = NEIGHBOURS;
        
        int n = 0;
        for (auto it = seed.begin(); it!=seed.end(); it++, n++) { // Note the seed ids and recall values
            Point &p1 = *it;
            p0.setneighbour(p1.id(),n);
            p0.setrecall(p1.recall(),n);
        }
        
        if (VERBOSE) cout << "seed " << p0.id() << ":"; for (auto it=seed.begin(); it != seed.end(); it++) cout << it->id() << " " ; cout << endl;
        
        // Generate tracklets of 3 points
        if (size<2) continue;
        vector<Point> nextneighbours;
        for (auto it1=seed.begin(); it1 != seed.end()-1; it1++)
        {
            Point &p1 = *it1;
            nextneighbours.push_back(p1);
            for (auto it2=it1+1; it2 != seed.end(); it2++) {
                Point &p2 = *it2;
                bool ok = checkTracklet(p0,p1,p2);
                if (ok) {
                    if (REF) cout << p2.id() << ": R3 OK " << p2.recall() << endl;
                    nextneighbours.push_back(p2);
                }
                else
                    if (REF) cout << p2.id() << ": R3 NOK " << p2.recall() << endl;
                
            }
        }
        
        if (VERBOSE) cout << "nextneighbours " << p0.id() << ":"; for (auto it=nextneighbours.begin(); it != nextneighbours.end(); it++) cout << it->id() << " " ; cout << endl;
        
        sort(nextneighbours.begin(),nextneighbours.end(),Point::sortId);
        set<Point> s(nextneighbours.begin(),nextneighbours.end()); // Remove double entries
        
        if (VERBOSE) cout << "set " << p0.id() << ":"; for (auto it=s.begin(); it != s.end(); it++) cout << it->id() << " " ; cout << endl;
        
        int j = 0;
        for (auto it = seed.begin(); it!=seed.end(); it++, j++) { // Note the seed ids and recall values
            if (j>=NEIGHBOURS) break;
            if (it->recall()<=0.0) continue;
            p0.setneighbour(it->id(),j);
            p0.setrecall(it->recall(),j);
        }
        
    }
}

bool Tracker::checkTracklet(Point &p0,Point &p1)
{
    p1.setrecall(0.0);
    double recall = Recall2(p0,p1)[0];
    bool ok = recall>THRESHOLD2;
    if (ok) {
        p1.setrecall(recall);
        n4++;
    }
    n2++;
    return ok;
}

bool Tracker::checkTracklet(Point &p0,Point &p1, Point &p2)
{
    p2.setrecall(0.0);
    double recall = Recall3(p0,p1,p2)[0];
    bool ok = recall>THRESHOLD3;
    if (ok) {
        p2.setrecall(recall);
        n1++;
    }
    n3++;
    return ok;
}

int Tracker::findTracks(int nhits, float *x, float *y, float *z, int* labels)
{
    std::clock_t c_start = std::clock();
    
    Point *p = new Point[nhits];
    vector<Point> points;
    points.reserve(nhits);
    map<int,Point*> hitmap;

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
    sort(points.begin(),points.end(),Point::sortRz);
    
    // Search neighbouring hits
    cout << "Searching neighbours..." << endl;
    
    kNearestNeighbour(points);

    for (int i=0;i<nhits;i++) { // Note the neighbour ids and recall values
        for (int j=0;j<NEIGHBOURS;j++) {
            p[points[i].id()].setneighbour(points[i].neighbour(j),j);
            p[points[i].id()].setrecall(points[i].recall(j),j);
        }
    }
 
    if (VERBOSE) {
        for (int i=0;i<nhits;i++) {
            cout << p[i].id() << "(" ;
            for (int j=0;j<NEIGHBOURS;j++) cout << p[i].neighbour(j) << "/" << p[i].recall(j) << " ";
            cout << ") " << endl;
        }
        cout << endl;
    }
 

    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;

    // fill the hitmap
    for (int i=0;i<nhits;i++) hitmap[points[i].id()] = &points[i];

    // Swimmer
    cout << "Find tracklets..." << endl;
    vector<vector<int> > tracklet;
    vector<vector<int> > shortpath;
    
    while (hitmap.size()>0) {
        int n = 0;
        int neighbour = -1;
        vector<Point> pvec;
        auto it = hitmap.begin();
        Point *p0 = it->second;
        pvec.push_back(*p0);
        hitmap.erase(it++);
        while (it != hitmap.end()) { // Follow the path until there are no more neighbours
            neighbour = p0->neighbour(n);
            if (VERBOSE) cout << p0->id() << "->" << neighbour << endl;
            if (neighbour < 0 || neighbour >= nhits) break;
            auto it = hitmap.find(neighbour);
            if (n<NEIGHBOURS-1 && it==hitmap.end()) { // hit is already assigned
                if (VERBOSE) cout <<  "->" << neighbour << endl;
                n++;  // try an alternative neighbour
                continue;
            }
            if (it==hitmap.end()) break;
            Point *p1 = it->second;  // copy the point into the tracklet vector
            pvec.push_back(*p1);
            hitmap.erase(it++);
            n = 0;
            p0 = p1;
        }
        
        sort(pvec.begin(), pvec.end(), Point::sortId); // Sort the hits acording to the Id
        vector<int> tmpvec;
        for (auto ip=pvec.begin();ip!=pvec.end();ip++) tmpvec.push_back(ip->id()); // Note the hit indices
        if (VERBOSE) print(tmpvec);
        if (pvec.size() >= TRACKLET)
            tracklet.push_back(tmpvec);
        else
            shortpath.push_back(tmpvec);
        if (VERBOSE) cout << "---" << endl;
    }

    cout << endl << "Number of tracklets   : " << tracklet.size() << endl;
    cout << "Number of short tracks: " << shortpath.size() << endl;
    if (tracklet.size() == 0) exit(0);
    
    // Sort the tracklet vector according to the tracklet length
    
    //sort(tracklet.begin(), tracklet.end(), sortFunc);
    
    // Print out the tracks vector
    if (VERBOSE) {
        cout << "Tracklets:" << endl;
        for (auto it = tracklet.begin(); it != tracklet.end(); ++it) print(*it);
    }
    
    // Print out the shot tracks vector
    if (VERBOSE) {
        cout << "Short Tracklets:" << endl;
        for (auto it = shortpath.begin(); it != shortpath.end(); ++it) print(*it);
    }
    
    cout << endl << "Re-assign n/a hits to tracklets..." << endl;
    
    // Re-assign the short paths to tracklets
    int napoints = 0;
    int rapoints = 0;
    for (auto it = shortpath.begin(); it != shortpath.end(); ++it) napoints += (*it).size();

#ifdef REASSIGN
    for (auto it1 = tracklet.begin(); it1 != tracklet.end(); ++it1) {
        vector<int> &t1 = *it1;
        int seed1 = t1[t1.size()-1]; // last hit
        int seed2 = t1[t1.size()-2]; // 2nd last hit
        for (auto it2 = shortpath.begin(); it2 != shortpath.end(); ++it2) {
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
        for (auto it = tracklet.begin(); it != tracklet.end(); ++it) print(*it);
    }

#endif
    
    // Assign labels
    int na = 0;
    int nc = 0;
    for (int i=0;i<tracklet.size();i++) {
        int id = 0;
        for (int j=0;j<tracklet[i].size();j++) {
            int hit = tracklet[i][j];
            if (j==0) id = hit;
            labels[hit] = i+1;
            p[hit].setlabel(i+1);
            if (hit == id++) nc ++;
            na++;
        }
    }
    
    cout << "Number of hits                      : " << nhits << endl;
    cout << "Number of assigned points           : " << na << endl;
    cout << "Number of correctly assigned points : " << nc << endl;
    cout << "Number of unassigned points         : " << napoints << endl;
    cout << "Number of reassigned points : " << rapoints << endl;
    cout << "Distance <" << DISTANCE << ": " << nd <<endl;
    cout << "Radius   <" << DELTAR << ": " << nr <<endl;
    cout << "Theta    <" << DELTATHE << ": " << nt <<endl;
    cout << "Phi      <" << DELTAPHI << ": " << np <<endl;
    cout << "PhiNN    <" << DELTANN << ": " << nx <<endl;
    cout << "Threshold2 <" << THRESHOLD2 << ": " << n4 <<endl;
    cout << "Threshold3 <" << THRESHOLD3 << ": " << n1 <<endl;
    cout << "Recalls2 : " << n2 << endl;
    cout << "Recalls3 : " << n3 << endl;

    delete [] p;
    
    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    return (int) tracklet.size();
}
