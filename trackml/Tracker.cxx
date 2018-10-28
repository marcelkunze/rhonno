// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"

// Assign track labels to hits (x,y,z)
// Hits are sorted according to their distance from the origin
// The hit pair quality is assessed by neural network function recall2
// Starting from the innermost hit a seeding hit pair is identified and extended to a tracklet
// The tracklet quality is assessed by neural network function recall3
// The tracklets are extended until no further suiting hit is found
// Corresponding labels are assigned to hits

#include "XMLP.h"
#include <ctime>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <set>

using namespace std;

//Structure for storing promising triples of hits
struct triple {
    int x, y, z;   // hit ids
    double r;      // quallity
    triple() {}
    triple(int a, int b, int c, double v=0) : x(a), y(b), z(c), r(v) {}
};

bool operator<(const triple&a, const triple&b) {
    return a.r < b.r;
}

bool operator==(const triple&a, const triple&b) {
    return a.x==b.x && a.y==b.y && a.z==b.z;
}

long Tracker::seedsok(0),Tracker::seedstotal(0);
long Tracker::trackletsok(0),Tracker::trackletstotal(0);
long Tracker::nr(0),Tracker::nd(0),Tracker::np(0),Tracker::nt(0),Tracker::nx(0);
long Tracker::n1(0),Tracker::n2(0),Tracker::n3(0),Tracker::n4(0);

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

// CHeck whether the points in a vector belong to the same track
long Tracker::checkLabels(std::vector<Point> &p) {
    if (p.size()==0) return 0;
    int label = p[0].label();
    if (label<=0) return 0;
    long n = 0;
    for (auto it:p) {
        if (it.label()==label) n++;
    }
    return n;
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
    if (angle<0.5*M_PI) {
        if (angle>DELTANN) { nx++; return null; } // Check 0 deg.
    }
    else {
        if ((M_PI-angle)>DELTANN) { nx++; return null; } // Check 180 deg.
    }

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

// Preselect points in a cylinder defined by RMIN..RMAX and ZMIN..ZMAX
long Tracker::selectPoints(std::vector<Point> &points,std::vector<Point> &inner,std::vector<Point> &outer,double rmin,double rmax,double zmin,double zmax) {
    long n = 0;
    for (auto it=points.begin(); it!=points.end();it++) {
        Point &p = *it;
        double z = p.z();
        double r = p.rz();
        if (z>ZMAX || z<ZMIN) { outer.push_back(p); continue; }
        if (r>RMAX || r<RMIN) { outer.push_back(p); continue; }
        inner.push_back(p);
        n++;
    }
    
    return n;
}

//id==49
#define TBD id==77
#define REF false

// Select points wrt. a reference point
long Tracker::selectPoints(std::vector<Point> &points, std::vector<Point> &good, std::vector<Point> &bad, Point &ref, double deltar, double deltathe, double distance)
{
    int id = ref.id();
    if (TBD&&REF) cout << "---------------------------> id: " << id << endl;

    int knn = 0;
    vector<Point> neighbours;
    auto it = points.begin();
    
    while (it++ != points.end()) {
        if (knn++ > MAXKNN) break; // Max. number of neighbouring points reached
 
        Point &p = *it;
        if (id==p.id()) continue;

        if (TBD&&REF) cout << p.id() << endl;
            
        double dr = abs(p.r()-ref.r()); // Check the radial distance
        if (dr > deltar) { if (TBD&&REF) cout << "R " << dr << endl; nr++; bad.push_back(p); continue; }
            
        double dt = abs(p.theta()-ref.theta()); // Check the elongation
        if (dt > deltathe) { if (TBD&&REF) cout << "T " << dt << endl; nt++; bad.push_back(p); continue; }
            
        double d = p.distance(ref); // Check the euclidean distance
        if (d > distance*p.r()) { if (TBD&&REF) cout << "D " << d << endl; nd++; bad.push_back(p); continue; }
            
        //double angle = acos(Point::dot(p0,p1)); // Check angle between (p0,p1)
        //if (angle > DELTAPHI) { if (TBD&&REF) cout << "P " << angle << endl; np++; continue; }
            
        good.push_back(p);
            
        if (TBD&&REF) {
                cout << p.id() << ": R " << dr << " T " << dt << " D " << d;
        }
            
    }
    
    if (TBD&&REF) cout << "<--------------------------- id: " << id << endl;
    
    return good.size();
}

// Look for seeding points using a KNN algorithm and a neural network to identify hit pairs
long Tracker::findSeeds(Point &p0, std::vector<Point> &points, std::vector<Point> &seed)
{
    if (points.size()<2) return 0;
    
    vector<Point> neighbours,bad;
    selectPoints(points,neighbours,bad,p0,DELTAR,DELTATHE,DISTANCE);
    
    sort(neighbours.begin(),neighbours.end(),Point::sortDist);
    
    // Generate seeding points
    for (auto it:neighbours)
    {
        Point &p1 = it;
        double recall = checkTracklet(p0,p1); // Search for hit pairs
        if (recall > 0) {
            if (REF) cout << p0.id() << " " << p1.id() << ": R2 OK " << recall << endl;
            p1.setrecall(recall);
            seed.push_back(p1);
        }
        else
            if (REF) cout << p0.id() << " " << p1.id() << ": R2 NOK " << recall << endl;
    }
    
    sort(seed.begin(),seed.end(),Point::sortRecall);
    
    int n = 0;
    for (auto &it:seed)
    {
        if (n>=NEIGHBOURS) break;
        Point &p1 = it;
        p0.setneighbour(p1.id(),n);
        p0.setrecall(p1.recall(),n);
        n++;
    }
    
    if (VERBOSE) { cout << "seed " << p0.id() << ": "; for (auto &it:seed) cout << it.id() << " " ; cout << endl; }
    
    long size = seed.size();
    seedstotal += seed.size();
    seedsok += checkLabels(seed);
    
    return size;
}

// Generate tracklets of 3 points wrt. the first point in seed
long Tracker::findTriples(Point &p0, std::vector<Point> &seed,std::vector<triple> &triples)
{
    long size = seed.size();
    if (size<2) return 0;

    triple t;
    t.x = p0.id();
    
    for (auto it1=seed.begin(); it1 != seed.end()-1; it1++)
    {
        Point &p1 = *it1;
        t.y = p1.id();
        for (auto it2=it1+1; it2 != seed.end(); it2++) {
            Point &p2 = *it2;
            double recall = checkTracklet(p0,p1,p2);
            if (recall > 0) {
                t.z = p2.id();
                t.r = recall;
                triples.push_back(t);
                if (REF) cout << t.x << " " << t.y << " " << t.z << ": R3 OK " << recall << endl;
            }
            else
                if (REF) cout << t.x << " " << t.y << " " << t.z << ": R3 NOK " << recall << endl;
            
        }
    }
    
    if (VERBOSE) { cout << "triples " << p0.id() << ":"; for (auto &it:triples) cout << " (" << it.x << "," << it.y << "," << it.z << ":" << it.r << ")"; cout << endl; }
    
    return triples.size();
}

double Tracker::checkTracklet(Point &p0,Point &p1)
{
    double recall = Recall2(p0,p1)[0];
    bool ok = recall>THRESHOLD2;
    recall = ok ? recall : -recall;
    if (ok) n4++;
    n2++;
    return recall;
}

double Tracker::checkTracklet(Point &p0,Point &p1, Point &p2)
{
    double recall = Recall3(p0,p1,p2)[0];
    bool ok = recall>THRESHOLD3;
    recall = ok ? recall : -recall;
    if (ok) n1++;
    n3++;
    return recall;
}

// Find tracks from points
int Tracker::findTracks(int nhits,float *x,float *y,float *z,int* labels,int *volumes,int *layers,int *modules,float *weights)
{
    std::clock_t c_start = std::clock();
    
    Point *p = new Point[nhits];
    vector<Point> points;
    points.reserve(nhits);
    map<int,Point*> hitmap;
    
    // Set up a cache for the point coordinates
    //cout << "Set up points cache..." << endl;
    for (int i=0;i<nhits;i++) {
        //labels[i] = 0;
        p[i] = Point(x[i],y[i],z[i],i,labels[i]);
        p[i].setvolume(volumes[i]);
        p[i].setlayer(layers[i]);
        p[i].setmodule(modules[i]);
        p[i].setweight(weights[i]);
        points.push_back(p[i]);
    }
    
    // Select points in a cylinder around the origin
    vector<Point> selection, other;
    long n = selectPoints(points,selection,other,RMIN,RMAX,ZMIN,ZMAX);
    if (n==0) return 0;
    
    // Sort the hits according to distance from origin
    //for (int i=0;i<nhits;i++) points[i] = p[i];
    cout << "Sorting " << nhits << " hits..." << endl;
    sort(selection.begin(),selection.end(),Point::sortRz);
    
    // Search neighbouring hits
    cout << "Searching neighbours..." << endl;
    
    for (auto &it: selection) {
        vector<Point> seed;
        vector<triple> triples;
        findSeeds(it,selection,seed);
        findTriples(it,seed,triples);
        // Transfer the result
        for (auto &it:seed) { // Note the neighbour ids and recall values
            int id = it.id();
            for (int j=0;j<NEIGHBOURS;j++) {
                p[id].setneighbour(it.neighbour(j),j);
                p[id].setrecall(it.recall(j),j);
            }
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
    for (auto it=selection.begin(); it!=selection.end();it++) {
        Point &p = *it;
        int id = p.id();
        hitmap[id] = &p;
    }
    
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
        if (pvec.size() >= TRACKLET) {
            tracklet.push_back(tmpvec);
            trackletstotal += pvec.size();
            trackletsok += checkLabels(pvec);
        }
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
        for (auto &it:tracklet) print(it);
    }
    
    // Print out the shot tracks vector
    if (VERBOSE) {
        cout << "Short Tracklets:" << endl;
        for (auto &it:shortpath) print(it);
    }
    
    cout << endl << "Re-assign n/a hits to tracklets..." << endl;
    
    // Re-assign the short paths to tracklets
    int napoints = 0;
    int rapoints = 0;
    for (auto &it:shortpath) napoints += it.size();
    
#ifdef REASSIGN
    for (auto &it1:tracklet) {
        vector<int> &t1 = it1;
        int seed1 = t1[t1.size()-1]; // last hit
        int seed2 = t1[t1.size()-2]; // 2nd last hit
        for (auto &it2:shortpath) {
            vector<int> &t2 = it2;
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
    cout << "Recalls3 : " << n3 << endl << endl;

    double perc = 0.0;
    if (seedstotal>0) perc = 100. * seedsok / (double) seedstotal;
    cout << "findSeeds: " << perc << "%" << endl;
    
    perc = 0.0;
    if (trackletstotal>0) perc = 100. * trackletsok / (double) trackletstotal;
    cout << "findTracks: " << perc << "%" << endl << endl;
    
    delete [] p;
    
    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    return (int) tracklet.size();
}
