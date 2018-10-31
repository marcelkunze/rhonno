// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

// Assign track labels to hits (x,y,z)
// Hits are sorted according to their distance from the origin
// The hit pair quality is assessed by neural network function recall2
// Starting from the innermost hit a seeding hit pair is identified and extended to a tracklet
// The tracklet quality is assessed by neural network function recall3
// The tracklets are extended until no further suiting hit is found
// Corresponding labels are assigned to hits

#include "Tracker.h"
#include "Point.h"
#include "XMLP.h"
#include <ctime>
#include <iostream>
#include <utility>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <set>

using namespace std;

bool operator<(const triple&a, const triple&b) {
    return a.r < b.r;
}

bool operator==(const triple&a, const triple&b) {
    return a.x==b.x && a.y==b.y && a.z==b.z;
}

//does hits a and b correspond to the same particle?
int Tracker::samepart(int a, int b) {
    long long aa = truth_part[a];
    long long bb = truth_part[b];
    return aa == bb && aa;
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
double* Tracker::recall2(Point &p1, Point &p2)
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
double* Tracker::recall3(Point &p1, Point &p2, Point &p3)
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
#define TBD true
#define REF true

// Select points wrt. a reference point
long Tracker::selectPoints(std::vector<Point> &points, std::vector<Point> &good, std::vector<Point> &bad, Point &ref, double deltar, double deltathe, double distance)
{
    int id = ref.id();
    if (TBD&&REF) cout << "---------------------------> id: " << id << endl;

    int knn = 0;
    vector<Point> neighbours;
    auto it = points.begin();
    
    while (it != points.end()) {
        if (knn++ > MAXKNN) break; // Max. number of neighbouring points reached
 
        Point &p = *it++;
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
    
    if (TBD&&REF) cout << endl << "<--------------------------- id: " << id << endl;
    
    return good.size();
}

// Look for seeding points using a KNN search and a neural network to identify hit pairs
std::vector<pair<int,float> > Tracker::findSeeds(Point &p0, std::vector<Point> &neighbours)
{
    vector<pair<int,float> > seed;
    
    paths.add(p0.id());
    
    // Generate seeding points
    for (auto it:neighbours)
    {
        Point &p1 = it;
        double recall = checkTracklet(p0,p1); // Search for hit pairs
        if (recall > 0) {
            if (REF) cout << p0.id() << " " << p1.id() << ": R2 OK " << recall << endl;
            seed.push_back(make_pair(p1.id(),(float)recall));
            paths.add(p0.id(),p1.id(),1000*recall);
            // Add double hits
            if (p0.twin()>0) {
                int twin = p0.twin();
                Point &t = points[twin];
                recall = checkTracklet(t,p1);
                seed.push_back(make_pair(twin,(float)recall));
                if (twin>p0.id())
                    paths.add(p0.id(),twin,1000*recall);
                else
                    paths.add(twin,p0.id(),1000*recall);

                if (REF) cout << " Added double hit " << twin << endl;
            }
        }
        else
            if (REF) cout << p0.id() << " " << p1.id() << ": R2 NOK " << recall << endl;
    }
    
    long size = seed.size();
    seedstotal += size;
    //seedsok += checkLabels(seed);
    
    return seed;
}

// Look for seeding points by hit pair combinations in the innnermost layers
void Tracker::findSeeds()
{
    const pair<int, int> start_list[5] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}};
    
    for (int i = 0; i < 5; i++) {
        int tube1 = start_list[i].first;
        for (auto &a : tubePoints[tube1]) {
            int tube2 = start_list[i].second;
            vector<Point> &b = tubePoints[tube2];
            std::vector<Point> neighbours,bad;
            selectPoints(b,neighbours,bad,a,DELTAR,DELTATHE,DISTANCE); // preselection of candidates b wrt a
            sort(neighbours.begin(),neighbours.end(),Point::sortDist);
            if (neighbours.size()>MAXKNN) neighbours.resize(MAXKNN); // k nearest neighbours
            vector<pair<int,float> > seed = findSeeds(a,neighbours);
            long n = seed.size();
            if (n>0&&VERBOSE) {
                cout << n << " seeds from " << a.id() << " (Tube: " << tube1 << "->" << tube2 << "):" << endl;
                for (auto it: seed) cout << it.first << "(" << it.second << ") ";
                cout << endl;
            }
        }
    }
    
    if (VERBOSE) {
        cout << paths << endl;
    }

}

//Find pairs using a neural network
vector<pair<int, int> > Tracker::findPairs() {
    
    const int n = 50;//How many pairs of layers to consider. Roughly proportional to run-time, and setting this to 30 gave practically the same score (less than 0.0002 reduction)
    pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
    
    vector<pair<int, int> > pairs;
    for (int i = 0; i < n; i++) {
        for (auto a : tube[start_list[i].first]) {
            for (auto b : tube[start_list[i].second]) {
                double recall = recall2(p[a],p[b])[0];
                if (recall > THRESHOLD2)
                    pairs.push_back(make_pair(a, b));
            }
        }
    }
    return pairs;
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
    double recall = recall2(p0,p1)[0];
    bool ok = recall>THRESHOLD2;
    recall = ok ? recall : -recall;
    if (ok) n4++;
    n2++;
    return recall;
}

double Tracker::checkTracklet(Point &p0,Point &p1, Point &p2)
{
    double recall = recall3(p0,p1,p2)[0];
    bool ok = recall>THRESHOLD3;
    recall = ok ? recall : -recall;
    if (ok) n1++;
    n3++;
    return recall;
}

// Generate an index to address the detector layers 0..47
int Tracker::getLayer(int volume_id, int layer_id) {
    
    int itopo[48] = {10,9,8,7,6,5,4,0,1,2,3,11,12,13,14,15,16,17,34,32,30,28,26,24,18,19,20,21,36,38,40,42,44,46,35,33,31,29,27,25,22,23,37,39,41,43,45,47};
    int metai_list[9][7] = {{0,1,2,3,4,5,6},{7,8,9,10,-1,-1,-1},{11,12,13,14,15,16,17},{18,19,20,21,22,23,-1},{24,25,26,27,-1,-1,-1},{28,29,30,31,32,33,-1},{34,35,36,37,38,39,-1},{40,41,-1,-1,-1,-1,-1},{42,43,44,45,46,47}};

    if (volume_id<7 || volume_id>18) return -1;
    if (volume_id <= 9)
        volume_id -= 7;
    else if (volume_id <= 14)
        volume_id -= 9;
    else
        volume_id -= 10;
    
    layer_id = layer_id/2-1;
    
    int index = itopo[metai_list[volume_id][layer_id]];
    return index;
}

void Tracker::readTubes() {
    
    for (int i = 0; i < points.size(); i++) {
        tube[points[i].layer()].push_back(i);
        //if (VERBOSE) cout << "Point " << i << " layer:" << points[i].layer() << endl;
    }
 
    for (int i = 0; i < 48; i++) {
        for (auto it : tube[i]) {
            tubePoints[i].push_back(points[it]);
            if (VERBOSE) {
                if (layer[i].type == Disc)
                    cout << "Disc " << i << ": Point " << it << " layer:" << points[it].layer() << endl;
                else
                    cout << "Tube " << i << ": Point " << it << " layer:" << points[it].layer() << endl;
            }
        }
    }

    for (int i = 0; i < 48; i++) {
        if (layer[i].type == Disc)
            sort(tube[i].begin(), tube[i].end(), z_cmp);
        else
            sort(tube[i].begin(), tube[i].end(), r_cmp);
    }
    
    for (int i = 0; i < 48; i++) {
        if (layer[i].type == Tube)
            sort(tubePoints[i].begin(), tubePoints[i].end(), Point::sortZ);
        else
            sort(tubePoints[i].begin(), tubePoints[i].end(), Point::sortRz);
    }

    // Filter double hits
    for (int i = 0; i < 48; i++) {
        vector<Point> &pvec = tubePoints[i];
        auto it = pvec.begin();
        while (it++ != pvec.end()) {
            Point &p0 = *(it-1);
            int id0 = p0.id();
            Point &p1 = *it;
            int id1 = p1.id();
            double d = p0.distance(p1);
            if (d<0.005) {
                p0.settwin(id1);
                ntwins++;
                if (VERBOSE) cout << "Twin " << id0 << "," << id1 << ":" << d << endl;
            }
        }
    }

}

bool Tracker::z_cmp(const int a, const int&b) { return points[a].z() < points[b].z(); }
bool Tracker::r_cmp(const int&a, const int&b) { return points[a].rz() < points[b].rz(); }

// Find tracks from points
int Tracker::findTracks(int nhits,float *x,float *y,float *z,int* layers,int* labels,int *truth)
{
    std::clock_t c_start = std::clock();
    
    p = new Point[nhits];
    points.reserve(nhits);
    map<int,Point*> hitmap;
    
    // Set up a cache for the point coordinates
    //cout << "Set up points cache..." << endl;
    for (int i=0;i<nhits;i++) {
        //labels[i] = 0;
        p[i] = Point(x[i],y[i],z[i],i,labels[i],truth[i]);
        p[i].setlayer(layers[i]);
        points.push_back(p[i]);
    }
    
    readTubes();
   
    // Search neighbouring hits
    cout << "Searching seeds..." << endl;
    findSeeds();
    
    // transfer the results
    for (auto &ip : points) {
        int a = ip.id();
        int n = 0;
        map<int,int> const &path = paths.connections(a);
        if (path.size()==0) continue;
        for (auto &it : path ) {
            if (n>NEIGHBOURS) break;
            int id = it.first;
            float recall = 0.001*it.second;
            ip.setneighbour(id,n);
            ip.setrecall(recall,n);
            p[a].setneighbour(id,n);
            p[a].setrecall(recall,n);
            n++;
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
    
    // Analyze the graph

    // fill the hitmap
    for (auto it=points.begin(); it!=points.end();it++) {
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
    cout << "Number of double hits               : " << ntwins << endl;
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

// functions to read trackml data

void Tracker::readBlacklist(string base_path,int filenum) {
    if (filenum < 1000) return;
    char file[1000];
    sprintf(file, "%s/event%09d-blacklist_particles.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open blacklist\n"); return; }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    long long particle_id;
    while (fscanf(fp, "%lld", &particle_id) == 1) {
        blacklist.insert(particle_id);
    }
    fclose(fp);
}

void Tracker::readTruth(string base_path,int filenum) {
    if (filenum < 1000) return;
    char file[1000];
    sprintf(file, "%s/event%09d-truth.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open truth\n"); return; }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    while (1) {
        int hit_id;
        long long particle_id;
        double tx, ty, tz, tpx, tpy, tpz, weight;
        if (fscanf(fp, "%d,%lld,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &hit_id, &particle_id, &tx, &ty, &tz, &tpx, &tpy, &tpz, &weight) != 9) break;
        if (!particle_id || blacklist.count(particle_id)) continue;
        
        truth_tracks[particle_id].push_back(hit_id);
        truth_pos[hit_id] = point(tx, ty, tz);
        truth_mom[hit_id] = point(tpx, tpy, tpz)*Bfield;
        truth_weight[hit_id] = weight;
        part_weight[particle_id] += weight;
        truth_part[hit_id] = particle_id;
        
        map<long long, int>::iterator it = partIDmap.find(particle_id);
        if( it==partIDmap.end() ){
            cout<<"Particle ID not found in map!!!"<<endl;
            cout<<"ID= "<<hit_id<<" hit "<<hit_id<<" iterator at ID "<<it->first<<endl;
            exit(0);
        }
        
        int newID = it->second;
        if( newID < 0 || newID>= particles.size() ){
            cout<<"Mapped particle ID is wrong!!!"<<endl;
            cout<<"ID= "<<hit_id<<" new ID "<<newID<<endl;
            exit(0);
        }

        Particle &p = particles[newID];
        p.hit.push_back(hit_id);

    }
    fclose(fp);
    
    cout << truth_tracks.size() << " particles with truth" << endl;
}

void Tracker::readParticles(string base_path,int filenum) {
    char file[1000];
    sprintf(file, "%s/event%09d-particles.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open particles\n"); return; }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    while (1) {
        long long id;
        int type;
        point p, m;
        int hits;
        int q;
        if (fscanf(fp, "%lld,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d", &id, &type, &p.x, &p.y, &p.z, &m.x, &m.y, &m.z, &q, &hits) == -1) break;
        if (blacklist.count(id)) continue;
        start_pos[id] = p;
        start_mom[id] = m*Bfield;
        part_q[id] = -q;
        part_hits[id] = hits;
        
        Particle part;
        partIDmap[ (long long) id ] = (int)particles.size();
        part.id = id;
        part.type = type;
        part.x = p.x;
        part.y = p.y;
        part.z = p.z;
        part.r = 0;
        part.px = m.x;
        part.py = m.y;
        part.pz = m.z;
        part.q = q;
        part.hits = hits;
        particles.push_back(part);
    }
    fclose(fp);
    cout << start_pos.size() << " particles" << endl;
}


bool Tracker::track_cmp(int a, int b) {
    point&ma = truth_mom[a];
    point&mb = truth_mom[b];
    double va = ma*ma, vb = mb*mb;
    if (fabs((va-vb)/va) > 1e-5) return va > vb;
    return (truth_pos[b]-truth_pos[a])*ma > 0;
}

//Sort the hits in each track chronologically
void Tracker::sortTracks() {
    int fails = 0, goods = 0;
    for (auto &p : truth_tracks) {
        vector<int> v;
        for (int hit_id : p.second) {
            v.push_back(hit_id);
        }
        sort(v.begin(), v.end(), track_cmp);
        for (int i = 0; i < v.size(); i++)
            p.second[i] = v[i];
        int bad = 0;
        for (int i = 2; i < v.size(); i++) {
            point&a = truth_pos[p.second[i-2]];
            point&b = truth_pos[p.second[i-1]];
            point&c = truth_pos[p.second[i]];
            if ((c.z-b.z)*(b.z-a.z) < 0) {
                fails++;
                bad++;
            }
            else goods++;
        }
    }
}


void Tracker::initOrder() {
    int handpicked[48] = {};
    int c = 0;
    for (int i = 0; i < 4; i++) handpicked[c++] = 7+i;
    for (int i = 0; i < 7; i++) handpicked[c++] = 7-1-i;
    for (int i = 0; i < 7; i++) handpicked[c++] = 11+i;
    
    for (int i = 0; i < 4; i++) handpicked[c++] = 24+i;
    for (int i = 0; i < 2; i++) handpicked[c++] = 40+i;
    for (int i = 0; i < 6; i++) {
        handpicked[c++] = 24-1-i;
        handpicked[c++] = 40-1-i;
    }
    
    for (int i = 0; i < 6; i++) {
        handpicked[c++] = 28+i;
        handpicked[c++] = 42+i;
    }
    for (int i = 0; i < 48; i++) {
        topo[i] = handpicked[i];
        itopo[topo[i]] = i;
        //cout << itopo[i] << "," ;
    }
}

//init layer geometries
void Tracker::initLayers() {
    double avgz1[2][7] = {{-1500,-1300,-1100,-960,-820,-700,-600},
        { 600, 700, 820, 960, 1100, 1300, 1500}};
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < 7; i++) {
            layer[k*11+i].minr = 30;
            layer[k*11+i].maxr = 176.5;
            layer[k*11+i].avgz = avgz1[k][i];
            layer[k*11+i].type = Disc;
        }
    double avgz2[2][6] = {{-2950,-2550,-2150,-1800,-1500,-1220},
        { 1220, 1500, 1800, 2150, 2550, 2950}};
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < 6; i++) {
            layer[k*10+i+18].minr = 240;
            layer[k*10+i+18].maxr = 701;
            layer[k*10+i+18].avgz = avgz2[k][i];
            layer[k*10+i+18].type = Disc;
            
            layer[k*8+i+34].minr = 755;
            layer[k*8+i+34].maxr = 1018;
            layer[k*8+i+34].avgz = avgz2[k][i];
            layer[k*8+i+34].type = Disc;
        }
    
    double avgr1[4] = {32.3, 72.1, 116.1, 172.1};
    double avgr2[4] = {260.3, 360.2, 500.2, 660.2};
    double avgr3[2] = {820.2, 1020.2};
    
    for (int i = 0; i < 4; i++) {
        layer[i+7].minz =-491;
        layer[i+7].maxz = 491;
        layer[i+7].avgr = avgr1[i];
        layer[i+7].type = Tube;
    }
    for (int i = 0; i < 4; i++) {
        layer[i+24].minz =-1084;
        layer[i+24].maxz = 1084;
        layer[i+24].avgr = avgr2[i];
        layer[i+24].type = Tube;
    }
    for (int i = 0; i < 2; i++) {
        layer[i+40].minz =-1084;
        layer[i+40].maxz = 1084;
        layer[i+40].avgr = avgr3[i];
        layer[i+40].type = Tube;
    }
    Layer layer2[48];
    for (int i = 0; i < 48; i++) layer2[i] = layer[i];
    for (int i = 0; i < 48; i++) layer[i] = layer2[topo[i]];
    
    layer[0].var0 = 1e-3;
    layer[1].var0 = 5e-4;
    for (int i = 2; i < 18; i++) layer[i].var0 = 3e-4;
    for (int i = 18; i < 22; i++) layer[i].var0 = 5e-2;
    for (int i = 22; i < 48; i++) layer[i].var0 = i%2 || i == 22 ? 9 : 0.1;
    
    for (int i = 0; i < 4; i++) layer[i].var1 = 0.5;
    for (int i = 4; i < 18; i++) layer[i].var1 = 5;
    for (int i = 18; i < 24; i++) layer[i].var1 = 7;
    for (int i = 24; i < 48; i++) layer[i].var1 = i%2 ? 19 : 11;
}

void Tracker::readHits(string base_path, int filenum) {
    initOrder();
    
    char file[1000];
    if (filenum >= 1000)
        sprintf(file, "%s/event%09d-hits.csv", base_path.c_str(), filenum);
    else
        sprintf(file, "%s/test/event%09d-hits.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open hits\n");
        exit(1);
    }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    
    // For one indexing
    hits.push_back(point(0,0,0));
    polar.push_back(point(0,0,0));
    meta.push_back(point(0,0,0));
    metai.push_back(0);
    
    int layers[9] = {7,4,7,6,4,6,6,2,6};
    int metai_list[9][7];
    int c = 0;
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < layers[i]; j++)
            metai_list[i][j] = c++;
    }
    cout << "Detectors: " << c << endl;
    
    for (int i = 0; i < 48; i++) {
        layer[i].minr = layer[i].minz = 1e9;
        layer[i].maxr = layer[i].maxz =-1e9;
    }
    
    while (1) {
        long long hit_id;
        double tx, ty, tz;
        int volume_id, layer_id, module_id;
        if (fscanf(fp, "%lld,%lf,%lf,%lf,%d,%d,%d", &hit_id, &tx, &ty, &tz, &volume_id, &layer_id, &module_id) == -1) break;
        if (hit_id != hits.size()) cout << "Hit id's not as expected" << endl;
        meta.push_back(point(volume_id, layer_id, module_id));
        if (volume_id <= 9)
            volume_id -= 7;
        else if (volume_id <= 14)
            volume_id -= 9;
        else
            volume_id -= 10;
        
        int mi = itopo[metai_list[volume_id][layer_id/2-1]];
        
        hits.push_back(point(tx, ty, tz));
        track_hits[hit_id] = point(tx, ty, tz);
        polar.push_back(point(sqrt(tx*tx+ty*ty), atan2(ty,tx), tz));
        metai.push_back(mi);
        
        double r = sqrt(tx*tx+ty*ty);
        Layer&l = layer[metai_list[volume_id][layer_id/2-1]];
        l.minr = min(l.minr, r);
        l.avgr += r;
        l.maxr = max(l.maxr, r);
        l.minz = min(l.minz, tz);
        l.avgz += tz;
        l.maxz = max(l.maxz, tz);
        l.count++;
        //cerr << tz << ' ' << r << endl;
    }
    fclose(fp);
    cout << hits.size() << " hits" << endl;
    
    initLayers();
    
    for (int hit_id = 1; hit_id < hits.size(); hit_id++) {
        metai_weight[truth_part[hit_id]][metai[hit_id]] += truth_weight[hit_id];
    }
    
    map<double, double> mir[48], mar[48];
    for (int i = 1; i < hits.size(); i++) {
        int mi = metai[i];
        if (layer[mi].type != Disc) continue;
        double&mir_ = mir[mi][polar[i].z];
        double&mar_ = mar[mi][polar[i].z];
        if (!mir_) mir_ = 1e9;
        mir_ = min(mir_, polar[i].x);
        mar_ = max(mar_, polar[i].x);
    }
    
    map<double, int> zi[48];
    for (int mi = 0; mi < 48; mi++) {
        if (layer[mi].type != Disc) continue;
        int k = 0;
        for (auto p : mir[mi]) {
            double mir_ = mir[mi][p.first]-1e-5;
            double mar_ = mar[mi][p.first]+1e-5;
            z_minr[mi][k] = mir_;
            z_maxr[mi][k] = mar_;
            disc_z[mi][k] = p.first;
            zi[mi][p.first] = k++;
        }
    }
    
    metaz.resize(hits.size());
    metaz[0] = 0;
    for (int i = 1; i < hits.size(); i++) {
        int mi = metai[i];
        metaz[i] = meta[i].z;
        if (layer[mi].type == Disc)
            metaz[i] = zi[mi][hits[i].z];
    }

}

// Volumes: 7,8,9, 12,13,14, 16,17,18
// Volume 8 is innermost, it has modules 2,4,6,8 concentric cylinders outwards

void Tracker::readDetectors(string base_path) {
    char file[1000];
    sprintf(file, "%s/detectors.csv", base_path.c_str());
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open detectors\n");
        exit(1);
    }
    
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    
    Detector d;
    while (fscanf(fp, "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &d.volume_id, &d.layer_id, &d.module_id, &d.c.x, &d.c.y, &d.c.z,
                  // &d.rx.x, &d.rx.y, &d.rx.z,
                  // &d.ry.x, &d.ry.y, &d.ry.z,
                  // &d.rz.x, &d.rz.y, &d.rz.z,
                  &d.rx.x, &d.ry.x, &d.rz.x,
                  &d.rx.y, &d.ry.y, &d.rz.y,
                  &d.rx.z, &d.ry.z, &d.rz.z,
                  &d.d, &d.minw, &d.maxw, &d.h, &d.cell_w, &d.cell_h) == 21) {
        //if (d.module_id >= 10000 || d.layer_id >= 1000 || d.volume_id >= 100) cout << "What!?" << endl;
        detectors[d.volume_id*10000000+d.layer_id*10000+d.module_id] = d;
    }
    fclose(fp);
}

void Tracker::readCells(string base_path,int filenum) {
    char file[1000];
    if (filenum >= 1000)
        sprintf(file, "%s/event%09d-cells.csv", base_path.c_str(), filenum);
    else
        sprintf(file, "%s/test/event%09d-cells.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open cells\n");
        exit(1);
    }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    int hit_id, ch0, ch1;
    double value;
    while (fscanf(fp, "%d,%d,%d,%lf", &hit_id, &ch0, &ch1, &value) == 4) {
        hit_cells[hit_id].push_back(make_pair(make_pair(ch0, ch1), value));
        //cout << hit_id << ' ' << ch0 << ' ' << ch1 << ' ' << value << endl;
    }
    fclose(fp);
}

point Tracker::normalize(point a) {
    point ret = a*(1./dist(a));
    if (ret.z < 0) ret = ret*-1;
    return ret;
}

//Calculate direction of each hit with cell's data
void Tracker::initHitDir() {
    for (int hit_id = 1; hit_id < hits.size(); hit_id++) {
        point m = meta[hit_id];
        Detector&d = detectors[int(m.x)*10000000+int(m.y)*10000+int(m.z)];
        
        //if (!hit_cells[hit_id].size()) cout << "Hit with zero cells" << endl;
        //if (metai[hit_id] < 18) continue;
        
        //Use linear regression for direction
        double mx = 0, my = 0, mw = 0;
        auto&cells = hit_cells[hit_id];
        for (auto&c : cells) {
            double w = c.second;
            double x = c.first.first*d.cell_w;
            double y = c.first.second*d.cell_h;
            mw += w;
            mx += x*w;
            my += y*w;
        }
        mx /= mw;
        my /= mw;
        double mxx = 0, mxy = 0, myy = 0;
        for (auto&c : cells) {
            double w = c.second;
            double x = c.first.first*d.cell_w-mx;
            double y = c.first.second*d.cell_h-my;
            mxx += x*x*w;
            myy += y*y*w;
            mxy += x*y*w;
        }
        //Find eigenvector with minimum eigenvalue
        double a = mxx-myy, b = 2*mxy;
        double x = a+b+sqrt(a*a+b*b);
        double y =-a+b+sqrt(a*a+b*b);

        //Analytical formula for z
        double z = 2*d.d*(fabs(x)/d.cell_w+fabs(y)/d.cell_h+1e-8);
        x *= (cells.size()*1.-1.3);//1.3 != 1 was adjusted empirically
        y *= (cells.size()*1.-1.3);
        point d1(x,y,z), d2(x,y,-z);
        d1 = d.rx*d1.x+d.ry*d1.y+d.rz*d1.z;
        d2 = d.rx*d2.x+d.ry*d2.y+d.rz*d2.z;
        hit_dir[hit_id][0] = normalize(d1);
        hit_dir[hit_id][1] = normalize(d2);
        
        for (int k = 0; k < 2; k++)
            if (hit_dir[hit_id][k]*hits[hit_id] < 0)
                hit_dir[hit_id][k] = hit_dir[hit_id][k]*-1;
        
    }
}

// Data initialization
Point* Tracker::p;
vector<Point> Tracker::points;
digraph<int> Tracker::paths;
long Tracker::seedsok(0),Tracker::seedstotal(0);
long Tracker::trackletsok(0),Tracker::trackletstotal(0);
unsigned long Tracker::nr(0),Tracker::nd(0),Tracker::np(0),Tracker::nt(0),Tracker::nx(0);
unsigned long Tracker::n1(0),Tracker::n2(0),Tracker::n3(0),Tracker::n4(0),Tracker::ntwins(0);
vector<point> Tracker::hits; //hit position
vector<Particle> Tracker::particles; //true tracks
map<long long,int> Tracker::partIDmap; // create particle ID->index map

vector<int> Tracker::tube[48]; // List of hits in each layer
vector<Point> Tracker::tubePoints[48]; // List of hits in each layer
map<long long, vector<int> > Tracker::truth_tracks; //truth hit ids in each track
map<long long, point> Tracker::track_hits; // Find points in hits
int Tracker::assignment[150000];
point Tracker::truth_pos[150000], Tracker::truth_mom[150000]; //truth position and momentum
double Tracker::truth_weight[150000]; //weighting of each hit
long long Tracker::truth_part[150000]; //particle this hit belongs to
set<long long> Tracker::blacklist;
map<long long, double> Tracker::part_weight; //weighting of each particle
map<long long, map<int, double> > Tracker::metai_weight; //weighting of each particle hit, also adding duplicates
map<long long, point> Tracker::start_pos; //start position
map<long long, point> Tracker::start_mom; //start momentum
map<long long, int> Tracker::part_q; //start charge
map<long long, int> Tracker::part_hits; // = truth_tracks[particle_id].size()
int Tracker::topo[48], Tracker::itopo[48]; //reordering of layers for approximate sorting
vector<point> Tracker::polar; //hit position in polar / cylindrical coordinates
vector<point> Tracker::meta; //volume_id / layer_id / module_id
vector<int> Tracker::metai, Tracker::metaz; //ordered layer id in [0,48), and classification of z for disc layers in [0,4)
double Tracker::disc_z[48][4];
Layer Tracker::layer[48];
double Tracker::z_minr[48][4], Tracker::z_maxr[48][4];
map<int, Detector> Tracker::detectors;
vector<std::pair<pair<int, int>, double> > Tracker::hit_cells[150000]; //pair<pair<ch0, ch1>, value>
point Tracker::hit_dir[150000][2]; //The two possible directions of the hit according to the cell's data for each hit

