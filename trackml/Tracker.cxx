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
#include "PolarModule.h"

#include <ctime>
#include <utility>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

void loadAdjacent(string,string);
void initDensity3(void);

// Find tracks from points
int Tracker::findTracks(int nhits,float *x,float *y,float *z,float *cx,float *cy,float *cz,int* volume,int* layer,int* module,int *hitid,long long *trackid,int* label)
{
    std::clock_t c_start = std::clock();
    
    _x = x;
    _y = y;
    _z = z;
    _cx = cx;
    _cy = cy;
    _cz = cz;
    _volume = volume;
    _layer = layer;
    _module = module;
    _hitid = hitid;
    _trackid = trackid;

    points.reserve(nhits);
    
    // Set up a cache for the point coordinates
    cout << "Reading hits..." << endl;
    for (int i=0;i<nhits;i++) {
        assignment[i] = 0;
        label[i] = 0;
        treePoint p = treePoint(x[i],y[i],z[i],cx[i],cy[i],cz[i],i,hitid[i],label[i]);
        p.setvolume(volume[i]);
        p.setlayer(layer[i]);
        p.setmodule(module[i]);
        p.settrackid(trackid[i]);
        points.push_back(p);
        hitIDmap[hitid[i]] = i; // Note the new index
    }
    
    // Sort the hits into the detector layers
    cout << "Setting up..." << endl;
    readDetectors(FILEPATH);
    initHitDir();
    readTubes();

    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    // Search neighbouring hits and generate a weighted directed graph
#ifdef PAIRFINDER
    cout << "Searching pairs..." << endl;
    long npairs = findPairs();
#endif

#ifdef SEEDFINDER
    cout << "Searching seeds..." << endl;
    long npairs = findSeeds();
#endif

#ifdef TRUTHFINDER
    cout << "Searching truth..." << endl;
    for (auto p : truepairs) {
        pairs.push_back(p);
        tracking.add(p.first,p.second,1.0);
    }
    long npairs = pairs.size();
#endif

    cout << npairs << " pairs" << endl;
    if (_verbose) {
        for (auto p : pairs) cout << "{" << p.first << "," << p.second << "}, ";
        cout << endl;
    }

    if (SCORE) scorePairs(pairs);

    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    
    // Search triples and add suiting combinations to the graph
    cout << "Searching triples..." << endl;

#ifdef TRIPLEFINDER
    long ntriples = findTriples();
    cout << ntriples << " triples" << endl;
    if (_verbose) {
        for (auto t: triples) cout << t.x << " " << t.y << " " << t.z << "(" << t.r << ") ";
        cout << endl;
    }
#endif

#ifdef TOPQUARKFINDER
   
    loadAdjacent(BASEPATH,"adjacency");
    initDensity3();
    
    PolarModule mod[48];
    for (int i = 0; i < 48; i++)
        mod[i] = PolarModule(i);

    // Transform to trackml indexing
    vector<pair<int,int> > toppairs;
    for (auto pa : pairs) {
        treePoint &p1 = points[pa.first];
        treePoint &p2 = points[pa.second];
        pair<int,int> pb = make_pair(p1.hitid(), p2.hitid());
        toppairs.push_back(pb);
    }
    
    vector<triple> straight_triples = findTriples(toppairs, mod, 1,0.5);
    for (auto&t : straight_triples) {
        triples.push_back(t);
        int id1 = hitIDmap[t.x];
        int id2 = hitIDmap[t.y];
        int id3 = hitIDmap[t.z];
        tracking.add(id2,id3,t.r);
    }
    
#endif

    if (SCORE) scoreTriples(triples);

    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;

    // Assemble tracklets from the seeds/pairs
    map<int,vector<int> > tracklet;
    
#ifdef SWIMMERFINDER
    cout << "Searching tracks with swimmer..." << endl;
    tracklet = swimmer();
#endif
    
#ifdef GRAPHFINDER
    cout << "Searching tracks by analyzing directed weighted graph..." << endl;
    tracklet = serialize(tracking);
#endif
    
    // Exvaluate the results
    cout << endl << "Number of tracklets   : " << tracklet.size() << endl;
    //cout << "Number of short tracks: " << shortpath.size() << endl;
    if (tracklet.size() == 0) exit(0);

    if (SCORE) scorePaths(tracklet);

    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    // Sort the tracklet vector according to the tracklet length
    //sort(tracklet.begin(), tracklet.end(), sortFunc);
    
    // Print out the tracks vector
    if (_verbose) {
        cout << "Tracklets:" << endl;
        for (auto &it : tracklet) print(it.second);
    }
    
    // Assign labels
    int na = 0;
    int nc = 0;
    for (auto it : tracklet) {
        auto v = it.second;
        int id = v[0];
        for (auto j : v) {
            na++;
            label[j] = it.first;
            points[j].setlabel(it.first);
            if (j == id++)
                nc ++;
            else
                id = j++;
        }
    }
    
    // Check the assignment
    
    if (SCORE) {
        map<int, int> map_assignment;
        for (int i = 1; i < hits.size(); i++)
            if (label[i]>0) map_assignment[i] = label[i];
        
        scoreAssignment(map_assignment);
        investigateAssignment(tracklet);
    }

    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;

    cout << "Number of hits                      : " << nhits << endl;
    cout << "Number of double hits               : " << ntwins << endl;
    cout << "Number of assigned points           : " << na << endl;
    cout << "Number of correctly assigned points : " << nc << endl;
    cout << "Distance <" << DISTANCE << ": " << nd <<endl;
    cout << "PhiNN    <" << DELTANN << ": " << np <<endl;
    cout << "Threshold2 " << THRESHOLD2 << ": " << n4 << " R2 OK" << endl;
    cout << "Threshold3 " << THRESHOLD3 << ": " << n1 << " R3 OK" << endl;
    cout << "Recalls2 : " << n2 << endl;
    cout << "Recalls3 : " << n3 << endl << endl;
    
    double perc = 0.0;
    if (seedstotal>0) perc = 100. * seedsok / (double) seedstotal;
    cout << "findSeeds: " << perc << "%" << endl;
    
    perc = 0.0;
    if (trackletstotal>0) perc = 100. * trackletsok / (double) trackletstotal;
    cout << "findTracks: " << perc << "%" << endl << endl;
    
    return (int) tracklet.size();
}


// Write path data to file
void Tracker::writeGraph(string filename, Graph<int> &g) {
    cout << "Writing graph data to " << filename << endl;
    ofstream path(filename);
    path << g;
    path.close();
}


// Read path data from file
void Tracker::readGraph(string filename, Graph<int> &g) {
    cout << "Reading graph data from " << filename << endl;
    ifstream input;
    input.open(filename);
    if (!input.is_open()) {
        cerr << "couldn't open filename" << endl;
    }
    else {
        input >> g;
        input.close();
    }
}


// Run through the graph and reconstruct the tracks from hits
map<int,vector<int> > Tracker::getTracks(Graph<int> &g) {
    return serialize(g);
}


// operations on triples
bool operator<(const triple&a, const triple&b) {
    return a.r < b.r;
}


bool operator==(const triple&a, const triple&b) {
    return a.x==b.x && a.y==b.y && a.z==b.z;
}


//does hits a and b correspond to the same particle?
int Tracker::samepart(treePoint &a, treePoint &b) {
    long long aa = a.trackid();
    long long bb = b.trackid();
    return aa == bb && aa;
}


// Print an integer hit id track vector
void Tracker::print(vector<int> const &input)
{
    for (unsigned int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
    cout << endl;
}


// Sort an hit id track vector by length
bool Tracker::sortFunc( const vector<int>& p1,
                       const vector<int>& p2 ) {
    return p1.size() > p2.size();
}

/*
// Sort two hits wrt distance
bool Tracker::sortByDistance(const int a,const int b) {
    return Point::sortByDistance(points[a],points[b]);
}
*/

// CHeck whether the points in a vector belong to the same track
long Tracker::checkLabels(std::vector<int> &ip)
{
    if (ip.size()==0) return 0;
    int label = points[ip[0]].label();
    if (label<=0) return 0;
    long n = 0;
    for (auto it:ip) {
        if (points[it].label()==label) n++;
    }
    return n;
}


// Generate an index to address the detector layers 0..47
int Tracker::getLayer(int volume_id, int layer_id) {
    
    const int itopo[LAYERS] = {10,9,8,7,6,5,4,0,1,2,3,11,12,13,14,15,16,17,34,32,30,28,26,24,18,19,20,21,36,38,40,42,44,46,35,33,31,29,27,25,22,23,37,39,41,43,45,47};
    const int metai_list[9][7] = {{0,1,2,3,4,5,6},{7,8,9,10,-1,-1,-1},{11,12,13,14,15,16,17},{18,19,20,21,22,23,-1},{24,25,26,27,-1,-1,-1},{28,29,30,31,32,33,-1},{34,35,36,37,38,39,-1},{40,41,-1,-1,-1,-1,-1},{42,43,44,45,46,47}};
    
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


// read the hit data and store it in a layer/phi cache
void Tracker::readTubes() {
    
    cout << "Reading tubes..." << endl;
    
    long nhits = points.size();
    for (int i = 0; i < nhits; i++) {
        int l = points[i].layer();
        if (l<0 || l>=LAYERS) continue;
        int m = points[i].module();
        if (m<0 || m>=MODULES) continue;
        float p = points[i].phi();
        float t = points[i].theta();
        int phi  = PHIFACTOR*(M_PI+p);
        int the  = THEFACTOR*(M_PI+t);
        if (phi>=PHIDIM) phi = PHIDIM-1;
        if (the>=THEDIM) the = THEDIM-1;
        tube[l][phi][the].push_back(i);
        if (_verbose) cout << "Point " << i << " phi:" << phi << " theta:" << the << " layer:" << l << " module: " << m << endl;
    }
    
    for (int i = 0; i < 48; i++) {
        if (layer[i].type == Disc) {
            for (int j=0;j<PHIDIM;j++)
                for (int k=0;k<THEDIM;k++)
                    sort(tube[i][j][k].begin(), tube[i][j][k].end(), z_cmp);
        }
        else {
            for (int j=0;j<PHIDIM;j++)
                for (int k=0;k<THEDIM;k++)
                    sort(tube[i][j][k].begin(), tube[i][j][k].end(), r_cmp);
        }
    }

    // Filter double hits
    cout << "Filter double hits..." << endl;
    for (int i = 0; i < LAYERS; i++) {
        for (int j = 0; j <PHIDIM; j++) {
            for (int k=0;k<THEDIM;k++) {
                auto pvec = tube[i][j][k];
                if (_verbose) {
                    long size = pvec.size();
                    if (size > 0) {
                        cout << "Tube " << i << " phi " << j << " theta " << k << " size: " << size << endl;
                        for (auto it : pvec) cout << it << ",";
                        cout << endl;
                    }
                }
                
                if (pvec.size()<2) continue;
                for (auto it1 = pvec.begin(); it1 != pvec.end()-1; it1++) {
                    int id1 = *it1;
                    treePoint &p1 = points[id1];
                    for (auto it2 = it1+1; it2 != pvec.end(); it2++) {
                        int id2 = *it2;
                        treePoint &p2 = points[id2];
                        double d = distance(id1,id2);
                        //if (_verbose) cout << "Distance " << id0 << "," << id1 << ":" << d << endl;
                        if (d<TWINDIST && checkModule(id1,id2)) { // short distance in different modules
                            if (id1<id2) {
                                p1.settwin(id2);
                                //assignment[id2] = -1;
                            }
                            else {
                                p2.settwin(id1);
                                //assignment[id1] = -11;
                            }
                            if (_verbose) cout << "Twin " << id1 << "," << id2 << ":" << d << endl;
                            ntwins++;
                        }
                    }
                }
            }
        }
    }
    
    // Prepare rhe hits in modules

    for (int i = 0; i < nhits; i++) {
        //if (assignment[i] > 0) continue; // Skip double hits
        int l = points[i].layer();
        if (l<0 || l>=LAYERS) continue;
        int m = points[i].module();
        if (m<0 || m>=MODULES) continue;
        int index = MODULES*l + m;
        module[index].push_back(i);
        modules[l].insert(index);
        //if (_verbose) cout << i << " " << l << " " << m << " " << index << endl;
    }
    
    for (int l =0;l<LAYERS;l++) {
        for (auto m : modules[l]) {
            if (module[m].size()==0) continue;
            if (layer[l].type == Disc)
                sort(module[m].begin(), module[m].end(), z_cmp);
            else
                sort(module[m].begin(), module[m].end(), r_cmp);
            
            if (_verbose) {
                //cout << "Module " << m << ":";
                //print(module[m]);
            }
        }
    }
    

}

bool Tracker::z_cmp(const int a, const int&b) { return points[a].z() < points[b].z(); }
bool Tracker::r_cmp(const int&a, const int&b) { return points[a].rz() < points[b].rz(); }
