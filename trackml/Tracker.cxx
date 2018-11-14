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


// Find tracks from points
int Tracker::findTracks(int nhits,float *x,float *y,float *z,int* layer,int* module,int* label,int *truth)
{
    std::clock_t c_start = std::clock();
    
    _x = x;
    _y = y;
    _z = z;
    
    points.reserve(nhits);
    
    // Set up a cache for the point coordinates
    cout << "Reading hits..." << endl;
    for (int i=0;i<nhits;i++) {
        assignment[i] = 0;
        label[i] = 0;
        treePoint p = treePoint(x[i],y[i],z[i],i,label[i],truth[i]);
        p.setlayer(layer[i]);
        p.setmodule(module[i]);
        points.push_back(p);
    }
    
    // Sort the hits into the detector layers
    cout << "Setting up..." << endl;
    readTubes();
    //setupCache();
    
    std::clock_t c_end = std::clock();
    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;
    
    // Search neighbouring hits and generate a weighted directed graph
#ifdef PAIRS
    cout << "Searching pairs..." << endl;
    long npairs = findPairs();
#endif

#ifdef SEEDS
    cout << "Searching seeds..." << endl;
    long npairs = findSeeds();
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

#ifdef TRIPLETS
    triples.clear();
    long ntriples = findTriples();
    cout << ntriples << " triples" << endl;
    if (_verbose) {
        for (auto t: triples) cout << t.x << " " << t.y << " " << t.z << "(" << t.r << ") ";
        cout << endl;
    }
    
    if (SCORE) scoreTriples(triples);
    
    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n" <<endl;

#endif
    
    // Assemble tracklets from the seeds/pairs
    map<int,vector<int> > tracklet;
    
#ifdef SWIMMER
    cout << "Searching tracks with swimmer..." << endl;
    tracklet = swimmer();
#endif
    
#ifdef GRAPH
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
int Tracker::samepart(int a, int b) {
    long long aa = points[a].truth();
    long long bb = points[b].truth();
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


// Sort two hits wrt distance
bool Tracker::sortDist(const int a,const int b) {
    return Point::sortDist(points[a],points[b]);
}


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


// Recall function for 2 points
double Tracker::checkTracklet(int p0,int p1)
{
    double recall = recall2(points[p0],points[p1])[0];
    bool ok = recall>THRESHOLD2;
    recall = ok ? recall : -recall;
    if (ok) n4++;
    n2++;
    return recall;
}


// Recall function for triple
double Tracker::checkTracklet(int p1,int p2,int p3)
{
    Point &point1 = points[p1];
    Point &point2 = points[p2];
    Point &point3 = points[p3];

    Point v1 = point2 - point1;
    Point v2 = point3 - point1;
    double angle = acos(Point::dot(v1,v2)); // Check angle between the last points
    if (angle<0.5*M_PI) {
        if (angle>DELTANN) { np++; return 0.0; } // Check 0 deg.
    }
    else {
        if ((M_PI-angle)>DELTANN) { np++; return 0.0; } // Check 180 deg.
    }
    
    double recall = recall3(point1,point2,point3)[0];
    bool ok = recall>THRESHOLD3;
    recall = ok ? recall : -recall;
    if (ok) n1++;
    n3++;
    return recall;
}


// Recall function for 2 points
double* Tracker::recall2(Point &p1, Point &p2)
{
    static XMLP net(NETFILE2);
    static float x[6]={0.,0.,0.,0.,0.,0.};
    
    x[0]     = p1.rz();     // rz1
    x[1]     = p1.phi();    // phi1
    x[2]     = p1.z();      // z1
    x[3]     = p2.rz();     // rz2
    x[4]     = p2.phi();    // phi2
    x[5]     = p2.z();      // z2
    
    return net.Recallstep(x);
}


// Recall function for 3 points
double* Tracker::recall3(Point &p1, Point &p2, Point &p3)
{
    static XMLP net(NETFILE3);
    static float x[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    
    x[0]     = p1.rz();     // rz1
    x[1]     = p1.phi();    // phi1
    x[2]     = p1.z();      // z1
    x[3]     = p2.rz();     // rz2
    x[4]     = p2.phi();    // phi2
    x[5]     = p2.z();      // z2
    x[6]     = p3.rz();     // rz3
    x[7]     = p3.phi();    // phi3
    x[8]     = p3.z();      // z3
    
    return net.Recallstep(x);
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
        //if (assignment[i] > 0) continue; // Skip doiuble hits
        int l = points[i].layer();
        if (l<0 || l>=LAYERS) continue;
        int m = points[i].module();
        if (m<0 || m>=MODULES) continue;
        int index = MODULES*l + m;
        module[index].push_back(i);
        modules[l].insert(index);
    }
    
    for (int i = 0; i < MODULES; i++) {
        if (layer[i].type == Disc)
            sort(module[i].begin(), module[i].end(), z_cmp);
        else
            sort(module[i].begin(), module[i].end(), r_cmp);
    }
    

}

bool Tracker::z_cmp(const int a, const int&b) { return points[a].z() < points[b].z(); }
bool Tracker::r_cmp(const int&a, const int&b) { return points[a].rz() < points[b].rz(); }


// functions to read trackml data (taken from topquark)

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
        if( newID < 0 || newID>= (int)particles.size() ){
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
        for (unsigned int i = 0; i < v.size(); i++)
            p.second[i] = v[i];
        int bad = 0;
        for (unsigned int i = 2; i < v.size(); i++) {
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
    int handpicked[LAYERS] = {};
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
    for (int i = 0; i < LAYERS; i++) {
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
    Layer layer2[LAYERS];
    for (int i = 0; i < LAYERS; i++) layer2[i] = layer[i];
    for (int i = 0; i < LAYERS; i++) layer[i] = layer2[topo[i]];
    
    layer[0].var0 = 1e-3;
    layer[1].var0 = 5e-4;
    for (int i = 2; i < 18; i++) layer[i].var0 = 3e-4;
    for (int i = 18; i < 22; i++) layer[i].var0 = 5e-2;
    for (int i = 22; i < LAYERS; i++) layer[i].var0 = i%2 || i == 22 ? 9 : 0.1;
    
    for (int i = 0; i < 4; i++) layer[i].var1 = 0.5;
    for (int i = 4; i < 18; i++) layer[i].var1 = 5;
    for (int i = 18; i < 24; i++) layer[i].var1 = 7;
    for (int i = 24; i < LAYERS; i++) layer[i].var1 = i%2 ? 19 : 11;
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
    
    for (int i = 0; i < LAYERS; i++) {
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
    
    for (unsigned int hit_id = 1; hit_id < hits.size(); hit_id++) {
        metai_weight[truth_part[hit_id]][metai[hit_id]] += truth_weight[hit_id];
    }
    
    map<double, double> mir[LAYERS], mar[LAYERS];
    for (unsigned int i = 1; i < hits.size(); i++) {
        int mi = metai[i];
        if (layer[mi].type != Disc) continue;
        double &mir_ = mir[mi][polar[i].z];
        double &mar_ = mar[mi][polar[i].z];
        if (!mir_) mir_ = 1e9;
        mir_ = min(mir_, polar[i].x);
        mar_ = max(mar_, polar[i].x);
    }
    
    map<double, int> zi[LAYERS];
    for (int mi = 0; mi < LAYERS; mi++) {
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
    for (unsigned int i = 1; i < hits.size(); i++) {
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
    for (unsigned int hit_id = 1; hit_id < hits.size(); hit_id++) {
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


void Tracker::scorePairs(vector<pair<int, int> >&pairs) {
    set<long long> found;
    for (auto p : pairs) {
        if (samepart(p.first, p.second)) {
            found.insert(truth_part[p.first]);
        }
    }
    double score = 0;
    for (long long p : found) {
        score += part_weight[p];
    }
    cout << "Score: " << score << " from " << pairs.size() << " pairs" << endl;
}


void Tracker::scoreTriples(vector<triple>&triples) {
    set<long long> found;
    for (auto p : triples) {
        if (samepart(p.x, p.y) && samepart(p.x, p.z)) {
            found.insert(truth_part[p.x]);
        }
    }
    double score = 0;
    for (long long p : found) {
        score += part_weight[p];
    }
    cout << "Score: " << score << " from " << triples.size() << " triples" << endl;
}


void Tracker::scorePaths(map<int,vector<int> > &paths) {
    int total_length = 0;
    map<long long, double> score_part;
    for (auto &path : paths) {
        
        total_length += path.second.size();
        map<long long, int> count;
        for (int i : path.second)
            count[truth_part[i]]++;
        long long part=0;
        int ma = 0;
        for (auto p : count) {
            if (p.second >= ma) {
                ma = p.second;
                part = p.first;
            }
        }
        double w = 0;
        for (int j : path.second) {
            if (truth_part[j] == part) w += truth_weight[j];
        }
        score_part[part] = max(score_part[part], w);
    }
    double score = 0;
    for (auto p : score_part) if (p.first) score += p.second;
    cout << "Score: " << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
}


void Tracker::scoreAssignment(map<int, int>& assignment) {
    map<int, int> track_length;
    for (auto p : assignment)
        track_length[p.second]++;
    
    double score = 0, score2 = 0;
    for (auto it : truth_tracks) {
        map<int, int> c;
        for (int i : it.second)
            if (assignment.count(i))
                c[assignment[i]]++;
        int pick = -1, maxlen = -1;
        for (auto p : c)
            if (p.second*2 > maxlen &&
                p.second*2 > track_length[p.first]) { pick = p.first; maxlen = p.second*2; }
        if (pick == -1) continue;
        
        for (int i : it.second)
            if (assignment[i] == pick) {
                if (maxlen > it.second.size()) score += truth_weight[i];
                else score2 += truth_weight[i];
            }
    }
    //"Final score: " should be the same as the official score (except for blacklisted electrons)
    cout << "Final score:  " << score << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
    cout << "Short score:  " << score+score2 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
}

//Investigate some stats on the assignment, to try to see of there are any more points to gather
void Tracker::investigateAssignment(map<int,vector<int> > &solution_paths) {
    cout << solution_paths.size() << endl;
    int good = 0;
    map<long long, int> goodc;
    double stolen = 0, missing = 0, between = 0, duplicate = 0, wrong = 0;
    for (int i = 1; i < solution_paths.size(); i++) {
        map<long long, int> c;
        long long bestp = -1;
        int mi = solution_paths[i].size();
        int missed = 0;
        for (int j : solution_paths[i]) {
            if (j <= 0) {missed++;continue;}
            c[truth_part[j]]++;
            if (c[truth_part[j]]*2 > mi) bestp = truth_part[j];
        }
        if (bestp != -1) {
            goodc[bestp]++;
            if (goodc[bestp] == 1) {
                set<int> found, m;
                for (int j : truth_tracks[bestp]) found.insert(j);
                int mi = 1e9, ma = -1;
                for (int j : solution_paths[i]) {
                    int k = abs(j);
                    if (truth_part[k] == bestp && j <= 0) stolen += truth_weight[k];
                    if (truth_part[k] == bestp) {
                        found.erase(k);
                        mi = min(mi, metai[k]);
                        ma = max(ma, metai[k]);
                        m.insert(metai[k]);
                    }
                }
                for (int j : found) {
                    missing += truth_weight[j];
                    if (metai[j] > mi && metai[j] < ma) between += truth_weight[j];
                    if (m.count(metai[j])) duplicate += truth_weight[j];
                }
                good++;
            }
            //log(-scorepathDensity(solution_paths[i])) << endl;
        } else {
            int bad = 0;
            for (int k = 1; k < solution_paths[i].size()-1; k++) {
                auto&v = solution_paths[i];
                if (v[k-1] > 0 && v[k] < 0 && v[k+1] > 0) bad = 1;
            }
            //cerr << bad*1./solution_paths[i].size() << endl;
            for (int j : solution_paths[i])
                if (j > 0)
                    wrong += truth_weight[j];
        }
    }
    cout << good << endl;
    cout << "Stolen: " << stolen << endl;
    cout << "Missing: " << missing << endl;
    cout << "Between: " << between << endl;
    cout << "Duplicate: " << duplicate << endl;
    cout << "Wrong: " << wrong << endl;
    double outside = 0;
    for (auto&p : part_weight)
        if (!goodc.count(p.first)) outside += p.second;
    cout << "Outside: " << outside << endl;
    
    
    map<int, int> map_assignment;
    for (int i = 1; i < solution_paths.size(); i++)
        for (int j : solution_paths[i])
            if (j > 0)
                map_assignment[j] = i;
    
    double completely = 0;
    for (auto&t : truth_tracks) {
        map<int, int> c;
        for (int i : t.second)
            if (map_assignment[i])
                c[map_assignment[i]]++;
        int ma = 0;
        for (auto&p : c) ma = max(ma, p.second);
        if (ma <= 2) completely += part_weight[t.first];
    }
    cout << "Completely lost: " << completely << endl;
}


//Score triple based on the deviation from a perfect helix, no prior that it should be straight
double Tracker::scoreTriple(int ai, int bi, int ci) {
    point center;
    double radius;
    circle(hits[ai], hits[bi], hits[ci], center, radius);
    
    point cb = hits[ci]-hits[bi];
    point ba = hits[bi]-hits[ai];
    double ang_cb = asin(dist(cb.x, cb.y)*.5/radius)*2;
    double ang_ba = asin(dist(ba.x, ba.y)*.5/radius)*2;
    if (radius != radius || fabs(radius) > 1e50) {
        ang_cb = dist(cb.x, cb.y);
        ang_ba = dist(ba.x, ba.y);
    }
    if (ba.z*cb.z < 0) ang_ba *= -1;
    
    //if (dist(cb.x, cb.y)*.5/radius > M_PI/2 || dist(ba.x, ba.y)*.5/radius > M_PI/2) return 1e9;
    //-radius*2e-5+
    double x = ba.z ? (fabs(cb.z*ang_ba/ba.z-ang_cb))*radius : 1e9;
    double y = ang_cb ? (fabs(cb.z*ang_ba/ang_cb-ba.z)) : 1e9;
    double score = min(x, y);//, fabs(cb.z-ba.z*ang_cb/ang_ba)));
    /*
     cout << endl;
     cout << truth_mom[bi]*part_q[truth_part[bi]] << endl;
     point rr = hits[bi]-center;
     if (cb.x*rr.y-cb.y*rr.x > 0) ang_cb *= -1;
     cout << point(-rr.y, rr.x, cb.z/ang_cb) << endl;*/
    return score;
}

// Data initialization

float* Tracker::_x;
float* Tracker::_y;
float* Tracker::_z;
bool Tracker::_verbose(false);
vector<treePoint> Tracker::points; // hit Points
vector<pair<int, int> > Tracker::pairs; // hit pair combinations
vector<triple> Tracker::triples; // hit triple combinations
Graph<int> Tracker::paths, Tracker::tracking; // graph to represent particle paths and tracking information
long Tracker::seedsok(0),Tracker::seedstotal(0);
long Tracker::trackletsok(0),Tracker::trackletstotal(0);
unsigned long Tracker::nd(0),Tracker::np(0),Tracker::nt(0);
unsigned long Tracker::n1(0),Tracker::n2(0),Tracker::n3(0),Tracker::n4(0),Tracker::ntwins(0);
vector<point> Tracker::hits; //hit position
vector<Particle> Tracker::particles; //true tracks
map<long long,int> Tracker::partIDmap; // create particle ID->index map
//vector<int> Tracker::knn[MAXDIM][MODULES];
vector<int> Tracker::module[LAYERS*MODULES]; // List of hits in each layer
set<int> Tracker::modules[LAYERS]; // List of modules in each layer
vector<int> Tracker::tube[LAYERS][PHIDIM][THEDIM]; // List of hits in each layer
map<long long, vector<int> > Tracker::truth_tracks; //truth hit ids in each track
map<long long, point> Tracker::track_hits; // Find points in hits
int Tracker::assignment[MAXDIM];
point Tracker::truth_pos[MAXDIM], Tracker::truth_mom[MAXDIM]; //truth position and momentum
double Tracker::truth_weight[MAXDIM]; //weighting of each hit
long long Tracker::truth_part[MAXDIM]; //particle this hit belongs to
set<long long> Tracker::blacklist;
map<long long, double> Tracker::part_weight; //weighting of each particle
map<long long, map<int, double> > Tracker::metai_weight; //weighting of each particle hit, also adding duplicates
map<long long, point> Tracker::start_pos; //start position
map<long long, point> Tracker::start_mom; //start momentum
map<long long, int> Tracker::part_q; //start charge
map<long long, int> Tracker::part_hits; // = truth_tracks[particle_id].size()
int Tracker::topo[LAYERS], Tracker::itopo[LAYERS]; //reordering of layers for approximate sorting
vector<point> Tracker::polar; //hit position in polar / cylindrical coordinates
vector<point> Tracker::meta; //volume_id / layer_id / module_id
vector<int> Tracker::metai, Tracker::metaz; //ordered layer id in [0,48), and classification of z for disc layers in [0,4)
double Tracker::disc_z[LAYERS][4];
Layer Tracker::layer[LAYERS];
double Tracker::z_minr[LAYERS][4], Tracker::z_maxr[LAYERS][4];
map<int, Detector> Tracker::detectors;
vector<std::pair<pair<int, int>, double> > Tracker::hit_cells[MAXDIM]; //pair<pair<ch0, ch1>, value>
point Tracker::hit_dir[MAXDIM][2]; //The two possible directions of the hit according to the cell's data for each hit

