#ifndef _TRACKER_H_
#define _TRACKER_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

//#define PAIRS
//#define SWIMMER

#define NETFILE2 "/Users/marcel/workspace/rhonno/trackml/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/trackml/XMLP3.net"

#define TRACKLET 2
#define TWINDIST 0.0051
#define MAXKNN 5
#define THRESHOLD2 0.90
#define THRESHOLD3 0.95
#define DISTANCE 1.0
#define DELTAR   0.5
#define DELTATHE 0.1
#define DELTANN  0.1

#define MAXDIM 120000
#define PHIDIM 13
#define PHIFACTOR 2
#define THEDIM 25
#define THEFACTOR 4

#define SCORE true

#include "Graph.h"
#include "Point.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

inline double dist(double x, double y) { return sqrt(x*x+y*y); }
inline double dist2(double x, double y) { return (x*x+y*y); }

//Basics for 3d coordinate representation
struct point {
    double x, y, z;
    point() {}
    point(double x, double y, double z) : x(x),y(y),z(z) {}
    inline point operator-(const point&p) {
        return point(x-p.x, y-p.y, z-p.z);
    }
    inline point operator+(const point&p) {
        return point(x+p.x, y+p.y, z+p.z);
    }
    inline double operator*(const point&p) {
        return x*p.x+y*p.y+z*p.z;
    }
    inline point operator*(double f) {
        return point(x*f, y*f, z*f);
    }
    inline
    double distance(const point &a)
    {
        return  sqrt((a.x - x) * (a.x - x) +
                     (a.y - y) * (a.y - y) +
                     (a.z - z) * (a.z - z));
    }
};

inline double dist(const point&p) { return sqrt(p.x*p.x+p.y*p.y+p.z*p.z); }

//Structure for storing promising triples of hits
struct triple {
    int x, y, z;   // hit ids
    double r;      // quallity
    triple() {}
    triple(int a, int b, int c, double v=0) : x(a), y(b), z(c), r(v) {}
};

//Geometry of layer
struct Layer {
    double minr, avgr, maxr;
    double minz, avgz, maxz;
    int count;
    int type;    
    double var0, var1;
};

struct Detector {
    int volume_id, layer_id, module_id;
    point c;
    point rx, ry, rz;
    double d, minw, maxw, h, cell_w, cell_h;
};

struct Particle // structure for truth particle info
{
    long long id;
    int type;
    double x;
    double y;
    double z;
    double r;
    double px;
    double py;
    double pz;
    double q;
    long hits;
    std::vector<int> hit;
};

class Tracker {
public:
    static int assignment[MAXDIM]; // hit hs been used
    static Graph<int> paths, tracking;
    static std::vector<int> tube[48][PHIDIM][THEDIM]; // List of hits in each layer
    static std::vector<treePoint> points; // hit Points
    static std::vector<point> hits; //hit position
    static std::vector<point> polar; //hit position in polar / cylindrical coordinates
    static std::vector<Particle> particles; //true tracks
    static std::map<long long,int> partIDmap; // create particle ID->index map
    static std::map<long long, std::vector<int> > truth_tracks; //truth hit ids in each track
    static std::map<long long, point> track_hits; // Find points in hits
    static std::vector<int> metai, metaz; //ordered layer id in [0,48), and classification of z for disc layers in [0,4)
    static std::vector<point> meta; //volume_id / layer_id / module_id
private:
    static std::vector<int> knn[MAXDIM][PHIDIM][THEDIM];
    static point truth_pos[MAXDIM], truth_mom[MAXDIM]; //truth position and momentum
    static double truth_weight[MAXDIM]; //weighting of each hit
    static long long truth_part[MAXDIM]; //particle this hit belongs to
    static std::set<long long> blacklist;
    static std::map<long long, double> part_weight; //weighting of each particle
    static std::map<long long, std::map<int, double> > metai_weight; //weighting of each particle hit, also adding duplicates
    static std::map<long long, point> start_pos; //start position
    static std::map<long long, point> start_mom; //start momentum
    static std::map<long long, int> part_q; //start charge
    static std::map<long long, int> part_hits; // = truth_tracks[particle_id].size()
    static int topo[48], itopo[48]; //reordering of layers for approximate sorting
    static double disc_z[48][4];
    static const int Tube = 0, Disc = 1;
    static constexpr double Bfield = 1673.0; //Empirical field strengh, to scale the momentum
    static Layer layer[48];
    static double z_minr[48][4], z_maxr[48][4];
    static std::map<int, Detector> detectors;
    static std::vector<std::pair<std::pair<int, int>, double> > hit_cells[MAXDIM]; //pair<pair<ch0, ch1>, value>
    static point hit_dir[MAXDIM][2]; //The two possible directions of the hit according to the cell's data for each hit

    static unsigned long nr, nd, np, nt, n1, n2, n3, n4, ntwins;
    static bool _verbose;
    static float *_x;
    static float *_y;
    static float *_z;
public:
    Tracker() {}
    inline
    static float distance(const int &a, const int &b)
    {
        return  sqrt((_x[a] - _x[b]) * (_x[a] - _x[b]) +
                     (_y[a] - _y[b]) * (_y[a] - _y[b]) +
                     (_z[a] - _z[b]) * (_z[a] - _z[b]));
    }
    inline
    static float radius(const int &a) { return sqrt(_x[a]*_x[a]+_y[a]*_y[a]+_z[a]*_z[a]); }
    static void verbose(bool verbose=true) {_verbose = verbose;}
    static int findTracks(int nhits,float *x,float *y,float *z,int *layer,int *label,int *truth);
    static void setupCache();
    static std::map<int,std::vector<int> > swimmer();
    static std::map<int,std::vector<int> >  getTracks(Graph<int> &g);
    static std::vector<std::pair<int,float> > findSeeds(int p,std::vector<int> &points);
    static std::vector<std::pair<int, int> > findSeeds();
    static std::vector<std::pair<int, int> > findPairs();
    static long findTriples(std::vector<std::pair<int,int> > seed, std::vector<triple> &triples);
    static long findTriples(int p0,int p1,std::vector<int> &points,std::vector<triple> &triples);
    static long addHits(int p0, int p1, int layer,int phi,int the,std::vector<triple> &triples);
    static long addHitsCached(int p0,int p1,int phi,std::vector<triple> &triples);
    inline
    static bool checkRadius(const int p0,const int p1) { double dr = abs(points[p0].r()-points[p0].r()); if (dr > DELTAR) { nr++; return false;} else return true; }
    inline
    static bool checkTheta(const int p0,const int p1) { float dt = fabs(points[p0].theta()-points[p1].theta()); if (dt > DELTATHE) { nt++; return false; } else return true; }
    inline
    static bool checkDistance(const int p0,const int p1) {
        float d = distance(p0,p1);
        if (d > DISTANCE*radius(p0)) {
            nd++;
            return false;
        }
        else return true;
    }
    static double checkTracklet(int p0,int p1);
    static double checkTracklet(int p0,int p1,int p2);
    static long checkLabels(std::vector<int> &p);
    static long checkTracks(std::map<int,std::vector<int> >  &tracks);
    static long seedstotal,seedsok;
    static long trackletstotal,trackletsok;
    static void readBlacklist(std::string base_path,int filenum);
    static void readTruth(std::string base_path,int filenum);
    static void readParticles(std::string base_path,int filenum);
    static void readHits(std::string base_path,int filenum);
    static void readCells(std::string base_path,int filenum);
    static void readDetectors(std::string base_path);
    static void readTubes();
    static void sortTracks();
    static int getLayer(int volume, int layer);
    static void scorePairs(std::vector<std::pair<int, int> >&pairs);
private:
    static void print(std::vector<int> const &input);
    static bool sortFunc( const std::vector<int>& p1,const std::vector<int>& p2 );
    static bool sortDist(const int a,const int b);
    static double* recall2(Point &p1, Point &p2);
    static double* recall3(Point &p1, Point &p2, Point &p3);
    static bool z_cmp(const int a, const int&b);
    static bool r_cmp(const int&a, const int&b);
    static int samepart(int a, int b);
    static bool track_cmp(int a, int b);
    static void initOrder();
    static void initLayers();
    static point normalize(point a);
    static void initHitDir();
};

#endif
