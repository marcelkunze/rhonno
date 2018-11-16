#ifndef _TRACKER_H_
#define _TRACKER_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#define SEEDS
//#define PAIRS
#define TRIPLETS
//#define SWIMMER
#define GRAPH

#define FILEPATH "/Users/marcel/workspace/train_sample/"
#define NETFILE1 "/Users/marcel/workspace/rhonno/trackml/XMLP1.net"
#define NETFILE2 "/Users/marcel/workspace/rhonno/trackml/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/trackml/XMLP3.net"

#define TRACKLET 2
#define TWINDIST 5.1
#define THRESHOLD2 0.50
#define THRESHOLD3 0.90
#define DISTANCE 180.0
#define DELTANN  0.08

#define MAXDIM 150000
#define LAYERS 48
#define MODULES 10000
#define PHIDIM 7
#define PHIFACTOR 1
#define THEDIM 7
#define THEFACTOR 1

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
    inline point operator-(const point &p) { return point(x-p.x, y-p.y, z-p.z); }
    inline point operator+(const point &p) { return point(x+p.x, y+p.y, z+p.z); }
    inline double operator*(const point &p) { return x*p.x+y*p.y+z*p.z; }
    inline point operator*(double f) { return point(x*f, y*f, z*f); }
    inline point cross(point f) const { return point(y*f.z-z*f.y, z*f.x-x*f.z, x*f.y-y*f.x); }
    inline
    double distance(const point &a)
    {
        return  sqrt((a.x - x) * (a.x - x) +
                     (a.y - y) * (a.y - y) +
                     (a.z - z) * (a.z - z));
    }

};

inline double dist(const point &p) { return sqrt(p.x*p.x+p.y*p.y+p.z*p.z); }

inline double dist3(point &a,point &b,point &c) {
    const point x = a-b;
    const point y = a-c;
    const point z = c-b;
    double d = dist( x.cross(y)) / dist(z);
    return d;
    
}

inline point normalize(point a) {
    point ret = a*(1./dist(a));
    if (ret.z < 0) ret = ret*-1;
    return ret;
}

//Find circle with center p, radius r, going through a, b, and c (in xy plane)
inline void circle(point&a, point&b, point&c, point&p, double &r) {
    double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
    double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
    double idet = .5/(ax*by-ay*bx);
    p.x = (aa*by-bb*ay)*idet;
    p.y = (ax*bb-bx*aa)*idet;
    p.z = 0;
    r = dist(p.x, p.y);
    p.x += c.x;
    p.y += c.y;
}

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
    static std::vector<std::pair<int, int> > pairs; // hit pair combinations
    static std::vector<triple> triples; // hit triple combinations
    static Graph<int> paths, tracking; // graph to represent particle paths and tracking information
    static std::vector<int> module[LAYERS*MODULES]; // List of hits in each module
    static std::set<int> modules[LAYERS]; // List of modules in each layer
    static std::vector<int> tube[LAYERS][PHIDIM][THEDIM]; // List of hits in each layer
    static std::vector<treePoint> points; // hit Points
    static std::vector<point> hits; //hit position
    static std::vector<point> polar; //hit position in polar / cylindrical coordinates
    static std::vector<Particle> particles; //true tracks
    static std::map<long long,int> partIDmap; // create particle ID->index map
    static std::map<long long, std::vector<int> > truth_tracks; //truth hit ids in each track
    static std::map<long long, point> track_hits; // Find points in hits
    static std::vector<int> metai, metaz; //ordered layer id in [0,48), and classification of z for disc layers in [0,4)
    static std::vector<point> meta; //volume_id / layer_id / module_id
    static std::vector<std::pair<std::pair<int, int>, double> > hit_cells[MAXDIM]; //pair<pair<ch0, ch1>, value>
    static point hit_dir[MAXDIM][2]; //The two possible directions of the hit according to the cell's data for each hit
private:
    //static std::vector<int> knn[MAXDIM][MODULES];
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
    static int topo[LAYERS], itopo[LAYERS]; //reordering of layers for approximate sorting
    static double disc_z[LAYERS][4];
    static const int Tube = 0, Disc = 1;
    static constexpr double Bfield = 1673.0; //Empirical field strengh, to scale the momentum
    static Layer layer[LAYERS];
    static double z_minr[LAYERS][4], z_maxr[LAYERS][4];
    static std::map<int, Detector> detectors;

    static unsigned long nd, np, nt, n1, n2, n3, n4, ntwins;
    static bool _verbose;
    static int *_hitid;
    static int *_volume;
    static int *_layer;
    static int *_module;
    static float *_cx;
    static float *_cy;
    static float *_cz;
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
    static int findTracks(int nhits,float *x,float *y,float *z,float *cx,float *cy,float *cz,int* volume,int *layer,int *module,int *label,int *truth,int *hitid);
    static std::map<int,std::vector<int> > swimmer();
    static std::map<int,std::vector<int> >  getTracks(Graph<int> &g);
    static long findSeedsPhiTheta();
    static long findSeeds();
    static long findPairs();
    static long findTriples();
    static long findTriples(int p0,int p1,std::vector<int> &points);
    static long addHits(int p0, int p1, int module,std::vector<triple> &triples);
    inline
    static bool checkModule(const int p0,const int p1) {return points[p0].module()!=points[p1].module();}
    inline
    static bool checkTheta(const int p0,const int p1) {
        int x1 = THEFACTOR*(M_PI+points[p0].theta());
        int x2 = THEFACTOR*(M_PI+points[p1].theta());
        return x1 == x2;
    }
    inline
    static bool checkPhi(const int p0,const int p1) {
        int x1 = PHIFACTOR*(M_PI+points[p0].phi());
        int x2 = PHIFACTOR*(M_PI+points[p1].phi());
        return x1 == x2;
    }
    static double checkTracklet(int p0,int p1);
    static double checkTracklet(int p0,int p1,int p2);
    static double* recall2(Point &p1, Point &p2);
    static double* recallPair(Point &p1, Point &p2);
    static double* recall3(Point &p1, Point &p2, Point &p3);
    static long checkLabels(std::vector<int> &p);
    static long checkTracks(std::map<int,std::vector<int> >  &tracks);
    static long seedstotal,seedsok;
    static long trackletstotal,trackletsok;
    static void writeGraph(std::string filename, Graph<int> &g);
    static void readGraph(std::string filename, Graph<int> &g);
    static void readBlacklist(std::string base_path,int filenum);
    static void readTruth(std::string base_path,int filenum);
    static void readParticles(std::string base_path,int filenum);
    static void readHits(std::string base_path,int filenum);
    static void readCells(std::string base_path,int filenum);
    static void readDetectors(std::string base_path);
    static void initHitDir();
    static void readTubes();
    static void sortTracks();
    static int getLayer(int volume, int layer);
    static int good_pair(int a, int b);
    static double dir_miss(int ai, int bi);
    static double scoreTripleLogRadius_and_HitDir(int ai,int bi,int ci,float* L);
    static bool getFeatures3(int ai, int bi, float *feature);
    static double wdistr(double r1, double dr, double az, double dz, double w);
    static double wdist(Point&a, Point&d, double w);
    static double zdist(Point&a, Point&b);
    static double zdist2(Point&a, Point&b);
    static void scorePairs(std::vector<std::pair<int, int> > &pairs);
    static void scoreTriples(std::vector<triple> &triples);
    static void scorePaths(std::map<int,std::vector<int> > &paths);
    static void scoreAssignment(std::map<int, int>& assignment);
    static void investigateAssignment(std::map<int,std::vector<int> > &solution_paths);
    static double scoreTriple(int ai, int bi, int ci);
    static void print(std::vector<int> const &input);
private:
    static bool sortFunc( const std::vector<int>& p1,const std::vector<int>& p2 );
    static bool sortDist(const int a,const int b);
    static bool z_cmp(const int a, const int&b);
    static bool r_cmp(const int&a, const int&b);
    static int samepart(int a, int b);
    static bool track_cmp(int a, int b);
    static void initOrder();
    static void initLayers();
};

#endif
