// Read the trackml data files and extract neural network training data
// M.Kunze, Heidelberg University, 2018

#include <iostream>
#include <set>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <stack>
#include <queue>

#define MAXPARTICLES 10
#define MAXHITS 150000

//Not doing much in practice
int debug = 0;
//Not the standard assert!
#define assert(a, m) {if (!(a)) {cout << m << endl; exit(0);}}

using namespace std;

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
};
ostream& operator<<(ostream&out, const point p) {
    out << p.x << ' ' << p.y << ' ' << p.z;
    return out;
}
inline double dist(const point&p) { return sqrt(p.x*p.x+p.y*p.y+p.z*p.z); }

//Structure for storing promising triples of hits
/*
struct triple {
    int x, y, z; //hit ids
    triple() {}
    triple(int a, int b, int c) : x(a), y(b), z(c) {}
};
bool operator<(const triple&a, const triple&b) {
    if (a.z != b.z) return a.z < b.z;
    if (a.y != b.y) return a.y < b.y;
    return a.x < b.x;
}
bool operator==(const triple&a, const triple&b) {
    return a.x==b.x && a.y==b.y && a.z==b.z;
}
*/
//Convert "dir" to polar (really cylindrical coordinates) coordinates, suffix p  (as in "refp") usually means polar coodinates throughout the code
point topolar(const point&dir, const point&ref, const point&refp) {
    return point(ref.x*dir.x+ref.y*dir.y, ref.x*dir.y-ref.y*dir.x, dir.z*refp.x);
}

#include "input.h"

//does hits a and b correspond to the same particle?
int samepart(int a, int b) {
    long long aa = truth_part[a];
    long long bb = truth_part[b];
    return aa == bb && aa;
}

#include "Tracker.h"

int main(int argc, char**argv) {
    //Read which event to run from arguments, default is event # 1000
    //Supports events 0-124 for test, or 1000-1099 for validation (small dataset)

    if (argc >= 2) {
        filenum = atoi(argv[1]);
        cout << "Running on event #" << filenum << endl;
    }
    ios::sync_with_stdio(false);
    cout << fixed;

    int eval = 0;
#if defined EVAL
    eval = 1;
#endif
    if (!eval) {
        readBlacklist();
        readTruth();
        sortTracks();
        readParticles();
    }
    readHits();

    long nParticles = truth_tracks.size();
    if (nParticles > MAXPARTICLES) nParticles = MAXPARTICLES;
    cout << "Particles: " << nParticles << endl;

    
    long nhits = hits.size();
    float x[nhits],y[nhits],z[nhits],weights[nhits];
    int labels[nhits],volumes[nhits],layers[nhits],modules[nhits];
    
    int i = 0;
    int n = 0;
    for (auto &track : truth_tracks) {
        if (n++ > MAXPARTICLES) break;
        vector<int> t = track.second;
        for (auto &id : t) {
            auto it = track_hits.find(id);
            if (it==track_hits.end()) continue;
            point &hit = it->second;
            x[i] = hit.x * 0.001; // in m
            y[i] = hit.y * 0.001; // in m;
            z[i] = hit.z * 0.001; // in m;
            i++;
        }
    }

    nhits = i-1;
    if (nhits > MAXHITS) nhits = MAXHITS;
    cout << "Hits: " << nhits << endl;
    
    cout << "Find tracks..." << endl;
    Tracker::findTracks((int)nhits,x,y,z,labels,volumes,modules,layers,weights);

}
