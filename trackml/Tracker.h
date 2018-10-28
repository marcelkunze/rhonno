#ifndef _TRACKER_H_
#define _TRACKER_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#define TRACKML

#ifdef TRACKML
#define NETFILE2 "/Users/marcel/workspace/rhonno/trackml/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/trackml/XMLP3.net"
#define RMIN  0.0
#define RMAX  0.2
#define ZMIN -0.7
#define ZMAX  0.7
#define TRACKLET 2
#define MAXKNN 150
#define THRESHOLD2 0.75
#define THRESHOLD3 0.95
#define DISTANCE 1.2
#define DELTAR   0.9
#define DELTATHE 0.1
#define DELTAPHI 0.1
#define DELTANN  0.1
#else
#define NETFILE2 "/Users/marcel/workspace/rhonno/RhoNNO/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/RhoNNO/XMLP3.net"
#define RMIN 0.0
#define RMAX 10.0
#define ZMIN 0.0
#define ZMAX 10.0
#define TRACKLET 2
#define MAXKNN 8
#define THRESHOLD 0.9
#define DISTANCE 0.5
#define DELTAR   1.0
#define DELTAPHI 0.04
#define DELTANN  0.02
#define NHITS 10
#define SIGMA 0.0
#endif

#define VERBOSE true

#include <cmath>
#include <algorithm>
#include <vector>

static unsigned long nr, nd, np, nt, nx, n1, n2, n3, n4;

class Point;
struct triple;

class Tracker {
private:
    static void print(std::vector<int> const &input);
    static bool sortFunc( const std::vector<int>& p1,const std::vector<int>& p2 );
    static double* Recall2(Point &p1, Point &p2);
    static double* Recall3(Point &p1, Point &p2, Point &p3);
public:
    Tracker() {}
    static int findTracks(int nhits, float *x, float *y, float *z, int* labels);
    static long findSeeds(Point &p,std::vector<Point> &points,std::vector<Point> &seeds);
    static long findTriples(Point &p,std::vector<Point> &points,std::vector<triple> &triples);
    static long selectPoints(std::vector<Point> &points, std::vector<Point> &inner, std::vector<Point> &outer, double rmin, double rmax, double zmin, double zmax);
    static long selectPoints(std::vector<Point> &points, std::vector<Point> &good, std::vector<Point> &bad, Point &ref, double deltar, double deltathe, double distance);
    static double checkTracklet(Point &p0,Point &p1);
    static double checkTracklet(Point &p0,Point &p1, Point &p2);
    static long checkLabels(std::vector<Point> &p);
    static long seedstotal,seedsok;
    static long trackletstotal,trackletsok;
    static long nr,nd,np,nt,nx,n1,n2,n3,n4;

};

#endif
