#ifndef _TRACKER_H_
#define _TRACKER_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#define TRACKML

#ifdef TRACKML
#define NETFILE2 "/Users/marcel/workspace/rhonno/trackml/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/trackml/XMLP3.net"
#define RMIN 0.0
#define RMAX 0.2
#define ZMIN 0.0
#define ZMAX 0.7
#define TRACKLET 2
#define NEIGHBOURS 3
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
#define NEIGHBOURS 3
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

class Point
{
private:
    int _id;                // Hit id
    int _label;             // Label
    int _truth;             // True label
    double _x, _y, _z;      // Cartesian coordinate
    double _r,_phi,_theta;  // Spherical coordinates
    double _rz;             // Cylindrical coordinates
    double _distance;       // Distance from test point
    int _neighbour[NEIGHBOURS];
    double _recall[NEIGHBOURS];
public:
    Point(void):_id(0),_label(0), _truth(0), _x(0),_y(0),_z(0),_r(0),_phi(0),_theta(0),_rz(0),_distance(0) { for (int i=0;i<NEIGHBOURS;i++) _neighbour[i] = _recall[i] = -1; }
    Point(const Point &p);
    Point(double x, double y, double z, int id=-1, int label=-1, int truth=-1);
    Point(float x, float y, float z, int id=-1, int label=-1, int truth=-1);
    Point operator+(const Point p) const { return Point(_x+p._x,_y+p._y,_z+p._z);}
    Point operator-(const Point p) const { return Point(_x-p._x,_y-p._y,_z-p._z);}
    bool operator<(const Point p) const { return _id<p._id;}
    static bool sortRad(const Point &a,const Point &b);
    static bool sortRz(const Point &a,const Point &b);
    static bool sortDist(const Point &a,const Point &b);
    static bool sortId(const Point &a,const Point &b);
    static bool sortRecall(const Point &a,const Point &b);
    static double angleBetween(const Point &a,const Point &b,const Point &c);
    static double dot(const Point &a,const Point &b);
    static Point circleCenter(const Point &p1,const Point &p2,const Point &p3);
    static double circleRadius(const Point &p1,const Point &p2,const Point &p3);
    static bool comparison(const Point &a,const Point &b);
    double distance(const Point &a);
    static double distance(const Point &a, const Point &b);
    static int classifyAPoint(std::vector<Point> &arr, int k, Point &p, int label);
    inline double x() const {return _x;}
    inline double y() const {return _y;}
    inline double z() const {return _z;}
    inline double r() const {return _r;}
    inline double theta() const {return _theta;}
    inline double phi() const {return _phi;}
    inline double rz() const {return _rz;}
    inline int id() const {return _id;}
    inline int label() const {return _label;}
    inline int truth() const {return _truth;}
    inline int neighbour(int i=0) const {if (i<NEIGHBOURS && i>=0) return _neighbour[i]; else return -1;}
    inline int* neighbours() {return _neighbour;}
    inline double recall(int i=0) const {if (i<NEIGHBOURS && i>=0) return _recall[i]; else return -1;}
    inline double* recalls() {return _recall;}
    //inline void setx(double x) { _x = x;}
    //inline void sety(double y) { _y = y;}
    //inline void setz(double z) { _z = z;}
    inline void setlabel(int label) { _label = label;}
    inline void settruth(int label) { _truth = label;}
    inline void setid(int id) { _id = id;}
    inline void setneighbour(int neighbour, int i=0) { if (i<NEIGHBOURS && i>=0) _neighbour[i] = neighbour;}
    inline void setrecall(double recall, int i=0) { if (i<NEIGHBOURS && i>=0) _recall[i] = recall;}
};

// Used to sort an array of points by increasing
// order of distance from origin
inline
bool Point::sortRad(const Point &a,const Point &b)
{
    return (a._r < b._r);
}

// Used to sort an array of points by increasing
// order of distance from origin in xy plane
inline
bool Point::sortRz(const Point &a,const Point &b)
{
    return (a._rz < b._rz);
}

// Used to sort an array of points by distance
inline
bool Point::sortDist(const Point &a,const Point &b)
{
    return (a._distance < b._distance);
}

// Used to sort an array of points by increasing id
inline
bool Point::sortId(const Point &a,const Point &b)
{
    return (a._id < b._id);
}

// Used to sort an array of points by decreasing recall
inline
bool Point::sortRecall(const Point &a,const Point &b)
{
    return (a._recall[0] > b._recall[0]);
}

// Calculate the angle between two vectors defined by (a,b) and (a,c)
inline
double Point::angleBetween(const Point &a,const Point &b,const Point &c)
{
    Point v1 = b - a;
    Point v2 = c - a;
    double angle = acos(dot(v1,v2));
    return angle;
}

// Calculate the scalar product between two vectors a and b
inline
double Point::dot(const Point &a,const Point &b)
{
    double r1 = sqrt(a._x*a._x + a._y*a._y + a._z*a._z);
    double r2 = sqrt(b._x*b._x + b._y*b._y + b._z*b._z);
    double product = a._x*b._x + a._y*b._y + a._z*b._z;
    return product / (r1 * r2);
}

// Calculate the circle center through 3 points
inline
Point Point::circleCenter(const Point &p1,const Point &p2,const Point &p3)
{
    Point center;
    center._x = 0.0;
    center._y = 0.0;
    center._z = 0.0;
    
    double offset = pow(p2._x,2) + pow(p2._y,2);
    double bc =   ( pow(p1._x,2) + pow(p1._y,2) - offset )/2.0;
    double cd =   (offset - pow(p3._x, 2) - pow(p3._y, 2))/2.0;
    double det =  (p1._x - p2._x) * (p2._y - p3._y) - (p2._x - p3._x)* (p1._y - p2._y);
    
    if (abs(det) < 0.0000001) { return center; }
    
    double idet = 1/det;
    
    center._x =  (bc * (p2._y - p3._y) - cd * (p1._y - p2._y)) * idet;
    center._y =  (cd * (p1._x - p2._x) - bc * (p2._x - p3._x)) * idet;
    
    return center;
}

// Calculate the circle radius through 3 points
inline
double Point::circleRadius(const Point &p1,const Point &p2,const Point &p3)
{
    Point center = circleCenter(p1,p2,p3);
    double radius = sqrt( pow(p2._x - center._x,2) + pow(p2._y-center._y,2));
    return radius;
}


// Used to sort an array of points by increasing
// order of distance
inline
bool Point::comparison(const Point &a,const Point &b)
{
    return (a._distance < b._distance);
}

inline
double Point::distance(const Point &a)
{
    _distance =  sqrt((a._x - _x) * (a._x - _x) +
                      (a._y - _y) * (a._y - _y) +
                      (a._z - _z) * (a._z - _z));
    return _distance;
}

inline
double Point::distance(const Point &a, const Point &b)
{
    double d =  sqrt((a._x - b._x) * (a._x - b._x) +
                     (a._y - b._y) * (a._y - b._y) +
                     (a._z - b._z) * (a._z - b._z));
    return d;
}

// This function finds classification of point p using
// k nearest neighbour algorithm. It assumes only two
// groups and returns 0 if p belongs to group 0, else
// 1 (belongs to group label).
inline
int Point::classifyAPoint(std::vector<Point> &arr, int k, Point &p, int label=1)
{
    // Fill distances of all points from p
    for (auto it=arr.begin(); it != arr.end(); ++it) it->_distance = p.distance(*it);

    // Sort the Points by distance from p
    std::sort(arr.begin(), arr.end(), comparison);
    
    // Now consider the first k elements and only
    // two groups
    int freq1 = 0;     // Frequency of group 0
    int freq2 = 0;     // Frequency of group 1
    for (int i = 0; i < k; i++)
    {
        if (arr[i]._label == 0)
            freq1++;
        else if (arr[i]._label == label)
            freq2++;
    }
    
    return (freq1 > freq2 ? 0 : 1);
}

static unsigned long nr, nd, np, nt, nx, n1, n2, n3, n4;

class Tracker {
private:
    static void print(std::vector<int> const &input);
    static bool sortFunc( const std::vector<int>& p1,const std::vector<int>& p2 );
    static double* Recall2(Point &p1, Point &p2);
    static double* Recall3(Point &p1, Point &p2, Point &p3);
public:
    Tracker() {}
    static int findTracks(int nhits, float *x, float *y, float *z, int* labels);
    static long selectPoints(std::vector<Point> &points, std::vector<Point> &selection, double rmin, double rmax, double zmin, double zmax);
    static long selectPoints(std::vector<Point> &points, std::vector<Point> &selection, Point &ref, double deltar, double deltathe, double distance);
    static void kNearestNeighbour(std::vector<Point> &points);
    static bool checkTracklet(Point &p0,Point &p1);
    static bool checkTracklet(Point &p0,Point &p1, Point &p2);
};

#endif
