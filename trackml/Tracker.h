#ifndef _TRACKER_H_
#define _TRACKER_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#define TRACKML

#ifdef TRACKML
#define NETFILE2 "/Users/marcel/workspace/rhonno/trackml/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/trackml/XMLP3.net"
#define NETFILE4 "/Users/marcel/workspace/rhonno/trackml/XMLP4.net"
#define TRACKLET 3
#define NEIGHBOURS 5
#define THRESHOLD 0.985
#define DISTANCE 0.5
#define DELTAR   1.0
#define DELTAPHI 0.1
#else
#define NETFILE2 "/Users/marcel/workspace/rhonno/RhoNNO/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/RhoNNO/XMLP3.net"
#define NETFILE4 "/Users/marcel/workspace/rhonno/RhoNNO/XMLP3.net"
#define TRACKLET 2
#define THRESHOLD 0.7
#define DISTANCE 1.0
#define DELTAR   1.0
#define DELTAPHI 0.04
#endif

#define VERBOSE true

#include <cmath>
#include <algorithm>

class Point
{
private:
    int _id;             // Hit id of point
    int _val;            // Group of point
    double _x, _y, _z;     // Cartesian coordinate of point
    double _r,_phi,_theta; // Spherical coordinates of point
    double _distance;    // Distance from test point
    int _neighbour;
    int _nextneighbour;
public:
    Point(void):_id(0),_val(0),_x(0),_y(0),_z(0),_r(0),_phi(0),_theta(0),_distance(0) {}
    Point(const Point &p);
    Point(double x, double y, double z, int id=-1, int val=-1);
    Point(float x, float y, float z, int id, int val);
    Point operator+(const Point p) const { return Point(_x+p._x,_y+p._y,_z+p._z);}
    Point operator-(const Point p) const { return Point(_x-p._x,_y-p._y,_z-p._z);}
    static bool sortDist(const Point &a,const Point &b);
    static bool sortId(const Point &a,const Point &b);
    static double angleBetween(const Point &a,const Point &b,const Point &c);
    static double dot(const Point &a,const Point &b);
    static Point circleCenter(const Point &p1,const Point &p2,const Point &p3);
    static double circleRadius(const Point &p1,const Point &p2,const Point &p3);
    static bool comparison(const Point &a,const Point &b);
    double distance(const Point &a);
    static double distance(const Point &a, const Point &b);
    static int classifyAPoint(Point arr[], int n, int k, Point p);
    inline double x() {return _x;}
    inline double y() {return _y;}
    inline double z() {return _z;}
    inline double r() {return _r;}
    inline double theta() {return _theta;}
    inline double phi() {return _phi;}
    inline int id() {return _id;}
    inline int val() {return _val;}
    inline int neighbour() {return _neighbour;}
    inline int nextneighbour() {return _nextneighbour;}
    inline void setx(double x) { _x = x;}
    inline void sety(double y) { _y = y;}
    inline void setz(double z) { _z = z;}
    inline void setval(int val) { _val = val;}
    inline void setid(int id) { _id = id;}
    inline void setneighbour(int neighbour) { _neighbour = neighbour;}
    inline void setnextneighbour(int nextneighbour) { _nextneighbour = nextneighbour;}
};

// Used to sort an array of points by increasing
// order of distance from origin
inline
bool Point::sortDist(const Point &a,const Point &b)
{
    return (a._r < b._r);
}

// Used to sort an array of points by increasing
// order of distance from origin
inline
bool Point::sortId(const Point &a,const Point &b)
{
    return (a._id < b._id);
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
    double d =  sqrt((a._x - _x) * (a._x - _x) +
                     (a._y - _y) * (a._y - _y) +
                     (a._z - _z) * (a._z - _z));
    return d;
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
// 1 (belongs to group 1).
inline
int Point::classifyAPoint(Point arr[], int n, int k, Point p)
{
    // Fill distances of all points from p
    for (int i = 0; i < n; i++) arr[i]._distance = p.distance(arr[i]);

    // Sort the Points by distance from p
    std::sort(arr, arr+n, comparison);
    
    // Now consider the first k elements and only
    // two groups
    int freq1 = 0;     // Frequency of group 0
    int freq2 = 0;     // Frequency of group 1
    for (int i = 0; i < k; i++)
    {
        if (arr[i]._val == 0)
            freq1++;
        else if (arr[i]._val == 1)
            freq2++;
    }
    
    return (freq1 > freq2 ? 0 : 1);
}

class Tracker {
private:
    static void print(std::vector<int> const &input);
    static bool sortFunc( const std::vector<int>& p1,const std::vector<int>& p2 );
    static double* Recall2(Point &p1, Point &p2);
    static double* Recall3(Point &p1, Point &p2, Point &p3);
    static double* Recall4(Point &p1, Point &p2, Point &p3, Point &p4);
public:
    Tracker() {}
    static int findTracks(int nhits, float *x, float *y, float *z, int* labels);
};

#endif
