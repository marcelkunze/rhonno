#ifndef _TRACKER_H_
#define _TRACKER_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

//#define TRACKML

#ifdef TRACKML
#define NETFILE2 "/Users/marcel/workspace/rhonno/trackml/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/trackml/XMLP3.net"
#define TRACKLET 2
#define THRESHOLD 0.7
#define DISTANCE 5.0
#define DELTAR   2.0
#define DELTAPHI 0.8
#else
#define NETFILE2 "/Users/marcel/workspace/rhonno/RhoNNO/XMLP2.net"
#define NETFILE3 "/Users/marcel/workspace/rhonno/RhoNNO/XMLP3.net"
#define TRACKLET 2
#define THRESHOLD 0.6
#define DISTANCE 0.2
#define DELTAR   0.2
#define DELTAPHI 0.02
#endif

#define VERBOSE true

#include <cmath>
#include <algorithm>

struct Point
{
    int id;             // Hit id of point
    int val;            // Group of point
    double x, y, z;     // Cartesian coordinate of point
    double r,phi,theta; // Spherical coordinates of point
    double distance;    // Distance from test point
};

// Used to sort an array of points by increasing
// order of distance from origin
inline bool sortDist(const Point &a,const Point &b)
{
    return (a.r < b.r);
}

// Used to sort an array of points by increasing
// order of distance from origin
inline bool sortId(const Point &a,const Point &b)
{
    return (a.id < b.id);
}


// Calculate the angle between two vectors defined by (a,b) and (a,c)
inline
double angleBetween(const Point &a,const Point &b,const Point &c)
{
    double x1 = b.x - a.x;
    double y1 = b.y - a.y;
    double z1 = b.z - a.z;
    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);

    double x2 = c.x - a.x;
    double y2 = c.y - a.y;
    double z2 = c.z - a.z;
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);

    double product =x1*x2 + y1*y2 + z1*z2;
    double angle = abs(product / r1 * r2);
    return angle;
}

// Calculate the circle center through 3 points
inline
Point circleCenter(const Point &p1,const Point &p2,const Point &p3)
{
    Point center;
    center.x = 0.0;
    center.y = 0.0;
    center.z = 0.0;
    
    double offset = pow(p2.x,2) + pow(p2.y,2);
    double bc =   ( pow(p1.x,2) + pow(p1.y,2) - offset )/2.0;
    double cd =   (offset - pow(p3.x, 2) - pow(p3.y, 2))/2.0;
    double det =  (p1.x - p2.x) * (p2.y - p3.y) - (p2.x - p3.x)* (p1.y - p2.y);
    
    if (abs(det) < 0.0000001) { return center; }
    
    double idet = 1/det;
    
    center.x =  (bc * (p2.y - p3.y) - cd * (p1.y - p2.y)) * idet;
    center.y =  (cd * (p1.x - p2.x) - bc * (p2.x - p3.x)) * idet;
    
    return center;
}

// Calculate the circle radius through 3 points
inline
double circleRadius(const Point &p1,const Point &p2,const Point &p3)
{
    Point center = circleCenter(p1,p2,p3);
    double radius = sqrt( pow(p2.x - center.x,2) + pow(p2.y-center.y,2));
    return radius;
}


// Used to sort an array of points by increasing
// order of distance
inline
bool comparison(const Point &a,const Point &b)
{
    return (a.distance < b.distance);
}

inline
double distanceBetween(const Point &a,const Point &b)
{
    double d =  sqrt((a.x - b.x) * (a.x - b.x) +
                     (a.y - b.y) * (a.y - b.y) +
                     (a.z - b.z) * (a.z - b.z));
    return d;
}

// This function finds classification of point p using
// k nearest neighbour algorithm. It assumes only two
// groups and returns 0 if p belongs to group 0, else
// 1 (belongs to group 1).
inline
int classifyAPoint(Point arr[], int n, int k, Point p)
{
    // Fill distances of all points from p
    for (int i = 0; i < n; i++) arr[i].distance = distanceBetween(arr[i],p);

    // Sort the Points by distance from p
    std::sort(arr, arr+n, comparison);
    
    // Now consider the first k elements and only
    // two groups
    int freq1 = 0;     // Frequency of group 0
    int freq2 = 0;     // Frequency of group 1
    for (int i = 0; i < k; i++)
    {
        if (arr[i].val == 0)
            freq1++;
        else if (arr[i].val == 1)
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
public:
    static int findTracks(int nhits, float *x, float *y, float *z, int* labels);
};

#endif
