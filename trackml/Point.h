#ifndef _POINT_H_
#define _POINT_H_
// Defines Point class for Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

class Point
{
protected:
    int _id;               // Hit id
    float _x, _y, _z;      // Cartesian coordinates
    float _rz,_phi,_theta; // Cylinder coordinates
    float _r;              // Spherical coordinates
    float _distance;       // Distance from test point
public:
    Point(void):_id(0),_x(0),_y(0),_z(0),_rz(0),_phi(0),_theta(0),_r(0),_distance(0) {}
    Point(const Point &p);
    Point(double x, double y, double z, int id=-1);
    Point(float x, float y, float z, int id=-1);
    Point operator+(const Point p) const { return Point(_x+p._x,_y+p._y,_z+p._z);}
    Point operator-(const Point p) const { return Point(_x-p._x,_y-p._y,_z-p._z);}
    bool operator<(const Point p) const { return _id<p._id;}
    static bool sortRad(const Point &a,const Point &b);
    static bool sortRz(const Point &a,const Point &b);
    static bool sortDist(const Point &a,const Point &b);
    static bool sortId(const Point &a,const Point &b);
    static float angleBetween(const Point &a,const Point &b,const Point &c);
    static float dot(const Point &a,const Point &b);
    static Point circleCenter(const Point &p1,const Point &p2,const Point &p3);
    static float circleRadius(const Point &p1,const Point &p2,const Point &p3);
    static bool comparison(const Point &a,const Point &b);
    float distance() { return _distance;}
    float distance(const Point &a);
    static float distance(const Point &a, const Point &b);
    inline float x() const {return _x;}
    inline float y() const {return _y;}
    inline float z() const {return _z;}
    inline float r() const {return _r;}
    inline float theta() const {return _theta;}
    inline float phi() const {return _phi;}
    inline float rz() const {return _rz;}
    inline int id() const {return _id;}
    inline void setid(int id) { _id = id;}
};

class treePoint : public Point
{
private:
    int _label;             // Label
    int _truth;             // True id
    int _layer;             // layer index (0...47)
    int _twin;              // Double hit index
    std::vector<int> _adjacent; // id of adjacent points
    std::vector<float> _recall;  // and recall values
public:
    treePoint() : Point {}, _label(0), _truth(0) {};
    treePoint(const treePoint &p);
    treePoint(double x, double y, double z, int id=-1, int label=-1, int truth=-1);
    treePoint(float x, float y, float z, int id=-1, int label=-1, int truth=-1);
    static bool sortRecall(const treePoint &a,const treePoint &b);
    static int classifyAPoint(std::vector<treePoint> &arr, int k, treePoint &p, int label);
    inline int label() const {return _label;}
    inline int truth() const {return _truth;}
    inline int layer() const {return _layer;}
    inline int twin() const {return _twin;}
    inline int neighbour(unsigned int i) const {if (i<_adjacent.size()) return _adjacent[i]; else return -1;}
    inline std::vector<int> &neighbours() {return _adjacent;}
    inline int recall(unsigned int i) const {if (i<_recall.size()) return _recall[i]; else return -1;}
    inline std::vector<float> &recalls() {return _recall;}
    inline void setlabel(int label) { _label = label;}
    inline void settruth(int truth) { _truth = truth;}
    inline void setlayer(int layer) { _layer = layer;}
    inline void settwin(int twin) { _twin = twin;}
    inline void setneighbour(int neighbour, double recall=-1.0) { _adjacent.push_back(neighbour); _recall.push_back(recall);}
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
bool treePoint::sortRecall(const treePoint &a,const treePoint &b)
{
    return (a.recall(0) > b.recall(0));
}

// Calculate the angle between two vectors defined by (a,b) and (a,c)
inline
float Point::angleBetween(const Point &a,const Point &b,const Point &c)
{
    Point v1 = b - a;
    Point v2 = c - a;
    float angle = acos(dot(v1,v2));
    return angle;
}

// Calculate the scalar product between two vectors a and b
inline
float Point::dot(const Point &a,const Point &b)
{
    float r1 = sqrt(a._x*a._x + a._y*a._y + a._z*a._z);
    float r2 = sqrt(b._x*b._x + b._y*b._y + b._z*b._z);
    float product = a._x*b._x + a._y*b._y + a._z*b._z;
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
    
    if (fabs(det) < 0.0000001) { return center; }
    
    double idet = 1/det;
    
    center._x =  (bc * (p2._y - p3._y) - cd * (p1._y - p2._y)) * idet;
    center._y =  (cd * (p1._x - p2._x) - bc * (p2._x - p3._x)) * idet;
    
    return center;
}

// Calculate the circle radius through 3 points
inline
float Point::circleRadius(const Point &p1,const Point &p2,const Point &p3)
{
    Point center = circleCenter(p1,p2,p3);
    float radius = sqrt( pow(p2._x - center._x,2) + pow(p2._y-center._y,2));
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
float Point::distance(const Point &a)
{
    _distance =  sqrt((a._x - _x) * (a._x - _x) +
                      (a._y - _y) * (a._y - _y) +
                      (a._z - _z) * (a._z - _z));
    return _distance;
}

inline
float Point::distance(const Point &a, const Point &b)
{
    float d =  sqrt((a._x - b._x) * (a._x - b._x) +
                     (a._y - b._y) * (a._y - b._y) +
                     (a._z - b._z) * (a._z - b._z));
    return d;
}

// This function finds classification of point p using
// k nearest neighbour algorithm. It assumes only two
// groups and returns 0 if p belongs to group 0, else
// 1 (belongs to group label).
inline
int treePoint::classifyAPoint(std::vector<treePoint> &arr, int k, treePoint &p, int label=1)
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

#endif

