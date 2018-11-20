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
    //int _hitid;            // TrackML hit id
    float _x,_y,_z;        // Cartesian coordinates
    float _cx,_cy,_cz;     // Directional cosines
    float _rz,_phi,_theta; // Cylinder coordinates
    float _r;              // Spherical coordinates
    float _distance;       // Distance from test point
public:
    Point(void):_id(0),_x(0),_y(0),_z(0),_cx(0),_cy(0),_cz(0),_rz(0),_phi(0),_theta(0),_r(0),_distance(0) {}
    Point(const Point &p);
    Point(double x,double y,double z,double cx=0.,double cy=0.,double cz=0.,int id=-1,int hitid=-1);
    Point(float x,float y,float z,float cx=0.,float cy=0.,float cz=0.,int id=-1,int hitid=-1);
    inline Point operator+(const Point p) const { return Point(_x+p._x,_y+p._y,_z+p._z);}
    inline Point operator-(const Point p) const { return Point(_x-p._x,_y-p._y,_z-p._z);}
    inline double operator*(const Point &p) { return _x*p._x + _y*p._y + _z*p._z; }
    inline Point operator*(double f) { return Point(_x*f, _y*f, _z*f); }
    inline Point scale(double f) { return Point(_x*f, _y*f, _z*f); }
    inline Point cross(Point p) const { return Point(_y*p._z - _z*p._y, _z*p._x - _x*p._z, _x*p._y - _y*p._x);}
    bool operator<(const Point p) const { return _id<p._id;}
    static bool sortByRadius(const Point &a,const Point &b);
    static bool sortByRz(const Point &a,const Point &b);
    static bool sortByDistance(const Point &a,const Point &b);
    static bool sortById(const Point &a,const Point &b);
    static float angleBetween(const Point &a,const Point &b,const Point &c);
    static Point normalize(Point a);
    static float dot(const Point &a,const Point &b);
    static Point cross(const Point &a,const Point &b);
    static void circle(Point&a, Point&b, Point&c, Point&p, double &r);
    static Point circleCenter(const Point &p1,const Point &p2,const Point &p3);
    static float circleRadius(const Point &p1,const Point &p2,const Point &p3);
    static bool comparison(const Point &a,const Point &b);
    static float norm(const Point &a);
    static float norm2(const Point &a);
    float distance() { return _distance;}
    float distance(const Point &a);
    static float distance(const Point &a, const Point &b);
    static float distance3(Point &a,Point &b,Point &c);
    inline float x() const {return _x;}
    inline float y() const {return _y;}
    inline float z() const {return _z;}
    inline float cx() { return _cx;}
    inline float cy() { return _cy;}
    inline float cz() { return _cz;}
    inline float r() const {return _r;}
    inline float theta() const {return _theta;}
    inline float phi() const {return _phi;}
    inline float rz() const {return _rz;}
    inline int id() const {return _id;}
    //inline int hitid() const {return _hitid;}
    inline void setid(int id) { _id = id;}
    //inline void sethitid(int id) { _hitid = id;}
    inline void setcx(float cx) { _cx = cx;}
    inline void setcy(float cy) { _cy = cy;}
    inline void setcz(float cz) { _cz = cz;}
    static Point distBetweenLines(Point &p1, Point &p2, Point &p3, Point &p4);
};

class treePoint : public Point
{
private:
    int _label;             // label
    long long _trackid;     // track id
    int _volume;            // volume
    int _layer;             // layer index (0...47)
    int _module;            // module index (1...3192)
    int _twin;              // Double hit index
    std::vector<int> _adjacent; // id of adjacent points
    std::vector<float> _recall;  // and recall values
public:
    treePoint() : Point {}, _label(0), _trackid(0) {};
    treePoint(const treePoint &p);
    treePoint(double x,double y,double z,double cx=0.,double cy=0.,double cz=0.,int id=-1,int label=-1,long long truth=-1);
    treePoint(float x,float y,float z,float cx=0.,float cy=0.,float cz=0.,int id=-1,int label=-1,long long truth=-1);
    static bool sortRecall(const treePoint &a,const treePoint &b);
    static int classifyAPoint(std::vector<treePoint> &arr, int k, treePoint &p, int label);
    inline int label() const {return _label;}
    inline long long trackid() const {return _trackid;}
    inline int volume() const {return _volume;}
    inline int layer() const {return _layer;}
    inline int module() const {return _module;}
    inline int twin() const {return _twin;}
    inline int neighbour(unsigned int i) const {if (i<_adjacent.size()) return _adjacent[i]; else return -1;}
    inline std::vector<int> &neighbours() {return _adjacent;}
    inline int recall(unsigned int i) const {if (i<_recall.size()) return _recall[i]; else return -1;}
    inline std::vector<float> &recalls() {return _recall;}
    inline void setlabel(int label) { _label = label;}
    inline void settrackid(long long trackid) { _trackid = trackid;}
    inline void setvolume(int volume) { _volume = volume;}
    inline void setlayer(int layer) { _layer = layer;}
    inline void setmodule(int module) { _module = module;}
    inline void settwin(int twin) { _twin = twin;}
    inline void setneighbour(int neighbour, double recall=-1.0) { _adjacent.push_back(neighbour); _recall.push_back(recall);}
    static bool sortByLayer(const treePoint &a,const treePoint &b);
};


// Used to sort an array of points by increasing
// order of distance from origin
inline
bool Point::sortByRadius(const Point &a,const Point &b)
{
    return (a._r < b._r);
}

// Used to sort an array of points by increasing
// order of distance from origin in xy plane
inline
bool Point::sortByRz(const Point &a,const Point &b)
{
    return (a._rz < b._rz);
}

// Used to sort an array of points by distance
inline
bool Point::sortByDistance(const Point &a,const Point &b)
{
    return (a._distance < b._distance);
}

// Used to sort an array of points by increasing id
inline
bool Point::sortById(const Point &a,const Point &b)
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

inline
Point Point::normalize(Point a) {
    Point ret = a*(1./Point::norm(a));
    if (ret.z() < 0) ret = ret*-1;
    return ret;
}

// Calculate the scalar product between two vectors a and b
inline
float Point::dot(const Point &a,const Point &b)
{
    double r1 = sqrt(a._x*a._x + a._y*a._y + a._z*a._z);
    double r2 = sqrt(b._x*b._x + b._y*b._y + b._z*b._z);
    double product = a._x*b._x + a._y*b._y + a._z*b._z;
    double result = product / (r1 * r2);
    if (result>1.0) result = 1.0;
    if (result<-1.0) result = -1.0;
    return result;
}

inline Point Point::cross(const Point &a,const Point &b) {
    return Point(a._y*b._z-a._z*b._y, a._z*b._x-a._x*b._z, a._x*b._y-a._y*b._x);
}

//Find circle with center p, radius r, going through a, b, and c (in xy plane)
inline void Point::circle(Point&a, Point&b, Point&c, Point&p, double &r) {
    p = Point(0.,0.,0.);
    r = 0.0;
    double ax = a.x()-c.x(), ay = a.y()-c.y(), bx = b.x()-c.x(), by = b.y()-c.y();
    double aa = ax*ax + ay*ay, bb = bx*bx + by*by;
    if (ax*by == ay*bx) return;
    double idet = .5/(ax*by - ay*bx);
    double x = (aa*by - bb*ay)*idet;
    double y = (ax*bb - bx*aa)*idet;
    double z = 0;
    r = sqrt(x*x + y*y);
    p = Point(x+c.x(),y+c.y(),z);
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
float Point::norm(const Point &a)
{
    return sqrt(a._x*a._x+a._y*a._y+a._z*a._z);
}

inline
float Point::norm2(const Point &a)
{
    return a._x*a._x+a._y*a._y+a._z*a._z;
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

// Distance from Point c to line defined by a and b
inline
float Point::distance3(Point &a,Point &b,Point &c) {
    const Point x = a-b;
    const Point y = a-c;
    const Point z = c-b;
    float d = norm( cross(x,y)) / norm(z);
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

// Used to sort an array of points by increasing layer
inline
bool treePoint::sortByLayer(const treePoint &a,const treePoint &b)
{
    return (a._layer < b._layer);
}


#endif

