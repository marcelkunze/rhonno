#include "Point.h"
#include <cmath>

Point::Point(double x, double y, double z, double cx, double cy, double cz, int id, int hitid)
{
    _id = id;
    _hitid = hitid;
    _x = x;
    _y = y;
    _z = z;
    _cx = cx;
    _cy = cy;
    _cz = cz;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
}

Point::Point(float x, float y, float z, float cx, float cy, float cz, int id, int hitid)
{
    _id = id;
    _hitid = hitid;
    _x = x;
    _y = y;
    _z = z;
    _cx = cx;
    _cy = cy;
    _cz = cz;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
}

Point::Point(const Point &p)
{
    _id = p._id;
    _hitid = p._hitid;
    _x = p._x;
    _y = p._y;
    _z = p._z;
    _cx = p._cx;
    _cy = p._cy;
    _cz = p._cz;
    _r = p._r;
    _rz = p._rz;
    _phi = p._phi;
    _theta = p._theta;
    _distance = p._distance;
}

treePoint::treePoint(double x, double y, double z, double cx, double cy, double cz, int id, int hitid, int label, long long trackid)
{
    _id = id;
    _hitid = hitid;
    _label = label;
    _trackid = trackid;
    _volume = -1;
    _layer = -1;
    _module = -1;
    _twin = -1;
    _x = x;
    _y = y;
    _z = z;
    _cx = cx;
    _cy = cy;
    _cz = cz;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
}

treePoint::treePoint(float x, float y, float z, float cx, float cy, float cz, int id, int hitid, int label, long long trackid)
{
    _id = id;
    _hitid = hitid;
    _label = label;
    _trackid = trackid;
    _volume = -1;
    _layer = -1;
    _module = -1;
    _twin = -1;
    _x = x;
    _y = y;
    _z = z;
    _cx = cx;
    _cy = cy;
    _cz = cz;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
}

treePoint::treePoint(const treePoint &p)
{
    _id = p._id;
    _hitid = p._hitid;
    _label = p._label;
    _trackid = p._trackid;
    _volume = p._volume;
    _layer = p._layer;
    _module = p._module;
    _twin = p._twin;
    _x = p._x;
    _y = p._y;
    _z = p._z;
    _cx = p._cx;
    _cy = p._cy;
    _cz = p._cz;
    _r = p._r;
    _rz = p._rz;
    _phi = p._phi;
    _theta = p._theta;
    _distance = p._distance;
    _adjacent.clear();
    _adjacent = p._adjacent;
    _recall = p._recall;
}

Point Point::distBetweenLines(Point &p1, Point &p2, Point &p3, Point &p4)
{
    Point u = p1 - p2;
    Point v = p3 - p4;
    Point w = p2 - p4;
    
    double a = dot(u,u);
    double b = dot(u,v);
    double c = dot(v,v);
    double d = dot(u,w);
    double e = dot(v,w);
    double D = a*c - b*b;
    double sD = D;
    double tD = D;
    
    double SMALL_NUM = 0.00000001;
    double tN(0.0), sN(0.0);
    
    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;       // force using point P0 on segment S1
        sD = 1.0;       // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {               // get the closest points on the infinite lines
        double sN = (b*e - c*d);
        double tN = (a*e - b*d);
        if (sN < 0.0) {   // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) { // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }
    
    if (tN < 0.0)    {        //tc < 0 => the t=0 edge is visible
        tN = 0.0;
        sN = -d;
        sD = a;
        // recompute sc for this edge
        if (-d < 0.0) sN = 0.0;
        if (-d > a) sN = sD;
        else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
            tN = tD;
            sN = (-d + b);
            sD = a;
            // recompute sc for this edge
            if ((-d + b) < 0.0) sN = 0.0;
            if ((-d + b) > a) sN = sD;
        }
    }
    
    // finally do the division to get sc and tc
    double sc(0.0), tc(0.0);
    if(abs(sN) < SMALL_NUM)
        sc = 0.0;
    else
        sc = sN / sD;
    
    if(abs(tN) < SMALL_NUM)
        tc = 0.0;
    else
        tc = tN / tD;

    //get the difference of the two closest points
    Point scu = u.scale(sc);
    Point tcv = v.scale(tc);
    Point dP = w + scu - tcv;  // = S1(sc) - S2(tc)
    double distance = norm(dP);
    Point outV = dP; // vector connecting the closest points
    
    Point closestP1 = p2+scu; // Closest point on object 1
    Point closestP2 = p4+tcv; // Closest point on object 2
    
    Point vertex = (closestP1+closestP2).scale(0.5);
    return vertex;
}
