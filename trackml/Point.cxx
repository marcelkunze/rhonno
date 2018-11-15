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

treePoint::treePoint(double x, double y, double z, double cx, double cy, double cz, int id, int label, int truth)
{
    _id = id;
    _hitid = -1;
    _label = label;
    _truth = truth;
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

treePoint::treePoint(float x, float y, float z, float cx, float cy, float cz, int id, int label, int truth)
{
    _id = id;
    _hitid = -1;
    _label = label;
    _truth = truth;
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
    _truth = p._truth;
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
