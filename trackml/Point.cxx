#include "Point.h"
#include <cmath>

Point::Point(double x, double y, double z, int id, int label, int truth)
{
    _id = id;
    _label = label;
    _truth = truth;
    _weight = 0.0;
    _volume = 0;
    _layer = 0;
    _module = 0;
    _x = x;
    _y = y;
    _z = z;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
    for (int i=0; i<NEIGHBOURS;i++) _neighbour[i] = _recall[i] = -1;
}

Point::Point(float x, float y, float z, int id, int label, int truth)
{
    _id = id;
    _label = label;
    _truth = truth;
    _weight = 0.0;
    _volume = 0;
    _layer = 0;
    _module = 0;
    _x = x;
    _y = y;
    _z = z;
    _r = sqrt(_x*_x+_y*_y+_z*_z);
    _rz = sqrt(_x*_x+_y*_y);
    _phi = atan2(_y,_x);
    _theta = acos(z/_r);
    _distance = 0.0;
    for (int i=0; i<NEIGHBOURS;i++) _neighbour[i] = _recall[i] = -1;
}

Point::Point(const Point &p)
{
    _id = p._id;
    _label = p._label;
    _truth = p._truth;
    _weight = p._weight;
    _volume = p._volume;
    _layer = p._layer;
    _module = p._module;
    _x = p._x;
    _y = p._y;
    _z = p._z;
    _r = p._r;
    _rz = p._rz;
    _phi = p._phi;
    _theta = p._theta;
    _distance = p._distance;
    for (int i=0; i<NEIGHBOURS;i++) {
        _neighbour[i] = p._neighbour[i];
        _recall[i] = p._recall[i];
    }
}

