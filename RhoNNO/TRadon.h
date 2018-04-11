//
//  radon.h
//  RhoNNOLib
//
//  Created by Marcel Kunze on 03.04.18.
//  Copyright © 2018 Marcel Kunze. All rights reserved.
//

#ifndef TRadon_h
#define TRadon_h

#include <TObject.h>
#include <TObjArray.h>

class TNtuple;

typedef struct { double x,y,z; } HIT;
typedef struct { double kappa,phi,gamma,sigma,density,x,y,z; std::vector<long> index; } RADON;
double  radon_hit_density(RADON *t);

struct Point {
    float x_, y_, z_, d_;
    Point() : x_(0), y_(0), z_(0) {}
    Point(float x, float y, float z) : x_(x), y_(y), z_(z), d_(sqrt(x_*x_ + y_*y_ + z_*z_)) {}
    bool operator< (Point const &other) const { return less(*this, other); }
    bool operator> (Point const &other) const { return greater(*this, other); }
    bool operator== (Point const &other) const { return equals(*this, other); }
    float x() const { return x_;}
    float y() const { return y_;}
    float z() const { return z_;}
    float d() const { return d_;}
    bool greater(Point const & a, Point const & b) const { return a.d_ > b.d_; }
    bool less(Point const & a, Point const & b) const { return a.d_ < b.d_; }
    bool equals(Point const & a, Point const & b) const { return a.d_ == b.d_; }
};

class TRadon : public TObject {
public:
    TRadon(double sigma=0.001, double threshold=10000.);
    std::vector<RADON>& Transform(std::vector<Point> &points);
    void GenerateTrack(std::vector<Point> &points, int np, double delta, double radius, double phi, double gamma, double error=0.0);
    void Draw (Option_t *option="");
    ~TRadon();
private:
    double getEta_g(RADON *t);
    double getTau_i(RADON *t);
    double getZ_g(RADON *t);
    double radon_hit_density(RADON *t);
    double sigma, threshold;
    TNtuple *nt1, *nt2;
    TObjArray Hlist;
    std::vector<Point> hits;
    std::vector<RADON> rt;
public:
    ClassDef(TRadon,1)    // Fuzzy Radon transform
};


#endif
