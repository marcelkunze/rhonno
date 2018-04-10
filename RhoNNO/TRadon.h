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
typedef struct { double kappa,phi,gamma,sigma,density,x,y,z; } RADON;
double  radon_hit_density(RADON *t);

struct Point {
    float x_, y_, z_, d_;
    Point(float x, float y, float z) : x_(x), y_(y), z_(z), d_(sqrt(x_*x_ + y_*y_ + z_*z_)) {}
    bool operator<(Point const& other) {
        return d_ < other.d_;
    }
    float x() { return x_;}
    float y() { return y_;}
    float z() { return z_;}
    float d() { return d_;}
};

class TRadon : public TObject {
public:
    TRadon();
    std::vector<RADON>& Transform(std::vector<Point> &hits);
    void GenerateTrack(std::vector<Point> &hits, int np, double delta, double radius, double phi, double gamma, double sigma=0.0);
    void Draw (Option_t *option="");
    ~TRadon();
private:
    double getEta_g(RADON *t);
    double getTau_i(RADON *t);
    double getZ_g(RADON *t);
    double radon_hit_density(RADON *t);
    TNtuple *nt1;
    TObjArray Hlist;
    std::vector<RADON> nt2;
public:
    ClassDef(TRadon,1)    // Fuzzy Radon transform
};


#endif
