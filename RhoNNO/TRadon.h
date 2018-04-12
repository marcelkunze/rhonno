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
#include <TVector3.h>

class TNtuple;

typedef struct { double x,y,z; } HIT;
typedef struct { double kappa,phi,gamma,sigma,density,x,y,z; std::vector<long> index; } RADON;
double  radon_hit_density(RADON *t);

class TRadon : public TObject {
public:
    TRadon(double sigma=0.001, double threshold=10000.);
    std::vector<RADON>& Transform(std::vector<TVector3> &points);
    void GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error=0.0);
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
    std::vector<TVector3> hits;
    std::vector<RADON> rt;
public:
    ClassDef(TRadon,1)    // Fuzzy Radon transform
};


#endif
