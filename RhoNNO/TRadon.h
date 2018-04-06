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

#define TAUMAX 0.5*M_PI

typedef struct { double x,y,z; } HIT;
typedef struct { double kappa,phi,gamma,sigma,density,x,y,z; } RADON;
double  radon_hit_density(RADON *t);

class TRadon : public TObject {
public:
    TRadon();
    TNtuple* Transform(TNtuple *hits);
    void GenerateTrack(TNtuple *nt, int np, double delta, double radius, double phi, double gamma, double sigma=0.0);
    void Draw (Option_t *option="");
    ~TRadon();
private:
    double getEta_g(RADON *t);
    double getTau_g(RADON *t);
    double getZ_g(RADON *t);
    double radon_hit_density(RADON *t);
    TNtuple *nt1, *nt2;
    TObjArray Hlist;
public:
    ClassDef(TRadon,1)    // Fuzzy Radon transform
};


#endif
