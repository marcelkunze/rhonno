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

typedef struct { float x,y,z; } HIT;
typedef struct { float kappa,phi,gamma,sigma,density,x,y,z; } RADON;
float  radon_hit_density(RADON *t);

class TRadon : public TObject {
public:
    TRadon();
    TNtuple* Transform(TNtuple *hits);
    void GenerateTrack(TNtuple *nt, int np, float delta, float radius, float phi, float gamma);
    ~TRadon();
private:
    float getEta_g(RADON *t);
    float getTau_g(RADON *t);
    float getZ_g(RADON *t);
    float radon_hit_density(RADON *t);
    TNtuple *nt1, *nt2;
    TObjArray Hlist;
public:
    ClassDef(TRadon,1)    // Fuzzy Radon transform
};


#endif
