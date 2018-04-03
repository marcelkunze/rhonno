//
//  radon.h
//  RhoNNOLib
//
//  Created by Marcel Kunze on 03.04.18.
//  Copyright © 2018 Marcel Kunze. All rights reserved.
//

#ifndef radon_h
#define radon_h

#define TAUMAX 0.5*M_PI

typedef struct { float x,y,z; } HIT;
typedef struct { float kappa,phi,gamma,sigma,density,x,y,z; } RADON;
float  radon_hit_density(RADON *t);

class TRadon {
public:
    TRadon();
    void transform(TNtuple *hits);
    ~TRadon();
private:
    float getEta_g(RADON *t);
    float getTau_g(RADON *t);
    float getZ_g(RADON *t);
    float radon_hit_density(RADON *t);
    TH2D h2[NHIST2];
    TNtuple *nt1, *nt2;
};


#endif /* radon_h */
