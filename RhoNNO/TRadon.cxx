// radon_function
//
// various functions for performing the
// fuzzy Radon transform

#include <stdlib.h>
#include <math.h>
#include <TNtuple.h>
#include <TH2D.h>
#include "RhoNNO/TRadon.h"

#include <string>
#include <iostream>
using namespace std;

ClassImp(TRadon)

#define DGAMMA 0.1
#define DKAPPA 0.25
#define DPHI M_PI/30.
#define SIGMA 0.005

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

TRadon::TRadon() : Hlist(0)
{
    nt1 = new TNtuple("RadonTransform","Radon Transform","kappa:phi:gamma:sigma:density:x:y:z");
    nt2 = new TNtuple("RadonCoordinates","Radon Coordinates","x:y:z");
 
    float gamma = 0.0;
    for (int i=0;i<25;i++,gamma+=DGAMMA) {
        string name = "Gamma=" + to_string(gamma);
        string title = "Radon density in r/Phi " + name;
        Hlist.Add(new TH2D(TString(name),TString(title),40,0.0,4.,30,0.0,3.0));
    }
}

TRadon::~TRadon() {
    Hlist.Write();
    nt1->Write();
    nt2->Write();
    delete nt1;
    delete nt2;
}

/*
 * Do a RADON transform of coordinates (all lengths in meter)
 * (1 GeV/c == 2.22m Radius at B=1.5T)
 */
TNtuple* TRadon::Transform(TNtuple *hits)
{
    
    float maxdensity=0,maxkappa=0,maxphi=0,maxgamma=0;
    float kappa,gamma,phi,density,sigma,d;
    long k=0,g=0,p=0,i=0;
    
    
    long ih = hits->GetEntries();
    
    sigma = SIGMA;
    gamma = 0.;
    for (g=0;g<25;g++) {
        gamma = gamma + DGAMMA;
        kappa = DKAPPA;
        for (k=0;k<20;k++) {
            kappa = kappa + DKAPPA;
            phi = 0.;
            for (p=0;p<30;p++) {
                phi = phi + DPHI;
                RADON t;
                t.kappa = kappa;
                t.phi   = phi;
                t.gamma = gamma;
                t.sigma = sigma;
                for (i=0,density=0.0;i<ih;i++) {
                    hits->GetEvent(i,1);
                    Float_t *x=hits->GetArgs();
                    t.x = x[0];
                    t.y = x[1];
                    t.z = x[2];
                    d =  radon_hit_density(&t);
                    density += d;
                    t.density = density;
                    if (d > 0.001) {
                        //printf("\nHit #%ld: %f %f %f",i,t.x,t.y,t.z);
                        printf("\nHit #%ld | r:%f k:%f p:%f g:%f : %f",i,1./kappa,kappa,phi,gamma,density);
                        nt1->Fill(t.kappa,t.phi,t.gamma,t.sigma,t.density,t.x,t.y,t.z);
                    }
                }
                if (density > 0.001 ) {
                    printf("\nr:%f k:%f p:%f g:%f : %f\n",1./kappa,kappa,phi,gamma,density);
                    TH2D *h = (TH2D *) Hlist[(Int_t) g];
                    h->Fill(1./kappa,phi,density);
                    if (density > maxdensity) {
                        maxdensity = density;
                        maxkappa = kappa;
                        maxphi   = phi;
                        maxgamma = gamma;
                    }
                }
            }
        }
    }
    if (maxkappa>0.0) printf("\nMax. Values r:%f k:%f p:%f g:%f : %f\n",1./maxkappa, maxkappa,maxphi,maxgamma,maxdensity);
    
    /*
     * Create a tuple with RADON coordinates
     */
    GenerateTrack(nt2,30,0.0125,1./maxkappa,maxphi,maxgamma);

    return nt2;
}

float   TRadon::getEta_g(RADON *t)
{
    float  sp,cp;
    sp    = sin(t->phi);
    cp    = cos(t->phi);
    return sqrt((((t->kappa * t->x) + sp) * ((t->kappa * t->x) + sp)) + (((t->kappa * t->y) - cp) * ((t->kappa * t->y) - cp)));
}

float   TRadon::getTau_g(RADON *t)
{
    float  sp,cp;
    sp    = sin(t->phi);
    cp    = cos(t->phi);
    return signum(t->kappa) * atan(((t->y * sp) + (t->x * cp)) / ((1 / t->kappa) + (t->x * sp) - (t->y * cp)));
}

float   TRadon::getZ_g(RADON *t)
{
    return ((t->gamma * getTau_g(t)) - t->z);
}

float   TRadon::radon_hit_density(RADON *t)
{
    float  eta_g,gamma,kappa,sigma,z_g,k2,s2,g2,factor, kappafunction,efunction,radon;
    eta_g = getEta_g(t);
    gamma = t->gamma;
    kappa = t->kappa;
    sigma = t->sigma;
    z_g   = getZ_g(t)  ;
    k2    = kappa * kappa;
    s2    = sigma * sigma;
    g2    = gamma * gamma;
    
    factor = 1. / (2 * M_PI * s2 * TAUMAX);
    kappafunction = kappa / sqrt(eta_g + k2 * g2);
    efunction     = -1. / (2 * s2) * (((eta_g - 1) * (eta_g - 1) / k2) + (eta_g * z_g * z_g / (eta_g + k2 * g2)));
    radon  =   factor * kappafunction * exp(efunction);
    
    return radon;
    
}

/*
 * Create a tuple with Track coordinates
 */
void TRadon::GenerateTrack(TNtuple *nt, int np, float delta, float radius, float phi, float gamma) {
    float tau = delta;
    for (int i=0; i<np; i++,tau+=delta)
    {
        Float_t X[3];
        X[0] = radius * ( sin(phi + (signum(radius) * tau)) - sin(phi));
        X[1] = radius * (-cos(phi + (signum(radius) * tau)) + cos(phi));
        X[2] = gamma * tau;
        nt->Fill(X);
    }
}

