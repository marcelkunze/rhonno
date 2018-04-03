// radon_function
//
// various functions for performing the
// fuzzy Radon transform

#include <stdlib.h>
#include <math.h>
#include <TNtuple.h>
#include <TH2D.h>
#include "RHONNO/radon.h"

#define NHIST2 20

TRadon::TRadon()
{
    nt1 = new TNtuple("Radon Transform","Radon Transform","kappa:phi:gamma:sigma:density:x:y:z");
    nt2 = new TNtuple("Radon Coordinates","Radon Coordinates","x:y:z");

    h2[0]=TH2D("Gamma=0","Radon density in r/Phi (Gamma=0)",40,0.0,4.,30,0.0,3.0);
    h2[1]=TH2D("Gamma=1","Radon density in r/Phi (Gamma=1)",40,0.0,4.,30,0.0,3.0);
    h2[2]=TH2D("Gamma=2","Radon density in r/Phi (Gamma=2)",40,0.0,4.,30,0.0,3.0);
    h2[3]=TH2D("Gamma=3","Radon density in r/Phi (Gamma=3)",40,0.0,4.,30,0.0,3.0);
    h2[4]=TH2D("Gamma=4","Radon density in r/Phi (Gamma=4)",40,0.0,4.,30,0.0,3.0);
    h2[5]=TH2D("Gamma=5","Radon density in r/Phi (Gamma=5)",40,0.0,4.,30,0.0,3.0);
    h2[6]=TH2D("Gamma=6","Radon density in r/Phi (Gamma=6)",40,0.0,4.,30,0.0,3.0);
    h2[7]=TH2D("Gamma=7","Radon density in r/Phi (Gamma=7)",40,0.0,4.,30,0.0,3.0);
    h2[8]=TH2D("Gamma=8","Radon density in r/Phi (Gamma=8)",40,0.0,4.,30,0.0,3.0);
    h2[9]=TH2D("Gamma=9","Radon density in r/Phi (Gamma=9)",40,0.0,4.,30,0.0,3.0);
    h2[10]=TH2D("Gamma=10","Radon density in r/Phi (Gamma=10)",40,0.0,4.,30,0.0,3.0);
    h2[11]=TH2D("Gamma=11","Radon density in r/Phi (Gamma=11)",40,0.0,4.,30,0.0,3.0);
    h2[12]=TH2D("Gamma=12","Radon density in r/Phi (Gamma=12)",40,0.0,4.,30,0.0,3.0);
    h2[13]=TH2D("Gamma=13","Radon density in r/Phi (Gamma=13)",40,0.0,4.,30,0.0,3.0);
    h2[14]=TH2D("Gamma=14","Radon density in r/Phi (Gamma=14)",40,0.0,4.,30,0.0,3.0);
    h2[15]=TH2D("Gamma=15","Radon density in r/Phi (Gamma=15)",40,0.0,4.,30,0.0,3.0);
    h2[16]=TH2D("Gamma=16","Radon density in r/Phi (Gamma=16)",40,0.0,4.,30,0.0,3.0);
    h2[17]=TH2D("Gamma=17","Radon density in r/Phi (Gamma=17)",40,0.0,4.,30,0.0,3.0);
    h2[18]=TH2D("Gamma=18","Radon density in r/Phi (Gamma=18)",40,0.0,4.,30,0.0,3.0);
    h2[19]=TH2D("Gamma=19","Radon density in r/Phi (Gamma=19)",40,0.0,4.,30,0.0,3.0);
}

/*
 * Do a RADON transform of coordinates (all lengths in meter)
 * (1 GeV/c == 2.22m Radius at B=1.5T)
 */
void TRadon::transform(TNtuple *hits)
//long ih,float *x,float *y,float *z)
{
    
    float maxdensity=0,maxkappa=0,maxphi=0,maxgamma=0;
    float kappa,gamma,phi,radius,density,sigma,d,tau;
    long k=0,g=0,p=0,i=0;
    
  
    long ih = hits->GetEntries();
    
    sigma = 0.001;
    for (g=0,gamma=0.;g<15;g++,gamma+=1.) {
        for (k=0,kappa=0.25;k<20;k++,kappa+=0.25) {
            for (p=0,phi=0.;p<30;p++,phi+=M_PI/30.) {
                for (i=0,density=0.0;i<ih;i++) {
                    RADON t;
                    Float_t X[7];
                    X[0] = t.kappa = kappa;
                    X[1] = t.phi   = phi;
                    X[2] = t.gamma = gamma;
                    X[3] = t.sigma = sigma;
                    hits->GetEvent(i,1);
                    Float_t *x=hits->GetArgs();
                    X[4] = t.x = x[0];
                    X[5] = t.y = x[1];
                    X[6] = t.z = x[2];
                    d =  radon_hit_density(&t);
                    density += d;
                    t.density = density;
                    if (d > 0.001) nt1->Fill(X);
                }
                if (density > 0.001 ) {
                    printf("\ng:%f k:%f p:%f : %f",gamma,kappa,phi,density);
                    h2[g].Fill(1./kappa,phi,density);
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
    printf("\nMax. Values g:%f k:%f p:%f : %f",maxgamma,maxkappa,maxphi,maxdensity);
    /*
     * Create a tuple with RADON coordinates
     */
    for (tau = 0.0; tau <= M_PI; tau = tau + 0.0125)
    {
        Float_t X[3];
        radius = 1./maxkappa;
        X[0] = radius * ( sin(maxphi + tau) - sin(maxphi));
        X[1] = radius * (-cos(maxphi + tau) + cos(maxphi));
        X[2] = maxgamma * tau;
        nt2->Fill(X);
    }
}

float   TRadon::getEta_g(RADON *t)
{
    float  kappa,phi,x,y,sp,cp;
    kappa = t->kappa;
    phi   = t->phi;
    x     = t->x;
    y     = t->y;
    sp    = sin(phi);
    cp    = cos(phi);
    return sqrt((((kappa * x) + sp) * ((kappa * x) + sp)) + (((kappa * y) - cp) * ((kappa * y) - cp)));
}

float   TRadon::getTau_g(RADON *t)
{
    float  kappa,phi,x,y,sp,cp;
    kappa = t->kappa;
    phi   = t->phi;
    x     = t->x;
    y     = t->y;
    sp    = sin(phi);
    cp    = cos(phi);
    return atan(((y * sp) + (x * cp)) / ((1 / kappa) + (x * sp) - (y * cp)));
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
