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

TRadon::TRadon() : Hlist(0)
{
    nt1 = new TNtuple("RadonTransform","Radon Transform","kappa:phi:gamma:sigma:density:x:y:z");
    nt2 = new TNtuple("RadonCoordinates","Radon Coordinates","x:y:z");

    for (int i=0;i<30;i++) {
        string name = "Gamma=" + to_string(i);
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
    float kappa,gamma,phi,radius,density,sigma,d,tau;
    long k=0,g=0,p=0,i=0;
    
  
    long ih = hits->GetEntries();
    
    sigma = 0.001;
    for (g=0,gamma=0.;g<30;g++,gamma+=1.) {
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
    printf("\nMax. Values g:%f k:%f p:%f : %f\n",maxgamma,maxkappa,maxphi,maxdensity);
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
    return nt2;
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
