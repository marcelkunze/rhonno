// radon_function
//
// various functions for performing the
// fuzzy Radon transform

#include <stdlib.h>
#include <math.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include "RhoNNO/TRadon.h"

#include <string>
#include <iostream>
#include <random>
using namespace std;

ClassImp(TRadon)

#define DGAMMA 0.1
#define NGAMMA 25
#define DKAPPA 0.1
#define NKAPPA 50
#define DPHI M_PI/30.
#define NPHI 30
#define SIGMA 0.00172
#define TAUMAX 0.5*M_PI

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

TRadon::TRadon() : Hlist(0)
{
    nt1 = new TNtuple("RadonTransform","Radon Transform","kappa:phi:gamma:sigma:density:x:y:z");
    nt2 = new TNtuple("RadonCoordinates","Radon Coordinates","x:y:z");
    
    float gamma = 0.0;
    for (int i=0;i<NGAMMA;i++,gamma+=DGAMMA) {
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
    
    double maxdensity=0,maxkappa=0,maxphi=0,maxgamma=0;
    double kappa,gamma,phi,density,sigma,d;
    long k=0,g=0,p=0,i=0;
    
    
    long ih = hits->GetEntries();
    
    sigma = SIGMA;
    gamma = 0.;
    for (g=0;g<NGAMMA;g++) {
        gamma = gamma + DGAMMA;
        kappa = DKAPPA;
        for (k=0;k<NKAPPA;k++) {
            kappa = kappa + DKAPPA;
            phi = 0.;
            for (p=0;p<NPHI;p++) {
                phi = phi + DPHI;
                RADON t;
                t.kappa = kappa;
                t.phi   = phi;
                t.gamma = gamma;
                t.sigma = sigma;
                long nhits = 0;
                for (i=0,density=d=0.0;i<ih;i++) {
                    hits->GetEvent(i,1);
                    Float_t *x=hits->GetArgs();
                    t.x = x[0];
                    t.y = x[1];
                    t.z = x[2];
                    d =  radon_hit_density(&t);
                    density += d;
                    t.density = density;
                    if (d >1.0) {
                        if (nhits==0) printf("\n\nr:%f k:%f p:%f g:%f",1./kappa,kappa,phi,gamma);
                        printf("\nHit #%ld %f,%f,%f| : %f",i,t.x,t.y,t.z,d);
                        nt1->Fill(t.kappa,t.phi,t.gamma,t.sigma,t.density,t.x,t.y,t.z);
                        nhits++;
                    }
                }
                if (density > 1.0 && nhits > 1) {
                    printf("\nDensity: %f",density);
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
    if (maxkappa>0.0) {
        /*
         * Create a tuple with RADON coordinates
         */
        float radius = 1./maxkappa;
        GenerateTrack(nt2,100,0.0125,radius,maxphi,maxgamma);
        printf("\nMax. Values r:%f k:%f p:%f g:%f : %f\n",radius, maxkappa,maxphi,maxgamma,maxdensity);
        
    }
    
    return nt2;
}

double   TRadon::getEta_g(RADON *t)
{
    double  sp,cp,x2,y2;
    sp = sin(t->phi);
    cp = cos(t->phi);
    x2 = (t->kappa * t->x + sp) * (t->kappa * t->x + sp);
    y2 = (t->kappa * t->y - cp) * (t->kappa * t->y - cp);
    return sqrt(x2+y2);
}

double   TRadon::getTau_i(RADON *t)
{
    double  sp,cp, nom, denom;
    sp  = sin(t->phi);
    cp  = cos(t->phi);
    nom = t->y * sp + t->x * cp;
    denom = (1.0 / t->kappa) + t->x * sp - t->y * cp;
    return signum(t->kappa) * atan(nom/denom);
}

double   TRadon::getZ_g(RADON *t)
{
    return ((t->gamma * getTau_i(t)) - t->z);
}

double   TRadon::radon_hit_density(RADON *t)
{
    double  eta_g,gamma,kappa,sigma,z_g,k2,s2,g2,factor, kappafunction,efunction,radon;
    eta_g = getEta_g(t);
    gamma = t->gamma;
    kappa = t->kappa;
    sigma = t->sigma;
    z_g   = getZ_g(t)  ;
    k2    = kappa * kappa;
    s2    = sigma * sigma;
    g2    = gamma * gamma;
    
    factor = 1. / (2 * M_PI * s2 * TAUMAX);
    kappafunction = abs(kappa) / sqrt(eta_g + k2 * g2);
    efunction     = -1.0 / (2 * s2) * (((eta_g - 1.0) * (eta_g - 1.0) / k2) + (eta_g * z_g * z_g / (eta_g + k2 * g2)));
    radon  =   factor * kappafunction * exp(efunction);
    
    return radon;
    
}

/*
 * Create a tuple with Track coordinates
 */

void TRadon::GenerateTrack(TNtuple *nt, int np, double delta, double radius, double phi, double gamma, double sigma) {
    default_random_engine generator;
    double tau = 0.1;
    for (int i=0; i<np; i++,tau+=delta)
    {
        Float_t X[3];
        X[0] = radius * ( sin(phi + tau) - sin(phi));
        X[1] = radius * (-cos(phi + tau) + cos(phi));
        X[2] = gamma * tau;
        if (sigma > 0.0) {
            normal_distribution<float> distribution0(X[0],sigma);
            X[0] = distribution0(generator);
            normal_distribution<float> distribution1(X[1],sigma);
            X[1] = distribution1(generator);
            normal_distribution<float> distribution2(X[2],sigma);
            X[2] = distribution2(generator);
        }
        nt->Fill(X);
        //        printf("\nHit #%d: %f %f %f",i,X[0],X[1],X[2]);
    }
}

void TRadon::Draw(Option_t *option) {
    // Draw the best track
    long ih = nt2->GetEntries();
    TPolyMarker3D *hitmarker = new TPolyMarker3D(ih);
    TPolyLine3D *connector = new TPolyLine3D(ih);
    
    for (Int_t i=0;i<ih-1;++i) {
        nt2->GetEvent(i,1);
        Float_t *x=nt2->GetArgs();
        hitmarker->SetPoint(i, x[0], x[1], x[2]);
        hitmarker->SetMarkerSize(0.5);
        hitmarker->SetMarkerColor(kRed);
        hitmarker->SetMarkerStyle(kFullDotLarge);
        hitmarker->Draw(option);
        connector->SetPoint(i, x[0], x[1], x[2]);
        connector->SetLineWidth(1);
        connector->SetLineColor(kRed);
        connector->Draw(option);
    }
    
}
