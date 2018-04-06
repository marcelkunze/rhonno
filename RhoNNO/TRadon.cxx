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
#define NGAMMA 20
#define DKAPPA 0.005
#define NKAPPA 100
#define DPHI M_PI/30.
#define NPHI 30
#define SIGMA 0.002
#define TAUMAX 0.5*M_PI
#define THRESHOLD 10000.

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

TRadon::TRadon() : Hlist(0)
{
    nt1 = new TNtuple("RadonTransform","Radon Transform","kappa:phi:gamma:sigma:density:x:y:z");
    nt2 = new TNtuple("RadonTrackParms","Radon Track Parameters","kappa:phi:gamma:sigma:density:x:y:z");
    
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
    
    double kappa,gamma,phi,density,sigma,d;
    long k=0,g=0,p=0,i=0;
    
    
    long ih = hits->GetEntries();
    
    sigma = SIGMA;
    gamma = 0.0;
    for (g=0;g<NGAMMA;g++,gamma += DGAMMA) {
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
                        nt1->Fill(t.kappa,t.phi,t.gamma,t.sigma,d,t.x,t.y,t.z);
                        nhits++;
                    }
                }
                if (density > 1.0 && nhits > 2) {
                    printf("\nDensity: %f",density);
                    TH2D *h = (TH2D *) Hlist[(Int_t) g];
                    h->Fill(1./kappa,phi,density);
                    nt2->Fill(t.kappa,t.phi,t.gamma,t.sigma,t.density,t.x,t.y,t.z);
                }
            }
        }
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
    return (signum(t->kappa)) * atan(nom/denom);
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
    double tau = 0.025;
    for (int i=0; i<np; i++,tau+=delta)
    {
        Float_t X[3];
        X[0] = radius * ( sin(phi + (signum(radius)) * tau) - sin(phi));
        X[1] = radius * (-cos(phi + (signum(radius)) * tau) + cos(phi));
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
    long ih = nt2->GetEntries();
    for (Int_t i=0;i<ih;++i) {
        nt2->GetEvent(i,1);
        Float_t *x=nt2->GetArgs();
        double kappa = x[0];
        double phi   = x[1];
        double gamma = x[2];
        double density = x[4];
        double radius = 1./kappa;
        printf("\nTrack candidate #%d r:%f k:%f p:%f g:%f : %f  ",i,radius,kappa,phi,gamma,density);
        
        if (density > THRESHOLD) {
            // Positive curvature
            printf("Density > %f", THRESHOLD);
            TNtuple *nt3 = new TNtuple("Track","Track","x:y:z");
            GenerateTrack(nt3,25,0.025,radius,phi,gamma);
            // Draw the track candidates
            long ih = nt2->GetEntries();
            TPolyMarker3D *hitmarker = new TPolyMarker3D(ih);
            TPolyLine3D *connector = new TPolyLine3D(ih);
            
            for (Int_t j=0;j<ih-1;++j) {
                nt3->GetEvent(j,1);
                Float_t *x=nt3->GetArgs();
                hitmarker->SetPoint(j, x[0], x[1], x[2]);
                hitmarker->SetMarkerSize(0.1);
                hitmarker->SetMarkerColor(kRed);
                hitmarker->SetMarkerStyle(kFullDotLarge);
                hitmarker->Draw(option);
                connector->SetPoint(j, x[0], x[1], x[2]);
                connector->SetLineWidth(1);
                connector->SetLineColor(kRed);
                connector->Draw(option);
            }
  
            // Negative curvature
            TNtuple *nt4 = new TNtuple("Track","Track","x:y:z");
            GenerateTrack(nt3,25,0.025,-radius,phi,gamma);
            // Draw the track candidates
            long ih2 = nt4->GetEntries();
            TPolyMarker3D *hitmarker2 = new TPolyMarker3D(ih2);
            TPolyLine3D *connector2 = new TPolyLine3D(ih2);
            
            for (Int_t j=0;j<ih-1;++j) {
                nt3->GetEvent(j,1);
                Float_t *x=nt3->GetArgs();
                hitmarker2->SetPoint(j, x[0], x[1], x[2]);
                hitmarker2->SetMarkerSize(0.1);
                hitmarker2->SetMarkerColor(kRed);
                hitmarker2->SetMarkerStyle(kFullDotLarge);
                hitmarker2->Draw(option);
                connector2->SetPoint(j, x[0], x[1], x[2]);
                connector2->SetLineWidth(1);
                connector2->SetLineColor(kRed);
                connector2->Draw(option);
            }
}
    }
}
