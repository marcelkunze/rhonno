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
#define DKAPPA 0.125
#define NKAPPA 25
#define DPHI M_PI/30.
#define NPHI 30
#define SIGMA 0.001
#define TAUMAX 0.5*M_PI
#define THRESHOLD 100000.

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

TRadon::TRadon() : Hlist(0)
{
    nt1 = new TNtuple("RadonTransform","Radon Transform","kappa:phi:gamma:sigma:density:x:y:z");
    
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
    delete nt1;
}

/*
 * Do a RADON transform of coordinates (all lengths in meter)
 * (1 GeV/c == 2.22m Radius at B=1.5T)
 */
std::vector<RADON>& TRadon::Transform(std::vector<Point> &hits)
{
    
    double kappa,gamma,phi,density,sigma,d;
    long k=0,g=0,p=0,i=0;
    
    sigma = SIGMA;
    gamma = 0.0;
    for (g=0;g<NGAMMA;g++,gamma += DGAMMA) {
        kappa = -1.125;
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
                density=d=0.0;
                vector<Point>::iterator it;
                for(it = hits.begin(); it != hits.end(); it++)    {
                    Point point=*it;
                    t.x = point.x();
                    t.y = point.y();
                    t.z = point.z();
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
                    nt2.push_back(t);  // Note the candidate track parameters
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

void TRadon::GenerateTrack(std::vector<Point> &hits, int np, double delta, double radius, double phi, double gamma, double sigma) {
    default_random_engine generator;
    double tau = 0.025;
    for (int i=0; i<np; i++,tau+=delta)
    {
        Float_t X,Y,Z;
        X = radius * ( sin(phi + (signum(radius)) * tau) - sin(phi));
        Y = radius * (-cos(phi + (signum(radius)) * tau) + cos(phi));
        Z = gamma * tau;
        if (sigma > 0.0) {
            normal_distribution<float> distribution0(X,sigma);
            X = distribution0(generator);
            normal_distribution<float> distribution1(Y,sigma);
            Y = distribution1(generator);
            normal_distribution<float> distribution2(Z,sigma);
            Z = distribution2(generator);
        }
        hits.push_back(Point(X,Y,Z));
    }
}

void TRadon::Draw(Option_t *option) {
    int i = 0;
    vector<RADON>::iterator radon;
    for(radon = nt2.begin(); radon != nt2.end(); radon++)    {
        RADON t=*radon;;
        double kappa = t.kappa;
        double phi   = t.phi;
        double gamma = t.gamma;
        double density = t.density;
        double radius = 1./kappa;
        printf("\nTrack candidate #%d r:%f k:%f p:%f g:%f : %f  ",i++,radius,kappa,phi,gamma,density);
        
        if (density > THRESHOLD) {
            printf("Density > %f", THRESHOLD);
            std::vector<Point> nt3;
            GenerateTrack(nt3,25,0.025,radius,phi,gamma);
            // Draw the track candidates
            TPolyMarker3D *hitmarker = new TPolyMarker3D(nt3.size());
            TPolyLine3D *connector = new TPolyLine3D(nt3.size());
            
            vector<Point>::iterator it;
            for(it = nt3.begin(); it != nt3.end(); it++)    {
                static int j = 0;
                Point p=*it;
                hitmarker->SetPoint(j,p.x(),p.y(),p.z());
                hitmarker->SetMarkerSize(0.1);
                hitmarker->SetMarkerColor(kRed);
                hitmarker->SetMarkerStyle(kFullDotLarge);
                hitmarker->Draw(option);
                connector->SetPoint(j++, p.x(),p.y(),p.z());
                connector->SetLineWidth(1);
                connector->SetLineColor(kRed);
                connector->Draw(option);
            }
        }
    }
}
