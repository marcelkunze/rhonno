// radon_function
//
// various functions for performing the
// fuzzy Radon transform

#include <stdlib.h>
#include <math.h>
#include <TNtuple.h>
#include <TVector3.h>
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
#define NGAMMA 40
#define DKAPPA 0.125
#define NKAPPA 25
#define DPHI M_PI/30.
#define NPHI 30
#define TAUMAX 0.5*M_PI
#define DRAWTRACK false

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

TRadon::TRadon(double sig, double thr) : sigma(sig), threshold(thr), Hlist(0)
{
    nt1 = new TNtuple("RadonTransform","Radon Transform","kappa:phi:gamma:sigma:density:x:y:z");
    nt2 = new TNtuple("RadonTrackParameters","Radon Track Parameters","kappa:phi:gamma:sigma:density:x:y:z");

    float gamma = -2.0;
    for (int i=0;i<NGAMMA;i++,gamma+=DGAMMA) {
        string name = "Gamma=" + to_string(gamma);
        string title = "Radon density in r/Phi " + name;
        Hlist.Add(new TH2D(name.data(),title.data(),40,0.0,4.,30,0.0,3.0));
    }
}

TRadon::~TRadon() {
    Hlist.Write();
    nt1->Write();
    delete nt1;
    nt2->Write();
    delete nt2;
}

/*
 * Do a RADON transform of coordinates (all lengths in meter)
 * (1 GeV/c == 2.22m Radius at B=1.5T)
 */
std::vector<RADON>& TRadon::Transform(std::vector<TVector3> &points)
{
    hits = points;
    double kappa,gamma,phi,density,d;
    long k=0,g=0,p=0;
    
    gamma = -2.0;
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
                long i = 0;
                long nhits = 0;
                density=d=0.0;
                vector<TVector3>::iterator it;
                for(it = hits.begin(); it != hits.end(); it++, i++)    {
                    TVector3 point=*it;
                    t.x = point.x();
                    t.y = point.y();
                    t.z = point.z();
                    d =  radon_hit_density(&t);
                    density += d;
                    t.density = density;
                    if (d > 1.0) {
                        // Suppress mirror points, as they yield the same radon value and do not belong to the same track
                        static int oldsx, oldsy, oldsz;
                        int sx = signum(t.x);
                        int sy = signum(t.y);
                        int sz = signum(t.z);
                        if (nhits==0) { oldsx = sx; oldsy = sy; oldsz = sz;}
                        if (sx==oldsx && sy==oldsy && sz==oldsz) {
                            t.index.push_back(i); // Note the indices of the corresponding hits
                            if (nhits==0) printf("\n\nr:%f k:%f p:%f g:%f",1./kappa,kappa,phi,gamma);
                            printf("\nHit #%ld %f,%f,%f| : %f",i,t.x,t.y,t.z,d);
                            nt1->Fill(t.kappa,t.phi,t.gamma,t.sigma,d,t.x,t.y,t.z);
                            oldsx = sx;
                            oldsy = sy;
                            oldsz = sz;
                        }
                        else
                            printf("\nMirror #%ld %f,%f,%f| : %f",i,t.x,t.y,t.z,d);
                        
                        nhits++;
                    }
                }
                if (density > 1.0 && nhits > 3) {
                    printf("\nDensity: %f",density);
                    TH2D *h = (TH2D *) Hlist[(Int_t) g];
                    h->Fill(1./kappa,phi,density);
                    rt.push_back(t);  // Note the candidate track parameters
                }
            }
        }
    }
    
    return rt;
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
    double  eta_g,gamma,kappa,error,z_g,k2,s2,g2,factor, kappafunction,efunction,radon;
    eta_g = getEta_g(t);
    gamma = t->gamma;
    kappa = t->kappa;
    error = t->sigma;
    z_g   = getZ_g(t)  ;
    k2    = kappa * kappa;
    s2    = error * error;
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

void TRadon::GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error) {
    default_random_engine generator;
    double tau = 0.025;
    for (int i=0; i<np; i++,tau+=delta)
    {
        Float_t X,Y,Z;
        X = radius * ( sin(phi + (signum(radius)) * tau) - sin(phi));
        Y = radius * (-cos(phi + (signum(radius)) * tau) + cos(phi));
        Z = gamma * tau;
        if (error > 0.0) {
            normal_distribution<float> distribution0(X,error);
            X = distribution0(generator);
            normal_distribution<float> distribution1(Y,error);
            Y = distribution1(generator);
            normal_distribution<float> distribution2(Z,error);
            Z = distribution2(generator);
        }
        points.push_back(TVector3(X,Y,Z));
    }
}

void TRadon::Draw(Option_t *option) {
    int n = 0;
    vector<RADON>::iterator radon;
    for(radon = rt.begin(); radon != rt.end(); radon++)    {
        
        RADON t=*radon;;
        double radius = 1./t.kappa;
        printf("\nTrack candidate #%d r:%f k:%f p:%f g:%f : %f  ",n++,radius,t.kappa,t.phi,t.gamma,t.density);
        printf("\nAssociated hits: ");
        for (int j = 0; j < t.index.size(); j++) {
            cout << t.index[j] << " ";
        }

        nt2->Fill(t.kappa,t.phi,t.gamma,t.sigma,t.density,t.x,t.y,t.z); // Store results as ROOT ntuple
        
        if (t.density > threshold) {
            printf("Density > %f", threshold);
            std::vector<TVector3> nt3;
            if (DRAWTRACK)
                GenerateTrack(nt3,25,0.025,radius,t.phi,t.gamma);
            else
                for (int i=0;i<t.index.size();i++) {
                    long index = t.index[i];
                    float x = hits[index].x();
                    float y = hits[index].y();
                    float z = hits[index].z();
                    TVector3 p(x,y,z);
                    nt3.push_back(p);
                }
                
            // Draw the track candidates
            TPolyMarker3D *hitmarker = new TPolyMarker3D(nt3.size());
            TPolyLine3D *connector = new TPolyLine3D(nt3.size());
  
            int j = 0;
            vector<TVector3>::iterator it;
            for(it = nt3.begin(); it != nt3.end(); it++, j++)    {
                TVector3 point=*it;
                hitmarker->SetPoint(j,point.x(),point.y(),point.z());
                hitmarker->SetMarkerSize(0.1);
                hitmarker->SetMarkerColor(kYellow);
                hitmarker->SetMarkerStyle(kFullDotLarge);
                hitmarker->Draw(option);
                connector->SetPoint(j, point.x(),point.y(),point.z());
                connector->SetLineWidth(1);
                connector->SetLineColor(kYellow);
                connector->Draw(option);
                printf("\n%d: x:%f y:%f z:%f d:%f  ",j,point.x(),point.y(),point.z(),point.Mag());
            }
        }
    }
}
