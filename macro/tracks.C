#include <random>
#include "Riostream.h"
#include "TVector3.h"

using namespace std;

#define NHITS 20
#define SIGMA 0.001

void GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error);
TVector3 PointWithError(TVector3 p, double sigmax, double sigmay, double sigmaz);

void tracks(long n=1000) {
    Float_t x,y,z;
    Int_t nlines = 0;
    auto f = TFile::Open("tracks.root","RECREATE");
    TH1F h1("h1","x distribution",100,-3,3);
    TH1F h2("h2","y distribution",100,-3,3);
    TH1F h3("h3","z distribution",100,-3,3);
    TH1F h4("h4","r distribution",100,-3,3);
    TH1F h5("h5","phi distribution",100,-4,4);
    TH1F h6("h6","theta distribution",100,-4,4);
    TH1F h7("h7","random r distribution",100,-4,4);
    TH1F h8("h8","random phi distribution",100,-4,4);
    TH1F h9("h9","random theta distribution",100,-4,4);
    TH1F h10("h10","p distribution",100,-10,10);
    TNtuple ntuple("tracks","training data","r1:phi1:theta1:r2:phi2:theta2:truth");
    TNtuple ntuple3("tracks3","training data","r1:phi1:theta1:r2:phi2:theta2:r3:phi3:theta3:truth");
    TRandom r;
    while (nlines++<n) {
        std::vector<TVector3> hits;
        int nhits = NHITS*r.Rndm()+3;
        Double_t p = 0.0;
        while (fabs(p)<0.1) p = 10.*(0.5-r.Rndm()); // 0.1...5 GeV
        h10.Fill(p);
        Double_t phi1 = 2.*(0.5-r.Rndm())*M_PI;
        Double_t gamma1 = 10.*(0.5-r.Rndm());
        Double_t phi2 = 2.*(0.5-r.Rndm())*M_PI;
        Double_t gamma2 = 10.*(0.5-r.Rndm());
        cout << "Generate track with " << nhits << " hits, p = " << p << endl;
        GenerateTrack(hits,nhits,0.0125,1./p,phi1,gamma1,SIGMA);

        TVector3 first(hits[0]);
        TVector3 last(hits[nhits-1]);
        float r0 = first.Mag();
        float r1 = last.Mag();
        float dr = abs(r1-r0);

        // Consider 2 consecutive points on the track
        nhits = hits.size();;
        for (int i=0; i<nhits-1; i++)    {
            TVector3 hit1 = hits[i];
            float r1 = sqrt(hit1.x()*hit1.x()+hit1.y()*hit1.y()+hit1.z()*hit1.z());
            float phi1 = atan2(hit1.y(),hit1.x());
            float theta1 = acos(hit1.z()/r1);
            TVector3 hit2 = hits[i+1];
            float r2 = sqrt(hit2.x()*hit2.x()+hit2.y()*hit2.y()+hit2.z()*hit2.z());
            float phi2 = atan2(hit2.y(),hit2.x());
            float theta2 = acos(hit2.z()/r2);
            ntuple.Fill(r1,phi1,theta1,r2,phi2,theta2,1.0); // true combination
            r2 = r0 + r.Rndm() * dr;// Randomize a point on the sphere
            phi2 = 2.*(0.5-r.Rndm())*M_PI;
            theta2 = r.Rndm()*M_PI;
            ntuple.Fill(r1,phi1,theta1,r2,phi2,theta2,0.0); // wrong combination
            h1.Fill(hit1.x());
            h2.Fill(hit1.y());
            h3.Fill(hit1.z());
            h4.Fill(r1);
            h5.Fill(phi1);
            h6.Fill(theta1);
        }
        
        // Consider 3 consecutive points on the track
        
        for (int i=0; i<nhits-2; i++)    {
            TVector3 hit1 = hits[i];
            float r1 = sqrt(hit1.x()*hit1.x()+hit1.y()*hit1.y()+hit1.z()*hit1.z());
            float phi1 = atan2(hit1.y(),hit1.x());
            float theta1 = acos(hit1.z()/r1);
            TVector3 hit2 = hits[i+1];
            float r2 = sqrt(hit2.x()*hit2.x()+hit2.y()*hit2.y()+hit2.z()*hit2.z());
            float phi2 = atan2(hit2.y(),hit2.x());
            float theta2 = acos(hit2.z()/r2);
            TVector3 hit3 = hits[i+2];
            float r3 = sqrt(hit3.x()*hit3.x()+hit3.y()*hit3.y()+hit3.z()*hit3.z());
            float phi3 = atan2(hit3.y(),hit3.x());
            float theta3 = acos(hit3.z()/r3);
            ntuple3.Fill(r1,phi1,theta1,r2,phi2,theta2,r3,phi3,theta3,1.0); //true combination
            r3 = r0 + r.Rndm() * dr;// Randomize a point on the sphere
            phi3 = 2.*(0.5-r.Rndm())*M_PI;
            theta3 = r.Rndm()*M_PI;
            ntuple3.Fill(r1,phi1,theta1,r2,phi2,theta2,r3,phi3,theta3,0.0); // wrong combination
            h1.Fill(hit1.x());
            h2.Fill(hit1.y());
            h3.Fill(hit1.z());
            h4.Fill(r1);
            h5.Fill(phi1);
            h6.Fill(theta1);
            h7.Fill(r3);
            h8.Fill(phi3);
            h9.Fill(theta3);
        }
        
    }
    
    printf(" generated %d tracks\n",nlines-1);
    //in.close();
    f->Write();
}

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

void GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error) {
    default_random_engine generator;
    double tau = 0.0125;
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

TVector3 PointWithError(TVector3 p, double sigmax, double sigmay, double sigmaz)
{
    default_random_engine generator;

    Float_t X,Y,Z;
    X = p.X();
    Y = p.Y();
    Z = p.Z();
    normal_distribution<float> distribution3(X,sigmax);
    X = distribution3(generator);
    normal_distribution<float> distribution4(Y,sigmay);
    Y = distribution4(generator);
    normal_distribution<float> distribution5(Z,sigmaz);
    Z = distribution5(generator);
    return TVector3(X,Y,Z);
}

