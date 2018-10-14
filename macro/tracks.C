#include <random>
#include "Riostream.h"
#include "TVector3.h"

using namespace std;

#define NHITS 10
#define SIGMA 0.0

void GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error);

void tracks(long n=5) {
    // read file $ROOTSYS/tutorials/tree/basic.dat
    // this file has 3 columns of float data
    //TString dir = gROOT->GetTutorialDir();
    //dir.Append("/tree/");
    //dir.ReplaceAll("/./","/");
    //ifstream in;
    //in.open(Form("%sbasic.dat",dir.Data()));
    Float_t x,y,z;
    Int_t nlines = 0;
    auto f = TFile::Open("tracks.root","RECREATE");
    TH1F h1("h1","x distribution",100,-3,3);
    TH1F h2("h2","y distribution",100,-3,3);
    TH1F h3("h3","z distribution",100,-3,3);
    TH1F h4("h4","r distribution",100,-3,3);
    TH1F h5("h5","phi distribution",100,-3,3);
    TH1F h6("h6","theta distribution",100,-3,3);
    TH1F h7("h7","p distribution",100,-10,10);
    TNtuple ntuple("tracks","training data","r1:phi1:theta1:r2:phi2:theta2:r3:phi3:theta3:truth");
    TRandom r;
    while (nlines++<n) {
        //in >> x >> y >> z;
        //if (!in.good()) break;
        std::vector<TVector3> t1,t2,t3;
        int nhits = NHITS * r.Rndm() + 3; // min. 3 hits
        Double_t p1 = 0.0;
        while (fabs(p1)<0.1) p1 = 10.*(0.5-r.Rndm()); // 0.1...5 GeV
        h7.Fill(p1);
        Double_t phi1 = 2.*(0.5-r.Rndm())*M_PI;
        Double_t gamma1 = 4.*(0.5-r.Rndm());
        cout << "Generate track with " << nhits << " hits, p1 = " << p1 << endl;
        GenerateTrack(t1,nhits,0.025,1./p1,phi1,gamma1,SIGMA);
        Double_t p2 = 0.0;
        while (fabs(p2)<0.1) p2 = 10.*(0.5-r.Rndm()); // 0.1...5 GeV
        h7.Fill(p2);
        Double_t phi2 = 2.*(0.5-r.Rndm())*M_PI;
        Double_t gamma2 = 4.*(0.5-r.Rndm());
        cout << "Generate track with " << nhits << " hits, p2 = " << p2 << endl;
        GenerateTrack(t2,nhits,0.025,1./p2,phi2,gamma2,SIGMA);
        Double_t p3 = 0.0;
        while (fabs(p3)<0.1) p3 = 10.*(0.5-r.Rndm()); // 0.1...5 GeV
        h7.Fill(p3);
        Double_t phi3 = 2.*(0.5-r.Rndm())*M_PI;
        Double_t gamma3 = 4.*(0.5-r.Rndm());
        cout << "Generate track with " << nhits << " hits, p3 = " << p3 << endl;
        GenerateTrack(t3,nhits,0.025,1./p2,phi3,gamma3,SIGMA);

        Double_t truth = 1.0;
        for (int i=0; i<nhits-2; i++)    {
            TVector3 hit1 = t1[i];
            float phi1 = atan2(hit1.y(),hit1.x());
            float r1 = sqrt(hit1.x()*hit1.x()+hit1.y()*hit1.y()+hit1.z()*hit1.z());
            float theta1 = acos(hit1.z()/r1);
            h1.Fill(hit1.x());
            h2.Fill(hit1.y());
            h3.Fill(hit1.z());
            h4.Fill(r1);
            h5.Fill(phi1);
            h6.Fill(theta1);
            for (int j=i+1; j<nhits-1; j++)    {
                TVector3 hit2 = t1[j];
                float phi2 = atan2(hit2.y(),hit2.x());
                float r2 = sqrt(hit2.x()*hit2.x()+hit2.y()*hit2.y()+hit2.z()*hit2.z());
                float theta2 = acos(hit2.z()/r2);
                for (int k=j+1; k<nhits; k++)    {
                    TVector3 hit3 = t1[k];
                    float phi3 = atan2(hit3.y(),hit3.x());
                    float r3 = sqrt(hit3.x()*hit3.x()+hit3.y()*hit3.y()+hit3.z()*hit3.z());
                    float theta3 = acos(hit3.z()/r3);
                    ntuple.Fill(r1,phi1,theta1,r2,phi2,theta2,r2,phi2,theta2,truth);
                }
            }
        }
        
        truth = 0.0;
        for (int i=0; i<nhits-2; i++)    {
            TVector3 hit1 = t1[i];
            float phi1 = atan2(hit1.y(),hit1.x());
            float r1 = sqrt(hit1.x()*hit1.x()+hit1.y()*hit1.y()+hit1.z()*hit1.z());
            float theta1 = acos(hit1.z()/r1);
            h1.Fill(hit1.x());
            h2.Fill(hit1.y());
            h3.Fill(hit1.z());
            h4.Fill(r1);
            h5.Fill(phi1);
            h6.Fill(theta1);
            for (int j=i+1; j<nhits-1; j++)    {
                TVector3 hit2 = t2[j];
                float phi2 = atan2(hit2.y(),hit2.x());
                float r2 = sqrt(hit2.x()*hit2.x()+hit2.y()*hit2.y()+hit2.z()*hit2.z());
                float theta2 = acos(hit2.z()/r2);
                for (int k=j+1; k<nhits; k++)    {
                    TVector3 hit3 = t3[k];
                    float phi3 = atan2(hit3.y(),hit3.x());
                    float r3 = sqrt(hit3.x()*hit3.x()+hit3.y()*hit3.y()+hit3.z()*hit3.z());
                    float theta3 = acos(hit3.z()/r3);
                    ntuple.Fill(r1,phi1,theta1,r2,phi2,theta2,r2,phi2,theta2,truth);
                }
            }
        }
    }
    
    printf(" generated %d tracks\n",nlines-1);
    //in.close();
    f->Write();
}

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

void GenerateTrack(std::vector<TVector3> &points, int np, double delta, double radius, double phi, double gamma, double error) {
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

