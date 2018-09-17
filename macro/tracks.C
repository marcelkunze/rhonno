#include <random>
#include "Riostream.h"
#include "TVector3.h"

using namespace std;

#define NHITS 10
#define SIGMA 0.001

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
    TH1F h1("h1","x distribution",100,-1,1);
    TH1F h2("h2","y distribution",100,-1,1);
    TH1F h3("h3","z distribution",100,-1,1);
    TH1F h4("h4","p distribution",100,-10,10);
    TNtuple ntuple("tracks","training data","x1:y1:z1:x2:y2:z2:cos:c1:c2:truth:p");
    TRandom r;
    while (nlines<n) {
        //in >> x >> y >> z;
        //if (!in.good()) break;
        std::vector<TVector3> t1,t2;
        int nhits = NHITS * r.Rndm() + 3; // min. 3 hits
        Double_t p1 = 0.0;
        while (fabs(p1)<0.1) p1 = 10.*(0.5-r.Rndm()); // 0.1...5 GeV
        h4.Fill(p1);
        Double_t phi1 = 2.*r.Rndm()*M_PI;
        Double_t gamma1 = 2.*(0.5-r.Rndm());
        cout << "Generate track with " << nhits << " hits, p1 = " << p1 << endl;
        GenerateTrack(t1,nhits,0.025,1./p1,phi1,gamma1,SIGMA);
        Double_t p2 = 0.0;
        while (fabs(p2)<0.1) p2 = 10.*(0.5-r.Rndm()); // 0.1...5 GeV
        h4.Fill(p2);
        Double_t phi2 = 2.*r.Rndm()*M_PI;
        Double_t gamma2 = 2.*(0.5-r.Rndm());
        cout << "Generate track with " << nhits << " hits, p2 = " << p2 << endl;
        GenerateTrack(t2,nhits,0.025,1./p2,phi2,gamma2,SIGMA);
        
        Double_t truth = 1.0;
        for(int i=0; i<nhits-1; i++)    {
           TVector3 hit1 = t1[i];
           Double_t c1 = hit1.CosTheta();
           for(int j=i+1; j<nhits; j++)    {
                TVector3 hit2 = t1[j];
                Double_t c2 = hit2.CosTheta();
                Double_t cos = hit1.Dot(hit2);
                if (nlines < 5) printf("x1=%8f, y1=%8f, z1=%8f x2=%8f, y2=%8f, z2=%8f t=%8f\n",hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),truth);
                h1.Fill(hit1.x());
                h2.Fill(hit1.y());
                h3.Fill(hit1.z());
                ntuple.Fill(hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),cos,c1,c2,truth,p1);
                ntuple.Fill(hit2.x(),hit2.y(),hit2.z(),hit1.x(),hit1.y(),hit1.z(),cos,c2,c1,truth,p1);
            }
        }
        for(int i=0; i<nhits-1; i++)    {
            TVector3 hit1 = t2[i];
            Double_t c1 = hit1.CosTheta();
            for(int j=i+1; j<nhits; j++)    {
                TVector3 hit2 = t2[j];
                Double_t c2 = hit2.CosTheta();
                Double_t cos = hit1.Dot(hit2);
                if (nlines < 5) printf("x1=%8f, y1=%8f, z1=%8f x2=%8f, y2=%8f, z2=%8f t=%8f\n",hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),truth);
                h1.Fill(hit1.x());
                h2.Fill(hit1.y());
                h3.Fill(hit1.z());
                ntuple.Fill(hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),cos,c1,c2,truth,p1);
                ntuple.Fill(hit2.x(),hit2.y(),hit2.z(),hit1.x(),hit1.y(),hit1.z(),cos,c2,c1,truth,p1);
            }
        }

        truth = 0.0;
        for(int i=0; i<nhits-1; i++)    {
            TVector3 hit1 = t1[i];
            Double_t c1 = hit1.CosTheta();
            for(int j=i+1; j<nhits; j++)    {
                TVector3 hit2 = t2[j];
                Double_t c2 = hit2.CosTheta();
                Double_t cos = hit1.Dot(hit2);
                if (nlines < 5) printf("x1=%8f, y1=%8f, z1=%8f x2=%8f, y2=%8f, z2=%8f t=%8f\n",hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),truth);
                h1.Fill(hit2.x());
                h2.Fill(hit2.y());
                h3.Fill(hit2.z());
                ntuple.Fill(hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),cos,c1,c2,truth,p1);
                ntuple.Fill(hit2.x(),hit2.y(),hit2.z(),hit1.x(),hit1.y(),hit1.z(),cos,c2,c1,truth,p1);
            }
        }
        for(int i=0; i<nhits-1; i++)    {
            TVector3 hit1 = t2[i];
            Double_t c1 = hit1.CosTheta();
            for(int j=i+1; j<nhits; j++)    {
                TVector3 hit2 = t1[j];
                Double_t c2 = hit2.CosTheta();
                Double_t cos = hit1.Dot(hit2);
                if (nlines < 5) printf("x1=%8f, y1=%8f, z1=%8f x2=%8f, y2=%8f, z2=%8f t=%8f\n",hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),truth);
                h1.Fill(hit2.x());
                h2.Fill(hit2.y());
                h3.Fill(hit2.z());
                ntuple.Fill(hit1.x(),hit1.y(),hit1.z(),hit2.x(),hit2.y(),hit2.z(),cos,c1,c2,truth,p1);
                ntuple.Fill(hit2.x(),hit2.y(),hit2.z(),hit1.x(),hit1.y(),hit1.z(),cos,c2,c1,truth,p1);
            }
        }

        nlines++;
    }
    printf(" generated %d tracks\n",2*nlines);
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

