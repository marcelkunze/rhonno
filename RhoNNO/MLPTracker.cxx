// Read HITS spatial data and try a XMLP network to segment the tracks...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include "RhoNNO/TXMLP.h"

#include <random>
#include <iostream>
#include <fstream>
using namespace std;

#define NHITS 5
#define SIGMA 0.001

#define THRESHOLD 0.5

// The user member function processes one event

std::vector<TVector3> hits;
std::vector<int> tracks[100];

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

Double_t* Recall(Double_t *invec);
int bestMatchingHit(TMatrixD &m, int row);
int bestHitPair(TMatrixD &m, int &row, int &col);
int findTracks(int nhits, float *x, float *y, float *z, int* labels);

void print(vector<int> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
}

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

int main(int argc, char* argv[]) {
    
    TFile output("MLPTracker.root","RECREATE");
    
    string filename("event");
    TNtuple nt1("Hits","NNO Tracking Data","x:y:z");
    
    if (argc > 1) filename = argv[1];
    
    ifstream infile(filename);
    if (infile) {
        cout << "Reading input file: " << filename << endl;
        double X,Y,Z;
        while (infile >> X >> Y >> Z) {
            X*=0.01; // transform to meter
            Y*=0.01;
            Z*=0.01;
            TVector3 point(X,Y,Z);
            hits.push_back(point);
            //cout << point.x() << "\t" << point.y() << "\t" << point.z() << "\t"<< point.d() << endl;
        }
    }
    else
    {
        // std::vector<TVector3>, int np, float delta tau, float radius, float phi, float gamma
        GenerateTrack(hits,NHITS,0.0125,1.0,M_PI/1.0,0.5,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/3.0,1.5,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/4.0,-1.2,SIGMA);
        GenerateTrack(hits,NHITS,0.0125,-1.0,M_PI/2.0,1.0,SIGMA);
    }
    
    // Sort the hits according to distance from origin
    
    //cout << "Sorting hits..." << endl;
    //reverse(hits.begin(),hits.end());
    
    // Initialize a 3D canvas and draw the hits
    TCanvas *c1 = new TCanvas("c1","NNO Tracking: XMLP",200,10,700,500);
    // create a pad
    TPad *p1 = new TPad("p1","p1",0.05,0.02,0.95,0.82,46,3,1);
    p1->SetFillColor(kBlack);
    p1->Draw();
    p1->cd();
    // creating a view
    TView *view = TView::CreateView(1);
    view->SetRange(-2,-2,-2,2,2,2); // draw in a 2 meter cube
    // Draw axis
    TAxis3D rulers;
    rulers.Draw();
    // draw hits as PolyMarker3D
    const int nhits = (int) hits.size();
    TPolyMarker3D *hitmarker = new TPolyMarker3D((UInt_t) nhits);
    vector<TVector3>::iterator it;
    for(it = hits.begin(); it != hits.end(); it++)    {
        static int i = 0;
        TVector3 p=*it;
        hitmarker->SetPoint(i++,p.x(),p.y(),p.z());
        nt1.Fill(p.x(),p.y(),p.z());
    }
    // set marker size, color & style
    hitmarker->SetMarkerSize(1.0);
    hitmarker->SetMarkerColor(kCyan);
    hitmarker->SetMarkerStyle(kStar);
    hitmarker->Draw();
    
    cout << endl << "Reconstruct Track data with a XMLP Network";
    cout << endl << " Number of hits:" << nhits << endl;
    
    float x[nhits],y[nhits],z[nhits];
    int labels[nhits];
    int nt;
    nt = findTracks(nhits,x,y,z,labels);
    
    cout << "Labels: ";
    for (int i=0;i<nhits;i++) cout << labels[i] << " ";
    
    //TBD: Fit the helix tracks from the hits in the containers
    
    cout << endl << "Writing..." << endl;
    c1->Write();
    
    return EXIT_SUCCESS;
}

int findTracks(int nhits, float *x, float *y, float *z, int* labels)
{
    for (int i=0;i<nhits;i++) labels[i] = -1; // Preset with no match
    
    Double_t in1[7], in2[7], *out;
    TMatrixD m(nhits,nhits);
    for(int i=0; i<nhits-1; i++)    {
        TVector3 hit1 = hits[i];
        Double_t c1 = hit1.CosTheta();
        in1[0] = hit1.x();
        in1[1] = hit1.y();
        in1[2] = hit1.z();
        in2[3] = hit1.x();
        in2[4] = hit1.y();
        in2[5] = hit1.z();
        for(int j=i+1; j<nhits; j++)    {
            TVector3 hit2 = hits[j];
            Double_t c2 = hit2.CosTheta();
            Double_t cos = hit1.Dot(hit2);
            in1[3] = hit2.x();
            in1[4] = hit2.y();
            in1[5] = hit2.z();
            in2[0] = hit2.x();
            in2[1] = hit2.y();
            in2[2] = hit2.z();
            in1[6] = cos;
            in2[6] = cos;
            m[i][j] = Recall(in1)[0];
            m[j][i] = Recall(in2)[0];
            if (m[i][j]< THRESHOLD) m[i][j] = 0.0;
            if (m[j][i]< THRESHOLD) m[j][i] = 0.0;
        }
    }
    m.Print();
    
    // Analyze the network
    // Sort out the tracks by following the network connections and fill the corresponding track hits into containers
    int ntracks = 0;
    
    int row, col;
    int seed = bestHitPair(m, row, col); // Look for seed
    while (seed > -1) {
        cout << "Best hit pair: (" << row << "," << col << ")" << endl;
        while (seed>-1) {
            labels[seed] = ntracks;
            tracks[ntracks].push_back(seed);
            seed = bestMatchingHit(m, seed);
            m.Print();
        }
        ntracks++;
        seed = bestHitPair(m, row, col); // Look for new seed
    }
    
    for(int i=0; i<ntracks; i++) {
        cout << "Track " << i << ":";
        print(tracks[i]);
        cout << endl;
    }
    
    return ntracks;
}

int bestMatchingHit(TMatrixD &m, int row)
{
    long nhits = m.GetNcols();
    int imax = -1;
    static int imaxold = -1;
    
    cout << "Seed hit: " << row << endl;
    Double_t max = 0.0;
    for(int i=0; i<nhits; i++)    {
        if (i!=row && m[row][i] > max) {
            max = m[row][i]; // Best matching hit
            imax = i;
        }
    }
    cout << "  Best matching hit: " << imax << endl;
    for(int i=0; i<nhits; i++) m[row][i] = -1.0;
    for(int i=0; i<nhits; i++) m[i][row] = -1.0;
    if (imax == imaxold) {
        //for(int i=0; i<nhits; i++) m[i][imax] = -1.0;
        //for(int i=0; i<nhits; i++) m[imax][i] = -1.0;
        return -1;
    }
    imaxold = imax;
    return imax;
}

int bestHitPair(TMatrixD &m, int &row, int &col)
{
    long nhits = m.GetNcols();
    col = -1;
    row = -1;
    
    Double_t max = -1.0;
    for(int i=0; i<nhits; i++)    {
        for(int j=0; j<nhits; j++)    {
            if (m[i][j] > max) {
                max = m[i][j]; // Best hit pair
                row = i;
                col = j;
            }
        }
    }
    return row;
}

Double_t* Recall(Double_t *invec)
{
    static TXMLP net("/Users/marcel/workspace/rhonno/RhoNNO/NNO0069.TXMLP");
    Float_t x[7],y[1];
    x[0]     = 2631.15    *    invec[0];    // x1
    x[1]     = 383.788    *    invec[1];    // y1
    x[2]     = 490.839    *    invec[2];    // z1
    x[3]     = 2631.15    *    invec[3];    // x2
    x[4]     = 383.788    *    invec[4];    // y2
    x[5]     = 490.839    *    invec[5];    // z2
    x[6]     = 126.061    *    invec[6];    // cos
    return net.Recall(x,y);
}

