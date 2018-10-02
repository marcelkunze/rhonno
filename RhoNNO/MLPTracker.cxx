// Read HITS spatial data and try a XMLP network to segment the tracks...

#include <TROOT.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TAxis3D.h>
#include <TPolyLine3D.h>
#include <TFile.h>
#include <TVector3.h>
#include "RhoNNO/TXMLP.h"

#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

#define NHITS 5
#define SIGMA 0.001

#define TRACKLET 3
#define DISTANCE 100
#define THRESHOLD 65

// The user member function processes one event

std::vector<TVector3> hits;
std::vector<int> tracks[20000];

#define signum(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

Double_t* Recall(Double_t *invec);
int bestMatchingHit(size_t nhits, int **m, int row);
int bestHitPair(size_t nhits, int **m, int &row, int &col);
int findTracks(int nhits, float *x, float *y, float *z, int* labels);

bool sortFunc( const vector<int>& p1,
              const vector<int>& p2 ) {
    return p1.size() > p2.size();
}

void print(vector<int> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
    cout << endl;
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

int getFarestHit(vector<int>, int* d);

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
    
    cout << "Sorting hits..." << endl;
    reverse(hits.begin(),hits.end());
    
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
    for (int i=0; i<nhits; i++) {
        TVector3 hit = hits[i];
        x[i] = hit.X();
        y[i] = hit.Y();
        z[i] = hit.Z();
    }
    nt = findTracks(nhits,x,y,z,labels);
    cout << "Number of tracks:" << nt << endl;
    for(int i=0; i<nt; i++) {
        cout << "Track " << i << ":";
        print(tracks[i]);
    }
    
    
    cout << "Labels: ";
    for (int i=0;i<nhits;i++) cout << labels[i] << " ";
    
    //TBD: Fit the helix tracks from the hits in the containers
    
    cout << endl << "Writing..." << endl;
    c1->Write();
    
    return EXIT_SUCCESS;
}

// Assign track labels to hits (x,y,z)
// The hit pair quality is assessed by the neural network
// The quality is noted in the hit pair matrix m[nhits][nhits]
int findTracks(int nhits, float *x, float *y, float *z, int* labels)
{
    for (int i=0;i<nhits;i++) labels[i] = -1; // Preset with no match
    
    Double_t in1[7], in2[7];
    
    // Allocate a nhits*nhits hit pair matrix as one continuous memory block
    int **m = new int*[nhits];
    if (nhits)
    {
        m[0] = new int[nhits * nhits];
        for (int i = 1; i < nhits; ++i)
            m[i] = m[0] + i * nhits;
    }
    
    // Allocate a nhits*nhits hit pair distance matrix as one continuous memory block
    int **d = new int*[nhits];
    if (nhits)
    {
        d[0] = new int[nhits * nhits];
        for (int i = 1; i < nhits; ++i)
            d[i] = d[0] + i * nhits;
    }
    
    for(int i=0; i<nhits-1; i++)    {
        in1[0] = x[i];
        in1[1] = y[i];
        in1[2] = z[i];
        in2[3] = x[i];
        in2[4] = y[i];
        in2[5] = z[i];
        d[i][i] = 0;
        m[i][i] = 0;
        // Polar coordinates
        //double phi = atan2(y[i],x[i]);
        //double r = sqrt(x[i]*x[i]+y[i]*y[i]);
        //double zz = z[i];
        for(int j=i+1; j<nhits; j++)    {
            in1[3] = x[j];
            in1[4] = y[j];
            in1[5] = z[j];
            in2[0] = x[j];
            in2[1] = y[j];
            in2[2] = z[j];
            double dist = sqrt((in1[0]-in2[0])*(in1[0]-in2[0]) + (in1[1]-in2[1])*(in1[1]-in2[1]) + (in1[2]-in2[2])*(in1[2]-in2[2]));
            in1[6] = dist;
            in2[6] = dist;
            d[i][j] = 1000.*dist; // Cache the distance between the points
            d[j][i] = 1000.*dist;
            d[j][j] = 0;
            m[j][j] = 0;
            m[i][j] = (int) 100. * Recall(in1)[0]; // Recall the hit pair matching quality
            m[j][i] = (int) 100. * Recall(in2)[0];
            if (m[i][j]< THRESHOLD) m[i][j] = 0; // Apply a cut on the quality
            if (m[j][i]< THRESHOLD) m[j][i] = 0;
        }
    }
    
    // Search neighbouring hits, the neural network recall identifies the hit belonging to a tracklet
    vector<vector<int>> tracklet;
    for(int i=0; i<nhits-1; i++)    {
        vector<int> tmpvec;
        tmpvec.push_back(i); //Note the row in the first place
        for(int j=i; j<nhits; j++)    {
            int dist = d[i][j];
            int recall  = (m[i][j]>m[j][i]) ? m[i][j]:m[j][i];
            if (dist < DISTANCE && recall>THRESHOLD) {
                tmpvec.push_back(j); // Note the columns with a good combination
            }
        }
        
        tracklet.push_back(tmpvec);
    }
    
    cout << "Number of tracklets: " << tracklet.size() << endl;
    
    // Sort the tracklet vector according to the tracklet length
    
    sort(tracklet.begin(), tracklet.end(), sortFunc);
    
    // Print out the sorted vector
    cout << "Sorted tracklets:" << endl;
    for( int i=0; i<tracklet.size(); i++ ) {
        for( int j=0; j<tracklet[i].size(); j++ ) {
            int row = tracklet[i][0];
            int col = tracklet[i][j];
            cout << tracklet[i][j] << "(" << m[row][col] << ") ";
        }
        cout << endl;
    }
    
    cout << "Seed: " << tracklet[0][0] << " length: " << tracklet[0].size() << endl;
    
    // Prune the tracklets by removing short tracks
    for (vector<vector<int>>::iterator it = tracklet.begin(); it != tracklet.end(); ++it) {
        vector<int> row = *it;
        if (row.size() < TRACKLET) { // Remove short tracklets
            tracklet.erase(it);
            *it--;
            continue;
        }
    }
    
    // Print out the pruned vector
    cout << "Pruned tracklets (remove short paths):" << endl;
    for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    
    // Prune the tracklets by removing rows with identical entries
    // Assemble tracks from the corresponding tracklets
    vector<vector<int>> track;
    for (vector<vector<int>>::iterator it = tracklet.begin(); it != tracklet.end(); ++it) {
        vector<int> row = *it;
        vector<int> tmpvec = row; // Vector to assemble the track
        for (vector<int>::iterator it2 = row.begin(); it2 != row.end(); ++it2) { // Search next seeding hits in remainder tracklet list
            int prune = *it2;
            for (vector<vector<int>>::iterator it3 = it+1; it3 != tracklet.end(); ++it3) {
                vector<int> nextrow = *it3;
                bool hitExists = find(nextrow.begin(),nextrow.end(),prune) != nextrow.end();
                if (hitExists) {
                    for (int j=0;j<nextrow.size();j++) tmpvec.push_back(nextrow[j]); // Append the hits to track before erasing the row
                    tracklet.erase(it3);
                    *it3--;
                    continue;
                }
            }
        }
        set<int> s( tmpvec.begin(), tmpvec.end() ); // Remove duplicates
        tmpvec.assign( s.begin(), s.end() );
        track.push_back(tmpvec);
    }
    
    // Print out the pruned vector
    cout << "Pruned tracklets:" << endl;
    for( int i=0; i<tracklet.size(); i++ ) print(tracklet[i]);
    
    // Print out the tracks vector
    cout << "Tracks:" << endl;
    for( int i=0; i<track.size(); i++ ) print(track[i]);
    
    if (nhits) delete [] m[0]; // Clean up the memory
    delete [] m;
    
    if (nhits) delete [] d[0]; // Clean up the memory
    delete [] d;
    
    for (int i=0;i<track.size();i++) {
        tracks[i] = track[i]; // Save the results
        for (int j=0;j<track[i].size();j++) {
            int hit = track[i][j];
            labels[hit] = i;
        }
    }
    return (int) track.size();
}

// Recall function on normalised network input
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
    x[6]     = 126.061    *    invec[6];    // dot
    return net.Recall(x,y);
}

