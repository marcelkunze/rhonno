// This example illustrates unsupervised training by use of a Growing Cell Structure (GCS).
// The input data file comprises 2D-vectors describing kinematics of three-body particle decays 
// (Dalitz-plot). It is the task of the GCS network to fit the input density distribution.
// (C) Marcel Kunze, Experimentalphysik 1, Bochum University
//
// For further information see
// http://www.ep1.ruhr-uni-bochum.de/~marcel/tutorial.html

#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TGraph2D.h>
#include <TGraphDelaunay2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include "RhoNNO/TGCS.h"

#include <iostream>
using namespace std;

TROOT root("root","Unsupervised training test");

#define MAXCELL 100
#define MAXEPOCH 200

int main(int argc, char* argv[]) {
    string filename("ppe.dat");
    if (argc > 1) filename = argv[1];
    
    cout << "Reading input file: " << filename << endl;
    FILE* F=fopen(filename.data(),"r");
    if (!F) {
        cerr << "File does not exist!" << endl;
        return EXIT_FAILURE;
    }
    
    TFile *f = new TFile("utrain.root","RECREATE");
    
    TNtuple tuple("Ntuple","PiPiEta Dalitzplot","m1:m2");
    while (!feof(F)) {
        float m1,m2;
        fscanf(F,"%f %f",&m1,&m2);
        tuple.Fill(m1,m2);
    }
    fclose(F);
    
    long nentries = tuple.GetEntries();
    
    cout << endl << "Training PiPiEta Dalitzplot with a GCS Network";
    cout << endl << "Number of points:" << nentries << endl;
    
    TGCS net( 2,       // 2 input nodes
             3,       // start with 3 cells
             MAXCELL,   // stop on 100 cells
             0.2,       // Learnstep of Winner-Cell
             0.02,      // Learnstep of Neighbours
             0.002,     // a_win_count
             10,        // connectors
             100,       // insert_step
             1000,      // delete_step
             "gcs.net"  // Network Filename is "gcs.net"
             );
    
    
    int EpC=0;
    while (EpC++<MAXEPOCH) {
        net.TrainEpoch(&tuple);
        printf("epoch nr.: %i, cells: %i\n",EpC,net.GetNumberOfCells());
    }
    
    // Fill a Graph
    unsigned int numberCells = net.GetNumberOfCells();
    TCanvas *c = new TCanvas("c","Graph2D example",-8000,-8000,8000,8000);
    c->Divide(2,1);
    TGraph2D *graph2d = new TGraph2D(numberCells);
    graph2d->SetName("Graph2D");
    graph2d->SetTitle("GCS: Voronoi Diagram");
    
    for (unsigned long i=0;i<numberCells;++i) {
        const TNeuralNetCell *cell = net.GetCell(i);
        const double *x = cell->GetVector();
        printf("\n Cell %ld: (%f,%f) \n",i,x[0],x[1]);
        graph2d->SetPoint(i,x[0],x[1],0.0);
    }
    
    c->cd(1);
    tuple.SetMarkerSize(0.4);
    tuple.SetMarkerStyle(20);
    tuple.Draw("m1:m2");
    c->Update();
    
    c->cd(2);
    net.Draw();
    c->Update();
    
    TGraphDelaunay2D *delaunay = new TGraphDelaunay2D(graph2d);
    delaunay->FindAllTriangles();
    delaunay->Write();
    
    c->Write();
    f->Write();
    
    return EXIT_SUCCESS;
}
