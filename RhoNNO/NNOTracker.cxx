// Read HITS spatial data and try a GNG network to segment the tracks...

#include "TROOT.h"
#include "TVector3.h"
#include "RhoNNO/TGNG.h"

#include <iostream>
using namespace std;

#define MAXEPOCH 100
#define DELETESTEP 500

// The user member function processes one event

#define     MAXHITS 1000
TVector3	hits[MAXHITS];
UInt_t      nhits=0, n=0;

int main(int argc, char* argv[]) {
    TString filename("event");
    if (argc > 1) filename = argv[1];
    
    cout << "Reading input file: " << filename << endl;
    FILE* F=fopen(filename,"r");
    if (!F) {
        cerr << "File does not exist!" << endl;
        return EXIT_FAILURE;
    }

	UInt_t i,j;
	double x,y,z;

    while(!feof(F)) {
        fscanf(F,"%lf %lf %lf",&x,&y,&z);
        if (nhits<MAXHITS) {
            hits[nhits].SetX(x);
            hits[nhits].SetY(y);
            hits[nhits++].SetZ(z);
        }
        else
            return EXIT_FAILURE;
    }
    fclose(F);
    
    //TBD: Initialize a 3D canvas and draw the hits
    // ...

    cout << endl << "Training Delphi TPC data with a GNG Network";
    cout << endl << " Number of hits:" << nhits << endl;
/*
// Funktioniert sehr gut:
gng net( 3,       // 3 input nodes 
        nhits,    // one cell per hit 
       0.5,       // Learnstep of Winner-Cell 
       0.05,      // Learnstep of Neighbours 
       0.8,       // a_win_count 
       0.01,      // a_edge_count 
       0.1,       // min_count
       3,         // connectors 
       10,        // insert_step
       500,       // delete_step 
       "gng.net"  // Network Filename is "gng.net" 
       ); 
*/

    TGNG net( 3,       // 3 input nodes
        nhits,    // 1 cell per hit
        0.5,       // Learnstep of Winner-Cell
        0.05,      // Learnstep of Neighbours
        0.8,       // a_win_count
        0.01,      // a_edge_count
        0.5,       // min_count
        3,         // connectors
        10,       // insert_step
        DELETESTEP,       // delete_step
       "gng.net"  // Network Filename is "gng.net" 
    );

	while (n++ < MAXEPOCH) {
            // Perform training with all hits
            for (i=0;i<nhits;i++) {
                Float_t x[3];
                hits[random()%nhits].GetXYZ(x);
                net.Learnstep(x);
            }

            // Show network evolution
            UInt_t numberCells = net.GetNumberOfCells();
            cout << endl << endl << "Epoch: " << n << endl << "Cells:" << numberCells;
            for (i=0;i<numberCells;++i) {
                const TNeuralNetCell *c = net.GetCell(i);
                const Double_t *x1 = c->GetVector();
                printf("\n Cell %d: X1(%f,%f,%f) \n",i,x1[0],x1[1],x1[2]);
                //TBD: Draw the cell location X1
                // ...
                UInt_t numberConnections = c->GetNumberOfConnections();
                for (j=0;j<numberConnections;j++) {
                    const TNeuralNetCell *u = c->GetConnectedCell(j);
                    const Double_t *x2 = u->GetVector();
                    const Int_t id = u->GetID();
                    printf("\t -> %d. X2(%f,%f,%f) \n",id,x2[0],x2[1],x2[2]);
                    //TBD: Draw the cell connections (X1->X2)
                    // ...
                }
            }
	}
    
    //TBD: Analyze the network
    // Sort out the tracks by following the network connections and fill the corresponding track hits into containers
    
    //TBD: Fit the helix tracks from the hits in the containers

	return EXIT_SUCCESS;
}

