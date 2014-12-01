// This example illustrates unsupervised training by use of a Growing Cell Structure (GCS).
// The input data file comprises 2D-vectors describing kinematics of three-body particle decays 
// (Dalitz-plot). It is the task of the GCS network to fit the input density distribution.
// (C) Marcel Kunze, Experimentalphysik 1, Bochum University
//
// For further information see
// http://www.ep1.ruhr-uni-bochum.de/~marcel/tutorial.html

#include "TROOT.h"
#include "RhoNNO/TGCS.h"

#include <iostream>
using namespace std;

TROOT root("root","Unsupervised training test");

NNO_INTYPE In[2];

#define MAXCELL 1000

int main(void) { 
    FILE* F=fopen("ppe.dat","r");
    
    TGCS GCS( 2,       // 2 input nodes 
	3,       // start with 3 cells 
	MAXCELL,   // stop on 100 cells 
	0.2,       // Learnstep of Winner-Cell 
	0.02,      // Learnstep of Neighbours 
	0.002,     // a_win_count 
	10,        // connectors
	100,       // insert_step
	500,       // delet_step 
	"gcs.net"  // Network Filename is "gcs.net" 
	); 
    
    int EpC=0;
    while (GCS.GetNumberOfCells()<MAXCELL) { 
	while (!feof(F)) {
	    int ip1,ip2;
	    fscanf(F,"%d %d",&ip1,&ip2);
	    In[0] = ip1;
	    In[1] = ip2;
	    GCS.Learnstep(In);             //perform a learnstep 
	} 
	rewind(F); //new epoch 
	printf("epoch nr.: %i, cells: %i\n",EpC++,GCS.GetNumberOfCells());
    } 
    fclose(F);
    return 0;
} 
