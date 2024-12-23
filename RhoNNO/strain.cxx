// Example for supervised training
//
// This sample file comprises supervised network training and recall
// examples for TMLP, TXMLP, TSGCS and TSGNG of the Neural Network
// Objects package (NNO).
// The program produces a simple training file in order to run the
// different neural networks contained in the nno-classlib.
// In this simple example the networks have to learn how to compute the
// parity  of a binary number. In supervised training the number occurs
// as inputvector with INNODES double precision values of 1 or 0.
// The training output contains the parity of the inputvector:
// 1 if even, -1 if odd.
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995-2001

const int INNODES = 5;
const int HIDNODES = 10;
const int CELLS = 100;

#define PLOTS false

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TList.h"
#include "TIterator.h"

#include "TNNK.h"
#include "TMLP.h"
#include "TXMLP.h"
#include "TSGCS.h"
#include "TSGNG.h"

#include <iostream>
using namespace std;

TROOT root("root","Supervised training test");
TApplication app("NNO Test program",0,0);


// Class to train a list of networks with parity data

class NetworkTrainer
{
public:
    enum PARITY {EVEN,ODD};
    
    NetworkTrainer(int nEpoch=100);
    virtual ~NetworkTrainer() {}
    void AddNetwork(VNeuralNet *net) { fNets.Add(net); }
    void Train();
    int MakeData(NNO_INTYPE*, bool randomize = false);
    void MakeFile(const char*,int, bool randomize = false);
    PARITY Parity(NNO_INTYPE*);
private:
    int fEpoch;
    TList fNets;
};


// Program to train 5 supervised networks: TNNK, MLP, XMLP, TGCS and TGNG
// (Parity problem)

int main() {
    
    cout << "NNO TEST PROGRAM" << endl;
    
    NetworkTrainer trainer(50); // Number of training epochs
    
    // Simple Multi-Layer-Perceptron.
    // Note that the output-layer-nodes of MLP contain no sigmoid function.
    // So the output range of MLP is not limited.
    
    VNeuralNet *net1 = new TMLP(
                                0.1,         // learningstep of neural weights of hidden neurons
                                0.01,        // learningstep of neural weights of output neurons
                                INNODES,     // number of inputnodes
                                HIDNODES,    // number of hiddennodes
                                1,           // number of outputnodes
                                1.0,         // average value range of inputdata
                                "mlp.net");  // network - filename (used by destructor)
    
    net1->SetMomentumTerm(0.4);
    net1->SetThreshold(0.0);
    trainer.AddNetwork(net1);
    
    // Multi-Layer-Perceptron training using xmlp.
    
    VNeuralNet *net2 = new TXMLP(
                                 3,		    // Layers  (2 hidden layer)
                                 1.0,		    // InputRange
                                 "xmlp.net",	    // NetFile
                                 INNODES,	    // inputnodes
                                 HIDNODES, 2   , 1,  // nodes ->  x-10-2-1 network
                                 0.1, 0.02, 0.01,    // learningsteps
                                 TNeuralNetParameters::TR_LINEAR_BEND,
                                 TNeuralNetParameters::TR_LINEAR_BEND,
                                 TNeuralNetParameters::TR_LINEAR);   // transferfunctions
    
    net2->SetMomentumTerm(0.0001);
    net2->SetThreshold(0.0);
    trainer.AddNetwork(net2);
    
    // Simple Supervised-Growing-Cell-Structure
    // Note that the output nodes of sgcs contain no sigmoid function and
    // no activity-threshold. Their output value range is not limited.
    // If the input vector does not match well to any cell, the output is
    // always near zero.
    
    VNeuralNet *net3 = new TSGCS(
                                 INNODES,   // number of inputnodes
                                 4,         // starting number of cells (dimension of connection topologie is Cells-1)
                                 1,         // number of outputnodes
                                 CELLS,     // maximum number of cells
                                 0.1,       // learningstep of winner cell
                                 0.002,     // learningstep of neighbour
                                 0.1,       // learningstep of neural weights
                                 0.01,      // decrement of error_count: err_count*=(1.0 - AErrCount)
                                 0.01,      // adjustment of squared deviation during learnig
                                 20,        // maximum number of allowed connections for one cell
                                 1000,      // cell insertion after 1000 Learningsteps
                                 0,         // don't remove cells
                                 "sgcs.net");  // network - filename (used by destructor)
    
    net3->SetThreshold(0.0);
    trainer.AddNetwork(net3);
    
    
    // Simple Supervised-Growing-Neural-Gas training
    // Note that the output nodes of sgng contain no sigmoid function and
    // no activity-threshold. Their output value range is not limited.
    // If the input vector does not match well to any cell, the output is
    // always near zero.
    
    VNeuralNet *net4 = new TSGNG(
                                 INNODES,      // number of inputnodes
                                 1,            // number of outputnodes
                                 CELLS,        // maximum number of cells
                                 0.1,          // learningstep of winner cell
                                 0.002,        // learningstep of neighbour
                                 0.1,          // learningstep of neural weights
                                 0.01,         // decrement of error_count: err_count*=(1.0 - AErrCount)
                                 0.01,         // decrement of edge_count: edge_count*=(1.0 - AEdgeCount)
                                 0.01,         // adjustment of squared deviation during learnig
                                 0.01,         // threshold for edge - removal
                                 20,           // maximum number of allowed connections for one cell
                                 1000,         // cell insertion afer 1000 Learningsteps
                                 100,          // attempt of edge removal afer 100 Learningsteps
                                 "sgng.net");  // network - filename (used by destructor)
    
    net4->SetThreshold(0.0);
    trainer.AddNetwork(net4);
    
    
    // Perform the network training
    
    cout << "Training phase..." << endl;
    trainer.Train();
    
    
    // Check network performance:
    // Read in the files and perform a recall
    
    cout << endl << "Test phase..." << endl;
    
    
    TMLP  Net1("mlp.net");
    cout << Net1.GetNetID() << ": ClassError = "
    << Net1.TestEpoch("parity.tst") << endl;
    
    TXMLP Net2("xmlp.net");
    cout << Net2.GetNetID() << ": ClassError = "
    << Net2.TestEpoch("parity.tst") << endl;
    
    TSGCS Net3("sgcs.net");
    cout << Net3.GetNetID() << ": ClassError = "
    << Net3.TestEpoch("parity.tst") << endl;
    
    TSGNG Net4("sgng.net");
    cout << Net4.GetNetID() << ": ClassError = "
    << Net4.TestEpoch("parity.tst") << endl;
    
    
    VNeuralNet *testNet = &Net4;
    
    cout << endl << "Classify a few samples with " << testNet->GetNetID() << endl;
    
    for (int i=0;i<10; i++) {
        NNO_INTYPE in[INNODES];
        int testData = trainer.MakeData(in,true);
        NetworkTrainer::PARITY parity = trainer.Parity(in);
        double *answer = testNet->Recall(in);
        if (parity == NetworkTrainer::EVEN)
            cout <<  testData << "\t(parity=EVEN) is classified " << answer[0] << endl;
        else
            cout <<  testData << "\t(parity=ODD)  is classified " << answer[0] << endl;
    }
    
    return EXIT_SUCCESS;
}


NetworkTrainer::NetworkTrainer(int nEpoch) : fEpoch(nEpoch)
{
    // Set up training and test files for the parity problem
    
    cout << "Generation of training and test files..." << endl;
    int ntrn = 1000;
    int ntst = 100;
    cout << "Training samples=" << ntrn << " ; Test samples=" << ntst << endl;
    MakeFile("parity.trn",ntrn,false); // create training file
    MakeFile("parity.tst",ntst,false); // create testing file
}

void NetworkTrainer::Train()
{
    // Perform training cycles for all networks in the list
    // based on the parity file.
    // Show the results using the default plotter
    
    TIter iter(&fNets);
    VNeuralNet *net;
    while ((net = (VNeuralNet *) iter.Next()) != NULL) {
        if (PLOTS) net->SetupPlots();
        for (int i=0;i<fEpoch;i++) {
            net->TrainEpoch("parity.trn");
            net->TestEpoch("parity.tst");
        }
        gSystem->Sleep(1000);
    }
}


NetworkTrainer::PARITY NetworkTrainer::Parity(NNO_INTYPE* Data) 
{ 
    // computes parity: 1= odd, 0=even
    
    int I;
    int count=0;
    for(I=0; I<INNODES; ++I) if ( Data[I] > 0.5 ) ++count;
    return (PARITY) (count&1);   // last bit of 'Count' contains parity
}

int NetworkTrainer::MakeData(NNO_INTYPE* in, bool randomize) 
{
    // There are two ways of operation (according to the state of the 'randomize' flag)
    // 1) Create all pattern and set inputs correspondingly (Shuffle the data set)
    // 2) Create random numbers
    
    const int ntrn = 1<<INNODES; // Number of patterns
    static int samples[ntrn];
    static bool initialized = false;
    if (!initialized) {
        
        for (int i=0;i<ntrn;i++) { // Generate all patterns
            samples[i] = i;
        }
        
        for (int j=0;j<ntrn;j++) { // Shuffle them
            int index1 = rand()%ntrn;
            int index2 = rand()%ntrn;
            int tmp = samples[index1];
            samples[index1] = samples[index2];
            samples[index2] = tmp;
        }
        
        initialized = true;
    }
    
    // Choose an item and set up the input vector
    
    int value = 0;
    
    if (randomize) {
        value = rand()%ntrn;
    }
    else {
        static int n = 0;
        value = samples[n++];
        if (n>=ntrn) n = 0;
    }
    
    // Generate the input vector
    
    for (int I=0; I<INNODES; ++I) {
        int mask = 1<<I;
        in[I] = (NNO_INTYPE) ((value&mask)>0);
    }
    
    return value;
}

void NetworkTrainer::MakeFile(const char* filename,int Records, bool randomize) 
{
    // Write a file that contains pairs of input and output vectors
    // for training of the parity
    
    NNO_INTYPE Data[INNODES]; // input of network
    NNO_OUTTYPE Out;          // output to be trained (parity of input)
    
    // Out is 1 if parity is even, else -1
    
    int I;
    FILE* file = fopen(filename,"wb");
    for (I=0;I<Records;++I) {
        MakeData(Data, randomize);
        if (Parity(Data)==EVEN) Out = 1.0; else Out = -1.0;
        fwrite(Data,sizeof(NNO_INTYPE),INNODES,file);
        fwrite(&Out,sizeof(NNO_OUTTYPE),1,file);
    }
    fclose(file);
}
