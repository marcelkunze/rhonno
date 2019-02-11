// NetworkTrainer
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// M.Kunze, Bochum University, 1997-2001
///
// Usage: NetworkTrainer <file> <start> <end>
// file          = steering file
// start         = epoch to start with (Default 1)
// end           = epoch to end with (Default 200)

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"

#include "NetworkTrainer.h"
#include "TDataServe.h"
#include "TNNK.h"
#include "TFD.h"
#include "TMLP.h"
#include "TXMLP.h"
#include "TSGNG.h"
#include "TSGCS.h"
#include "TGNG.h"
#include "TGCS.h"
#include "TLVQ.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdlib>
using namespace std;

TROOT root("trainer","Network Trainer");
TApplication *theApp = new TApplication("Network Trainer",0,0);

int main(int argc,char* argv[])
{
    string steeringFile;
    if (argc>1) steeringFile = argv[1];
    
    int startEpoch =                1;       // First epoch and
    if (argc>2) startEpoch = atoi(argv[2]);
    
    int stopEpoch   =             100;       // Last epoch to be trained
    if (argc>3) stopEpoch = atoi(argv[3]);
    
    NetworkTrainer *trainer = new NetworkTrainer(steeringFile,startEpoch,stopEpoch);
    trainer->SetupDataServer();
    trainer->SetupNetworks();
    trainer->Train();
    trainer->PrintOn();
    
    string file("Recall");
    file =  file + trainer->GetModel() + ".cpp";
    trainer->WriteSourceCode(file);
    
    return EXIT_SUCCESS;
}

NetworkTrainer::NetworkTrainer(string file,int se,int ee)
: fPidDataServer(0), fTrainingServer(0), fVectorsEpoch(0), fNet(0), fMomentum(0.0), 
fTrnMax(0), fTstMax(1000), fInNodes(NNODIMENSION), fHid1Nodes(10), fHid2Nodes(1), fOutNodes(1),
fCells(1000), fBalance(false), fPlots(false), fScale(1.0), fAutoScale(false),
fTransfer(TNeuralNetParameters::TR_FERMI), fTarget(0)
{
    fStartEpoch = se;
    fStopEpoch = ee;
    
    cout << endl << "NNO NETWORK TRAINER"
    << endl << "-------------------" << endl;
    
    fModel = "TMLP";
    
    int i;
    for (i=0;i<NNODIMENSION;i++) fInScale[i] = 1.0;
    for (i=0;i<NNODIMENSION;i++) fOutScale[i] = 1.0;
    for (i=0;i<NNODIMENSION;i++) fInMean[i] = 0.0;
    for (i=0;i<NNODIMENSION;i++) fOutMean[i] = 0.0;
    
    string steeringFile = file;
    if (steeringFile!="") ReadSteeringFile(steeringFile);
    
    cout << "\n\nTraining " << fModel << " epoch: " << fStartEpoch << " - " << fStopEpoch << endl;
}

NetworkTrainer::~NetworkTrainer() 
{
    delete fPidDataServer;
    delete fTrainingServer;
    delete fNet;
    //fFile->Close(); delete fFile;
}

// Perform a training cycle

void NetworkTrainer::SetupDataServer(string file)
{
    string path(fDataPath);
    
    if (fTree!="") {
        
        cout << endl << "Reading data files: " << endl;
        
        if (fPidDataServer==0)
            fPidDataServer = new TDataServe("NetworkData","Network training data",fInNodes,1);
        
        int i;
        for (i=0;i<fAll.GetSize();i++) {
            string fileName = path + (((TObjString*)fAll.At(i))->GetString()).Data();
            cout << "+file " << fileName << endl;
            fPidDataServer->TTreeDataRead(fileName,fTree,fInput,fOutput,fCut);
        }
        for (i=0;i<fPro.GetSize();i++) {
            string fileName = path + (((TObjString*)fPro.At(i))->GetString()).Data();
            cout << "+sample " << fileName << endl;
            fPidDataServer->TTreeDataRead(fileName,fTree,fInput,"1",fCut);
        }
        for (i=0;i<fCon.GetSize();i++) {
            string fileName = path + (((TObjString*)fCon.At(i))->GetString()).Data();
            cout << "-sample " << fileName << endl;
            fPidDataServer->TTreeDataRead(fileName,fTree,fInput,"0",fCut);
        }
    }
    else {
        string fileName = path + file;
        fFile = new TFile(fileName.data(),"READ");
        fPidDataServer = (TDataServe*)fFile->Get("PidData");
    }
    
    fPidDataServer->Init(0);
    fPidDataServer->MixTrn();
    fTrnMax = fPidDataServer->GetNumTrnvecs();
    if (fVectorsEpoch<=0 || fTrnMax<fVectorsEpoch+fTstMax) fVectorsEpoch = fTrnMax-fTstMax;
    
    // Set up a server cache for training data
    
    int i;
    int j;
    
    if (fAutoScale) {
        
        cout << endl << "Scale input layer: " << endl;
        
        
        string input(fInput);
        replace( input.begin(), input.end(), ':', ' ');
        input += " ";
        istringstream inStream(input);
        
        float *inScale  = fPidDataServer->GetInputScale();
        for (j=0;j<fInNodes;j++) {
            fInScale[j] = inScale[j];
            string token;
            inStream >> token;
            cout << token << "\t*\t" << fInScale[j] << endl;
        }
        
        /*	cout << endl << "Scale factors for output:" << endl;
         float *outScale = fPidDataServer->GetOutputScale();
         for (j=0;j<fOutNodes;j++) {
         fOutScale[j] = outScale[j];
         cout << fOutScale[j] << "\t";
         }
         */
    }
    
    // Assemble a new training data set
    
    cout << endl << endl << "Setting up the normalised Server " << endl;
    
    fTrainingServer = new TDataServe("NetworkData","Network training data",fInNodes,fOutNodes);
    
    // Read the original vectors
    for (i=0;i<fVectorsEpoch+fTstMax;i++) {
        float *inv  = (float *)fPidDataServer->GetInvecTrn(i);
        float *outv = (float *)fPidDataServer->GetOutvecTrn(i);
        for (j=0;j<fInNodes;j++) fInVector[j] = (float) fInScale[j]*inv[j];
        for (j=0;j<fOutNodes;j++) fOutVector[j] = (float) fOutScale[j]*outv[j];
        fTrainingServer->Putvec(fInVector,fOutVector);
        if (i>0&&i%10000==0) cout << i << endl;
    }
    
    // Reserve test vecs
    fTrainingServer->Init(fTstMax);
    fTrnMax = fTrainingServer->GetNumTrnvecs();
    fTstMax = fTrainingServer->GetNumTstvecs();
    cout << fTrnMax << " Training samples; " << fTstMax << " Test samples" << endl;
    
}

void NetworkTrainer::SetupNetworks()
{
    // Initialise Network
    
    if (fStartEpoch==1) {
        
        cout << "Initializing " << fModel << " network" << endl;
        fNetworkFile = "NNO0001." + fModel;
        
        if (fModel == "TFD") {
            fNet = new TXMLP(1,fScale,fNetworkFile,fInNodes,fOutNodes,0.1,fTransfer);
        }
        else if (fModel == "TMLP") {
            fNet = new TMLP(0.1,0.01,fInNodes,fHid1Nodes,fOutNodes,fScale,fNetworkFile,fTransfer);
        }
        else if (fModel == "TXMLP") {
            fNet = new TXMLP(3,fScale,fNetworkFile,fInNodes,fHid1Nodes,fHid2Nodes,fOutNodes,0.1,0.02,0.01,TNeuralNetParameters::TR_FERMI,TNeuralNetParameters::TR_FERMI,fTransfer);
        }
        else if (fModel == "TNNK") {
            string hidden;
            hidden += fHid1Nodes;
            if (fHid2Nodes>1) {
                hidden += ":";
                hidden += fHid2Nodes;
            }
            Text_t *hid = (char *) hidden.data();
            fNet = new TNNK(0.2,0.0,fMomentum,fInNodes,hid,fOutNodes,fNetworkFile);
        }
        else if (fModel == "TSGNG") {
            fNet = new TSGNG(fInNodes,fOutNodes,200,0.1,0.02,0.1,0.01,0.01,0.01,0.01,10,5000,500, fNetworkFile);
        }
        else if (fModel == "TSGCS") {
            fNet = new TSGCS(fInNodes,   // number of inputnodes
                             4,                // starting number of cells (dimension of connection topologie is Cells-1)
                             fOutNodes,        // number of outputnodes
                             fCells,           // maximum number of cells
                             0.1,              // learningstep of winner cell
                             0.002,            // learningstep of neighbour
                             0.1,              // learningstep of neural weights
                             0.01,             // decrement of error_count: err_count*=(1.0 - AErrCount)
                             0.01,             // adjustment of squared deviation during learnig
                             10,               // maximum number of allowed connections for one cell
                             1000,             // cell insertion after 1000 Learningsteps
                             0,                // don't remove cells
                             fNetworkFile);    // network - filename (used by destructor)
        }
        else if (fModel == "TGNG") {
            fNet = new TGNG(	fInNodes,   // number of inputnodes
                            fCells,           // maximum number of cells
                            0.2,              // learningstep of winner cell
                            0.02,             // learningstep of neighbour
                            0.1,              // learningstep of neural weights
                            0.01,             // decrement of error_count: err_count*=(1.0 - AErrCount)
                            0.01,             // adjustment of squared deviation during learnig
                            20,               // maximum number of allowed connections for one cell
                            1000,             // cell insertion
                            0,                // cell removal
                            fNetworkFile);    // network - filename (used by destructor)
        }
        else if (fModel == "TGCS") {
            fNet = new TGCS( fInNodes,   // number of inputnodes
                            fInNodes+1,	    // number of initial cells (Dimension+1)
                            fCells,		    // stop growth
                            0.2,		    // Learnstep of Winner-Cell
                            0.02,		    // Learnstep of Neighbours
                            0.002,		    // a_win_count
                            20,		    // connectors
                            500,		    // insert_step
                            0,		    // delet_step
                            fNetworkFile);      // network - filename (used by destructor)
        }
        else if (fModel == "TLVQ") {
            fNet = new TLVQ( fInNodes,   // number of inputnodes
                            fCells,		    // cells
                            0.2,		    // Learnstep of Winner-Cell
                            fNetworkFile);      // network - filename (used by destructor)
        }
        else {
            cout << "Error initializing network. Unknown model " << fModel << endl;
            exit(0);
        }
    }
    
    else {
        cout << "Loading " << fModel << " network " << fNetworkFile << endl;
        if (fModel == "TFD"){
            fNet = new TFD(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TMLP"){
            fNet = new TMLP(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TXMLP"){
            fNet = new TXMLP(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TNNK"){
            fNet = new TNNK(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TSGNG"){
            fNet = new TSGNG(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TSGCS"){
            fNet = new TSGCS(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TGNG"){
            fNet = new TGNG(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TGCS"){
            fNet = new TGCS(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else if (fModel == "TLVQ"){
            fNet = new TLVQ(Makename(fStartEpoch , fNetworkPath, fNetworkFile));
        }
        else {
            cout<<"Error initializing network. Unknown Model."<<endl;
            exit(0);
        }
    }
    
    fNet->SetThreshold(0.5);
    fNet->SetMomentumTerm(fMomentum);
    fNet->BalanceSamples(fBalance);
}

double NetworkTrainer::Train()
{
    // Prepare drawing
    if (fPlots) fNet->SetupPlots();
    
    double error = 0.0;
    double errbest = FLT_MAX;
    int epobest = 0;
    
    for (int epo=fStartEpoch; epo<=fStopEpoch; epo++){
        
        cout << endl << "epoch: " << epo << " / best epoch: " << epobest << endl;
        
        // Perform training for an epoch
        error = fNet->TrainEpoch(fTrainingServer);
        
        // Save the networks after each epoch
        string network =  Makename(epo , fNetworkPath, fNetworkFile);
        cout << "Saving: " << network << endl;
        fNet->Save(network);
        
        // Test and save the best network
        error = Test();
        if (error < errbest) {
            epobest = epo;
            errbest = error;
            network =  Makename(0 , fNetworkPath, fNetworkFile);
            cout << "Saving as best network " << network << " with error " << error << "%" << endl;
            fNet->Save(network);
            if (errbest<fTarget) break; // limit reached
        }
        
        // Adapt the learning rate, freeze network upon convergence
        if (fModel == "xmlp") {
            if (error<adaptThreshold) {
                cout << '\t' << "Adapt rate: ";
                for (int i=0; i<3; i++) {
                    TNeuralNetParameters &parm = ((TXMLP*)fNet)->GetPerceptron(i)->GetParameters();
                    parm.fLearnStep *= adaptRate;	// Decrease leraning rate
                    cout << parm.fLearnStep << " ";
                }
            }
        }
        
    }
    
    cout << "Finished training with best epoch: " << epobest << endl;
    
    return error;
}

double NetworkTrainer::Test()
{
    double tst = 100. * fNet->TestEpoch(fTrainingServer) / fTstMax; // ok in percent
    return tst;
}

// Print a summary

void NetworkTrainer::PrintOn()
{
    cout << endl;
    cout << "Number of trainvecs:        " << fTrnMax << endl;
    cout << "Vectors per epoch:          " << fVectorsEpoch << endl;
    cout << "Last epoch:                 " << fStopEpoch << endl;
}


// Generate names with epoch number extension

string NetworkTrainer::Makename(int z, string networkPath, string name){
    
    string number = to_string(z);
    if (z<10)
        name = networkPath + "NNO000" + number + "." + fModel;
    else if (z<100)
        name = networkPath + "NNO00" + number + "." + fModel;
    else if (z<1000)
        name = networkPath + "NNO0" + number + "." + fModel;
    else
        name = networkPath + "NNO" + number + "." + fModel;
    
    return name;
}

bool NetworkTrainer::ReadSteeringFile(string filename)
{
    const int len(1024);
    Char_t skip[len];
    
    ifstream  s(filename, ios::in);
    
    if (!s) {
        cout << endl << "NetworkTrainer: Could not open " << filename << endl;
        return false;
    }
    else
        cout << endl << "NetworkTrainer: Reading parameters from " << filename << endl;
    
    while (!s.eof()) {
        string key,model;
        s >> key;  // Get key
        
        if (key == "") { s.getline(skip,len); } // empty line
        
        else if (key == "target") {
            s >> fTarget;
            cout << "+Training target " << fTarget << endl;
        }
        
        else if (key == "start") {
            s >> fStartEpoch;
            cout << "+Start epoch " << fStartEpoch << endl;
        }
        
        else if (key == "stop") {
            s >> fStopEpoch;
            cout << "+Stop epoch " << fStopEpoch << endl;
        }
        
        else if (key == "epoch") {
            s >> fVectorsEpoch;
            cout << "+Epoch size " << fVectorsEpoch << endl;
        }
        
        else if (key == "test") {
            s >> fTstMax;
            cout << "+Test samples " << fTstMax << endl;
        }
        
        else if (key == "scale") {
            s >> fScale;
            cout << "+Scale all inputs by " << fScale << endl;
        }
        
        else if (key == "autoscale") {
            string choice;
            s >> choice;
            if (choice=="true") {
                AutoScale();
                cout << "+Scale inputs by mean values " << endl;
            }
        }
        
        else if (key == "balance") {
            string choice;
            s >> choice;
            if (choice=="true") {
                BalanceSamples();
                cout << "+Balance the positive and negative samples " << endl;
            }
        }
        
        else if (key == "plots") {
            string choice;
            s >> choice;
            if (choice=="true") {
                ShowControlPlots();
                cout << "+Produce control plots " << endl;
            }
        }
        
        else if (key == "momentum") {
            s >> fMomentum;
            cout << "+Set momentum term " << fMomentum << endl;
        }
        
        else if (key == "transfer") {
            string trans;
            s >> trans;
            
            if (trans=="TR_FERMI")
                fTransfer = TNeuralNetParameters::TR_FERMI;
            else if (trans=="TR_SIGMOID")
                fTransfer = TNeuralNetParameters::TR_SIGMOID;
            else if (trans=="TR_LINEAR")
                fTransfer = TNeuralNetParameters::TR_LINEAR;
            else if (trans=="TR_LINEAR_BEND")
                fTransfer = TNeuralNetParameters::TR_LINEAR_BEND;
            else if (trans=="TR_USER")
                fTransfer = TNeuralNetParameters::TR_USER;
            else
                cerr << "Do not know how to set ";
            
            cout << "+Transfer functions " << trans << endl;
        }
        
        else if (key == "tree") {
            s >> fTree;
            cout << "+tree " << fTree << endl;
        }
        
        else if (key == "datapath") {
            s >> fDataPath;
            cout << "+datapath " << fDataPath << endl;
            fDataPath += "/";
        }
        
        else if (key == "networkpath") {
            s >> fNetworkPath;
            cout << "+networkpath " << fNetworkPath << endl;
            fNetworkPath += "/";
        }
        
        else if (key == "input") {
            s >> fInput;
            cout << "+input " << fInput << endl;
        }
        
        else if (key == "output") {
            s >> fOutput;
            cout << "+output " << fOutput << endl;
        }
        
        else if (key == "pro") {
            string pro;
            s >> pro;
            fPro.Add(new TObjString(pro.data()));
            cout << "+pro " << pro << endl;
        }
        
        else if (key == "con") {
            string con;
            s >> con;
            fCon.Add(new TObjString(con.data()));
            cout << "+con " << con << endl;
        }
        
        else if (key == "file") {
            string all;
            s >> all;
            fAll.Add(new TObjString(all.data()));
            cout << "+file " << all << endl;
        }
        
        else if (key == "inscale") {
            cout << "inscale: " << endl;
            for (int i=0;i<fInNodes;i++) {
                s >> fInScale[i];
                cout << "\t" <<  fInScale[i] << endl;
            }
            continue;
        }
        
        else if (key == "outscale") {
            cout << "outscale: " << endl;
            for (int i=0;i<fOutNodes;i++) {
                s >> fOutScale[i];
                cout << "\t" <<  fOutScale[i] << endl;
            }
            continue;
        }
        
        else if (key == "cut") {
            string cut;
            s >> cut;
            fCut = cut;
            cout << "+cut " << fCut << endl;
        }
        
        else if (key == "fisher") {
            fModel = "TFD";
            s >> fInNodes >> fOutNodes;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            if (fOutNodes>NNODIMENSION) { cerr << "Too much output nodes:" << fOutNodes << endl; return false; }
            cout << fModel << " " << fInNodes << "-" << fOutNodes << endl;
        }
        
        else if (key == "mlp") {
            fModel = "TMLP";
            s >> fInNodes >> fHid1Nodes >> fOutNodes;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            if (fOutNodes>NNODIMENSION) { cerr << "Too much output nodes:" << fOutNodes << endl; return false; }
            cout << fModel << " " << fInNodes << "-" << fHid1Nodes << "-" << fOutNodes << endl;
        }
        
        else if (key == "xmlp") {
            fModel = "TXMLP";
            s >> fInNodes >> fHid1Nodes >> fHid2Nodes >> fOutNodes;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            if (fOutNodes>NNODIMENSION) { cerr << "Too much output nodes:" << fOutNodes << endl; return false; }
            cout << fModel << " " << fInNodes << "-" << fHid1Nodes << "-" << fHid2Nodes << "-" << fOutNodes << endl;
        }
        
        else if (key == "tnnk") {
            fModel = "TNNK";
            s >> fInNodes >> fHid1Nodes >> fHid2Nodes >> fOutNodes;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            if (fOutNodes>NNODIMENSION) { cerr << "Too much output nodes:" << fOutNodes << endl; return false; }
            cout << fModel << " " << fInNodes << "-" << fHid1Nodes << "-" << fHid2Nodes << "-" << fOutNodes << endl;
        }
        
        else if (key == "sgng") {
            fModel = "TSGNG";
            s >> fInNodes >> fCells >> fOutNodes;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            if (fOutNodes>NNODIMENSION) { cerr << "Too much output nodes:" << fOutNodes << endl; return false; }
            cout << fModel << " " << fInNodes << "-" << fCells << "-" << fOutNodes << endl;
        }
        
        else if (key == "sgcs") {
            fModel = "TSGCS";
            s >> fInNodes >> fCells >> fOutNodes;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            if (fOutNodes>NNODIMENSION) { cerr << "Too much output nodes:" << fOutNodes << endl; return false; }
            cout << fModel << " " << fInNodes << "-" << fCells << "-" << fOutNodes << endl;
        }
        
        else if (key == "gng") {
            fModel = "TGNG";
            s >> fInNodes >> fCells;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            cout << fModel << " " << fInNodes << " " << fCells << endl;
        }
        
        else if (key == "gcs") {
            fModel = "TGCS";
            s >> fInNodes >> fCells;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            cout << fModel << " " << fInNodes << " " << fCells << endl;
        }
        
        else if (key == "lvq") {
            fModel = "TLVQ";
            s >> fInNodes >> fCells;
            if (fInNodes>NNODIMENSION) { cerr << "Too much input nodes:" << fInNodes << endl; return false; }
            cout << fModel << " " << fInNodes << " " << fCells << endl;
        }
        
        
        //else if (key.find("#")!=string::npos) { s.getline(skip,len);  }
        else if (key.substr(0,1)=="#") {
            s.getline(skip,len);
        }
    }
    
    s.close();
    return kTRUE;
}

void NetworkTrainer::WriteSourceCode(string filename)
{
    TDatime theTime;
    string input(fInput);
    replace( input.begin(), input.end(), ':', ' ');
    input += " ";
    istringstream inStream(input);
    
    std::ofstream f(filename);
    f << "// " << fModel << " network trained with NNO NetworkTrainer at " << theTime.AsString() << endl;
    f << "// Input parameters  " << fInput << endl;
    f << "// Output parameters " << fOutput << endl;
    int i;
    f << "// Training files:" << endl;
    for (i=0;i<fAll.GetSize();i++) {
        string fileName = fDataPath + (((TObjString*)fAll.At(i))->GetString()).Data();
        f << "//" << fileName << endl;
    }
    for (i=0;i<fPro.GetSize();i++) {
        string fileName = fDataPath + (((TObjString*)fPro.At(i))->GetString()).Data();
        f << "//" <<  fileName << endl;
    }
    for (i=0;i<fCon.GetSize();i++) {
        string fileName = fDataPath + (((TObjString*)fCon.At(i))->GetString()).Data();
        f << "//" <<  fileName << endl;
    }
    f << endl;
    f << "#include \"RhoNNO/" << fModel << ".h\"" << endl << endl;
    f << "double* Recall(double *invec)" << endl;
    f << "{" << endl;
    f << "\tstatic " << fModel << " net(\"" << fModel << ".net\");" << endl;
    f << "\tfloat x[" << fInNodes << "];" << endl;
    for (i=0;i<fInNodes;i++) {
        f << "\tx[" << i << "] \t= " << fInScale[i] << "\t*\tinvec[" << i << "];" ;
        string token;
        inStream >> token;
        f << "\t// " << token << endl;
    }
    f << "\treturn net.Recall(x);" << endl;
    f << "}" << endl;
}

