// TXMLP
//
// Implementation of the Multi-Layer-Perceptron
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TXMLP.h"
#include "VNeuralNetPlotter.h"

ClassImp(TXMLP)

#include <iostream>
using namespace std;

TXMLP:: TXMLP(int layers,double inputRange,string netFile,int innodes,...)
: VSupervisedNet("XMLP",innodes,0,netFile) 
{    
    int I;
    if (layers<1) Errorf((char *)"(TXMLP) at least one layer neccessary");
    fParm.fLayers  = layers;
    fParm.fInScale = 1.0/inputRange;
    
    int*    nodes = new int   [fParm.fLayers];      TestPointer(nodes);
    double* step  = new double[fParm.fLayers];      TestPointer(step);
    TNeuralNetParameters::TRANSFER* func  = new TNeuralNetParameters::TRANSFER[fParm.fLayers];      TestPointer(func);
    fPerc           = new TPerceptron*[fParm.fLayers]; TestPointer(fPerc);
    
    va_list ap;
    va_start(ap,innodes);
    
    for (I=0;I<fParm.fLayers;++I) nodes[I] = va_arg(ap,int);
    for (I=0;I<fParm.fLayers;++I) step[I]  = va_arg(ap,double);
    for (I=0;I<fParm.fLayers;++I) func[I]  = (TNeuralNetParameters::TRANSFER) va_arg(ap,int);
    
    va_end(ap);
    
    fParm.fOutNodes = nodes[layers-1];
    
    for (I=0;I<fParm.fLayers;++I) {
        if (I==0)
            fPerc[I] = new TPerceptron(fParm.fInNodes,nodes[I],step[I],func[I],I);
        else
            fPerc[I] = new TPerceptron(fPerc[I-1],nodes[I],step[I],func[I],I);
        
        TestPointer(fPerc[I]);
    }
    fOut = new double[fParm.fOutNodes]; // As we did not know before, allocate here...
    TestPointer(fOut);
    
    delete[] nodes;
    delete[] step;
    delete[] func;
}

TXMLP:: TXMLP(int layers,double inputRange,string netFile,int innodes,int n0,int n1,int n2,double s0,double s1,double s2,
              TNeuralNetParameters::TRANSFER f0,TNeuralNetParameters::TRANSFER f1,TNeuralNetParameters::TRANSFER f2)
: VSupervisedNet("XMLP",innodes,0,netFile) 
{    
    int I;
    if (layers!=3) Errorf((char *)"(TXMLP) Constructor needs 3 Layers");
    fParm.fLayers  = layers;
    fParm.fInScale = 1.0/inputRange;
    
    int*    nodes = new int   [fParm.fLayers];      TestPointer(nodes);
    double* step  = new double[fParm.fLayers];      TestPointer(step);
    TNeuralNetParameters::TRANSFER* func  = new TNeuralNetParameters::TRANSFER[fParm.fLayers];      TestPointer(func);
    fPerc           = new TPerceptron*[fParm.fLayers];  TestPointer(fPerc);
    
    nodes[0] = n0;nodes[1] = n1;nodes[2] = n2;
    step[0]  = s0;step[1]  = s1;step[2]  = s2;
    func[0]  = f0;func[1]  = f1;func[2]  = f2;
    
    fParm.fOutNodes = nodes[layers-1];
    
    for (I=0;I<fParm.fLayers;++I) {
        if (I==0)
            fPerc[I] = new TPerceptron(fParm.fInNodes,nodes[I],step[I],func[I],I);
        else
            fPerc[I] = new TPerceptron(fPerc[I-1],nodes[I],step[I],func[I],I);
        
        TestPointer(fPerc[I]);
    }
    fOut = new double[fParm.fOutNodes]; // As we did not know before, allocate here...
    TestPointer(fOut);
    
    delete[] nodes;
    delete[] step;
    delete[] func;
}

void TXMLP:: AllocNet(void) 
{
    fPerc = new TPerceptron*[fParm.fLayers]; TestPointer(fPerc); // MK: Allocation
    int I;
    for (I=0;I<fParm.fLayers;++I) {
        if (I==0)
            fPerc[I] = new TPerceptron();
        else
            fPerc[I] = new TPerceptron(fPerc[I-1]);
        
        TestPointer(fPerc[I]);
    }
}

TXMLP::~TXMLP() 
{
    if (fFilename!="") if (fShouldSave) Save();
    int I;
    for (I=0;I<fParm.fLayers;++I) delete fPerc[I];
    delete[] fPerc;
}

void TXMLP::ReadBinary(void) 
{
    int I;
    fread(&fParm,sizeof(TNeuralNetParameters),1,fFile);
    AllocNet();
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->ReadBinary();
    }
}

void  TXMLP::ReadText(void) 
{
    int I;
    int layers;
    double scale;
    fscanf(fFile,"layers    %i\n",&layers);
    fParm.fLayers = layers;
    fscanf(fFile,"in_scale  %le\n",&scale);
    fParm.fInScale = scale;
    AllocNet();
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->ReadText();
    }
}

void TXMLP::WriteBinary(void) 
{
    int I;
    fwrite(&fParm,sizeof(TNeuralNetParameters),1,fFile);
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->WriteBinary();
    }
}

void  TXMLP::WriteText(void) 
{
    int I;
    fprintf(fFile,"layers    %i\n",fParm.fLayers);
    fprintf(fFile,"in_scale  %le\n",fParm.fInScale);
    for (I=0;I<fParm.fLayers;++I) {
        fPerc[I]->SetFile(fFile);
        fPerc[I]->WriteText();
    }
}


double* TXMLP::Recall(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    int I;
    
    // convert input - optimized version
    NNO_INTYPE* i = in;
    double* pi = fPerc[0]->fIn;
    const double scale = fParm.fInScale;
    
    // Unroll input scaling loop for better performance
    for (I = 0; I < fParm.fInNodes - 3; I += 4) {
        pi[I]   = i[I]   * scale;
        pi[I+1] = i[I+1] * scale;
        pi[I+2] = i[I+2] * scale;
        pi[I+3] = i[I+3] * scale;
    }
    // Handle remaining elements
    for(; I < fParm.fInNodes; ++I) {
        pi[I] = i[I] * scale;
    }
    
    // recallstep of each perceptron
    for (I=0;I<fParm.fLayers;++I) fPerc[I]->Recall();
    for (I=0;I<fParm.fOutNodes;++I) fOut[I] = fPerc[fParm.fLayers-1]->fOut[I];
    
    if (fPlotter) {
        bool good = true;
        if (out!=0) good = out[0]>fParm.fThreshold;
        fPlotter->AddTestSample(fOut[0],good);
    }
    
    return fOut;
}

double TXMLP::Train(NNO_INTYPE* in,NNO_OUTTYPE* trout) 
{
    int I,J;
    fShouldSave = true;
    
    // convert input - optimized version
    NNO_INTYPE* i = in;
    double* pi = fPerc[0]->fIn;
    const double scale = fParm.fInScale;
    
    // Unroll input scaling loop for better performance
    for (I = 0; I < fParm.fInNodes - 3; I += 4) {
        pi[I]   = i[I]   * scale;
        pi[I+1] = i[I+1] * scale;
        pi[I+2] = i[I+2] * scale;
        pi[I+3] = i[I+3] * scale;
    }
    // Handle remaining elements
    for(; I < fParm.fInNodes; ++I) {
        pi[I] = i[I] * scale;
    }
    
    // recallstep of each perceptron
    for (I=0;I<fParm.fLayers;++I) fPerc[I]->Recall();
    for (I=0;I<fParm.fOutNodes;++I) fOut[I] = fPerc[fParm.fLayers-1]->fOut[I];
    
    // Calculate squared error - optimized
    double S_Err = 0.0;
    double* d = fPerc[fParm.fLayers-1]->fDiffSrc;
    double* out = fPerc[fParm.fLayers-1]->fOut;
    NNO_OUTTYPE* tr_out = trout;
    
    for (J = 0; J < fParm.fOutNodes; ++J) {
        double diff = *tr_out++ - *out++;
        *d++ = diff;
        S_Err += diff * diff;
    }
    
    for (I=fParm.fLayers-1;I>=0;--I) fPerc[I]->Train();
    
    if (fPlotter) fPlotter->AddTrainSample(trout[0],trout[0]>fParm.fThreshold);
    
    return S_Err;
}

void TXMLP::SetMomentumTerm(double f)
{
    for (int I=0;I<fParm.fLayers;++I) {
        fParm.fMu = f;
        TNeuralNetParameters &parm = fPerc[I]->GetParameters();
        parm.fMu = f;
    }
}
