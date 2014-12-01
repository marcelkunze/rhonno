//////////////////////////////////////////////////////////////////////////
//									//
// TXMLP								//
//									//
// Implementation of the Multi-Layer-Perceptron				//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/TXMLP.h"
#include "RhoNNO/VNeuralNetPlotter.h"

ClassImp(TXMLP)

#include <iostream>
using namespace std;

TXMLP:: TXMLP(Int_t layers,Double_t inputRange,const char* netFile,Int_t innodes,...)
: VSupervisedNet("XMLP",innodes,0,netFile) 
{    
    Int_t I;
    if (layers<1) Errorf((char *)"(TXMLP) at least one layer neccessary");
    fParm.fLayers  = layers;
    fParm.fInScale = 1.0/inputRange;
    
    Int_t*    nodes = new Int_t   [fParm.fLayers];      TestPointer(nodes);
    Double_t* step  = new Double_t[fParm.fLayers];      TestPointer(step);
    TNeuralNetParameters::TRANSFER* func  = new TNeuralNetParameters::TRANSFER[fParm.fLayers];      TestPointer(func);
    fPerc           = new TPerceptron*[fParm.fLayers]; TestPointer(fPerc);
    
    va_list ap;
    va_start(ap,innodes);
    
    for (I=0;I<fParm.fLayers;++I) nodes[I] = va_arg(ap,Int_t);
    for (I=0;I<fParm.fLayers;++I) step[I]  = va_arg(ap,Double_t);
    for (I=0;I<fParm.fLayers;++I) func[I]  = (TNeuralNetParameters::TRANSFER) va_arg(ap,Int_t);
    
    va_end(ap);
    
    fParm.fOutNodes = nodes[layers-1];
    
    for (I=0;I<fParm.fLayers;++I) {
	if (I==0)
	    fPerc[I] = new TPerceptron(fParm.fInNodes,nodes[I],step[I],func[I],I);
	else
	    fPerc[I] = new TPerceptron(fPerc[I-1],nodes[I],step[I],func[I],I);

	TestPointer(fPerc[I]);
    }
    fOut = new Double_t[fParm.fOutNodes]; // As we did not know before, allocate here...
    TestPointer(fOut);
    
    delete[] nodes;
    delete[] step;
    delete[] func;
}

TXMLP:: TXMLP(Int_t layers,Double_t inputRange,const char* netFile,Int_t innodes,Int_t n0,Int_t n1,Int_t n2,Double_t s0,Double_t s1,Double_t s2,
	      TNeuralNetParameters::TRANSFER f0,TNeuralNetParameters::TRANSFER f1,TNeuralNetParameters::TRANSFER f2)
: VSupervisedNet("XMLP",innodes,0,netFile) 
{    
    Int_t I;
    if (layers!=3) Errorf((char *)"(TXMLP) Constructor needs 3 Layers");
    fParm.fLayers  = layers;
    fParm.fInScale = 1.0/inputRange;
    
    Int_t*    nodes = new Int_t   [fParm.fLayers];      TestPointer(nodes);
    Double_t* step  = new Double_t[fParm.fLayers];      TestPointer(step);
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
    fOut = new Double_t[fParm.fOutNodes]; // As we did not know before, allocate here...
    TestPointer(fOut);
    
    delete[] nodes;
    delete[] step;
    delete[] func;
}

void TXMLP:: AllocNet(void) 
{
    fPerc = new TPerceptron*[fParm.fLayers]; TestPointer(fPerc); // MK: Allocation
    Int_t I;
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
    Int_t I;
    for (I=0;I<fParm.fLayers;++I) delete fPerc[I];
    delete[] fPerc;
}

void TXMLP::ReadBinary(void) 
{
    Int_t I;
    fread(&fParm,sizeof(TNeuralNetParameters),1,fFile);
    AllocNet();
    for (I=0;I<fParm.fLayers;++I) {
	fPerc[I]->SetFile(fFile);
	fPerc[I]->ReadBinary();
    }
}

void  TXMLP::ReadText(void) 
{
    Int_t I;
    Int_t layers;
    Double_t scale;
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
    Int_t I;
    fwrite(&fParm,sizeof(TNeuralNetParameters),1,fFile);
    for (I=0;I<fParm.fLayers;++I) {
	fPerc[I]->SetFile(fFile);
	fPerc[I]->WriteBinary();
    }
}

void  TXMLP::WriteText(void) 
{
    Int_t I;
    fprintf(fFile,"layers    %i\n",fParm.fLayers);
    fprintf(fFile,"in_scale  %le\n",fParm.fInScale);
    for (I=0;I<fParm.fLayers;++I) {
	fPerc[I]->SetFile(fFile);
	fPerc[I]->WriteText();
    }
}


Double_t* TXMLP::Recall(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    Int_t I;
    
    // convert input
    NNO_INTYPE* i=in;
    Double_t* pi = fPerc[0]->fIn;
    for (I=0;I<fParm.fInNodes;++I) *pi++ = *i++ * fParm.fInScale;
    
    // recallstep of each perceptron
    for (I=0;I<fParm.fLayers;++I) fPerc[I]->Recall();
    for (I=0;I<fParm.fOutNodes;++I) fOut[I] = fPerc[fParm.fLayers-1]->fOut[I];

    if (fPlotter) {
	Bool_t good = kTRUE;
	if (out!=0) good = out[0]>fParm.fThreshold;
	fPlotter->AddTestSample(fOut[0],good);
    }

    return fOut;
}

Double_t TXMLP::Train(NNO_INTYPE* in,NNO_OUTTYPE* trout) 
{
    Int_t I,J;
    fShouldSave = kTRUE;

    // convert input
    NNO_INTYPE* i=in;
    Double_t* pi = fPerc[0]->fIn;
    for (I=0;I<fParm.fInNodes;++I) *pi++ = *i++ * fParm.fInScale;
    
    // recallstep of each perceptron
    for (I=0;I<fParm.fLayers;++I) fPerc[I]->Recall();
    for (I=0;I<fParm.fOutNodes;++I) fOut[I] = fPerc[fParm.fLayers-1]->fOut[I];
    
    Double_t S_Err = 0;
    Double_t*   d = fPerc[fParm.fLayers-1]->fDiffSrc;
    Double_t* out = fPerc[fParm.fLayers-1]->fOut;
    NNO_OUTTYPE* tr_out=trout;
    for (J=0;J<fParm.fOutNodes;++J) { 
	*d = *tr_out++ - *out++; 
	S_Err += *d * *d; 
	d++; 
    }

    for (I=fParm.fLayers-1;I>=0;--I) fPerc[I]->Train();

    if (fPlotter) fPlotter->AddTrainSample(trout[0],trout[0]>fParm.fThreshold);

    return S_Err;
}

void TXMLP::SetMomentumTerm(Double_t f)
{
    for (int I=0;I<fParm.fLayers;++I) {
	fParm.fMu = f;
	TNeuralNetParameters &parm = fPerc[I]->GetParameters();
	parm.fMu = f;
    }
}
