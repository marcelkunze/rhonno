//////////////////////////////////////////////////////////////////////////
//									//
// TPerceptron								//
//									//
// Implementation of the perceptron (Supervised Learning)		//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//	 
//////////////////////////////////////////////////////////////////////////	 

#include "TMath.h"
#include "RhoNNO/TPerceptron.h"

ClassImp(TPerceptron)

// Transferfunctions
void TransferFermi(Double_t in,Double_t* out,Double_t* deriv) 
{
    if (in < -10.0) in = -10.0;
    if (in >  10.0) in =  10.0;
    Double_t O = 1. / (1. + exp(-in));
    *out   = O;
    *deriv = O * (1. - O);
}

void TransferSigmoid(Double_t in,Double_t* out,Double_t* deriv) 
{
    if (in < -10.0) in = -10.0;
    if (in >  10.0) in =  10.0;
    Double_t O = 1. - 2. / (1. + exp(-in));
    *out   = O;
    *deriv = 2. * O * (1. - O);
}

void TransferLinear(Double_t in,Double_t* out,Double_t* deriv) 
{
    *out = in;
    *deriv = 1.0;
}

void TransferLinearBend(Double_t in,Double_t* out,Double_t* deriv) 
{
    if (in < -1.0) {
	*out = -0.9 + in * 0.1;
	*deriv = 0.1;
    } else
	if (in > 1.0) {
	    *out = 0.9 + in * 0.1;
	    *deriv = 0.1;
	} else {
	    *out = in;
	    *deriv = 1.0;
	}
}

TPerceptron::TPerceptron(Int_t inNodes,
			 Int_t outNodes,
			 Double_t learnStep,
			 TNeuralNetParameters::TRANSFER transferId,
			 Int_t perceptronId) 
{
    
    fParm.fInNodes      = inNodes;
    fParm.fOutNodes     = outNodes;
    fParm.fLearnStep    = learnStep;
    fParm.fTransferId   = transferId;
    fParm.fPerceptronId = perceptronId;
    fPrev     = 0;
    AllocNet();
    InitNet();
}

TPerceptron::TPerceptron(TPerceptron* prev,
			 Int_t outNodes,
			 Double_t learnStep,
			 TNeuralNetParameters::TRANSFER transferId,
			 Int_t perceptronId) 
{
    
    fParm.fInNodes      = prev->fParm.fOutNodes;
    fParm.fOutNodes     = outNodes;
    fParm.fLearnStep    = learnStep;
    fParm.fTransferId   = transferId;
    fParm.fPerceptronId = perceptronId;
    fPrev             = prev;
    AllocNet();
    InitNet();
}

TPerceptron::TPerceptron(void) 
{
    fPrev = 0;
}

TPerceptron::TPerceptron(TPerceptron* prev) 
{
    fPrev = prev;
}

void TPerceptron::AllocNet(void) 
{
    int I;
    fU       = new PerceptronUnit[fParm.fOutNodes]; TestPointer(fU);
    fOut     = new Double_t[fParm.fOutNodes]; TestPointer(fOut);
    fDiffSrc = new Double_t[fParm.fOutNodes]; TestPointer(fDiffSrc);
    if (fPrev == 0) {
	fIn      = new Double_t[fParm.fInNodes]; TestPointer(fIn);
	fDiffDst = 0;
    } else {
	fIn      = fPrev->fOut;
	fDiffDst = fPrev->fDiffSrc;
    }
    fUbound = &fU[fParm.fOutNodes];
    PerceptronUnit* up = fU;
    for (I=0;I<fParm.fOutNodes;++I) {
	up->fVector = new Double_t[fParm.fInNodes]; 
	TestPointer(up->fVector);
	up->fDelta = new Double_t[fParm.fInNodes]; 
	TestPointer(up->fDelta);
	up->fThreshold = 0;
	up->fID = I;
	++up;
    }
    switch (fParm.fTransferId) {
    case TNeuralNetParameters::TR_FERMI  :      Transfer=TransferFermi;      break;
    case TNeuralNetParameters::TR_SIGMOID:      Transfer=TransferSigmoid;      break;
    case TNeuralNetParameters::TR_LINEAR :      Transfer=TransferLinear;     break;
    case TNeuralNetParameters::TR_LINEAR_BEND:  Transfer=TransferLinearBend; break;
    default:             Transfer=0;
    }
}

void TPerceptron::InitNet(void) 
{
    PerceptronUnit* up;
    int J;
    for(up=fU;up<fUbound;++up) {
	for (J=0;J<fParm.fInNodes;++J) {
	    up->fVector[J] = Random();
	    up->fDelta[J]  = 0.0;
	}
	up->fThreshold = Random();
    }
}


TPerceptron ::~TPerceptron() 
{
    PerceptronUnit* up = fU;
    if (fU!=0)  {
	for(up=fU;up<fUbound;++up) {
	    delete[] up->fVector;
	    delete[] up->fDelta;
	}
	delete[] fU;
    }
    delete[] fOut;
    delete[] fDiffSrc;
    if (fPrev==0) delete[] fIn;
}


void TPerceptron::WriteBinary() 
{
    PerceptronUnit* up;
    fwrite(&fParm,sizeof(PerceptronBase),1,fFile);
    for(up=fU;up<fUbound;++up) {
	fwrite(up->fVector,sizeof(Double_t),fParm.fInNodes,fFile);
	fwritevar(up->fThreshold);
	fwritevar(up->fID);
    }
}

void TPerceptron::ReadBinary() 
{
    PerceptronUnit* up;
    fread(&fParm,sizeof(PerceptronBase),1,fFile);
    AllocNet();
    for(up=fU;up<fUbound;++up) {
	fread(up->fVector,sizeof(Double_t),fParm.fInNodes,fFile);
	freadvar(up->fThreshold);
	freadvar(up->fID);
    }
}

void TPerceptron::WriteText() 
{
    fprintf(fFile,"\nPerceptron ID %i\n",fParm.fPerceptronId);
    fprintf(fFile,"innodes     %i\n",fParm.fInNodes);
    fprintf(fFile,"outnodes    %i\n",fParm.fOutNodes);
    fprintf(fFile,"learn_step  %le\n",fParm.fLearnStep);
    fprintf(fFile,"transfer_id %i\n",fParm.fTransferId);
    PerceptronUnit* up;
    int I;
    for(up=fU;up<fUbound;++up) {
	fprintf(fFile,"\n");
	fprintf(fFile,"unit number      %i\n",up->fID);
	fprintf(fFile,"threshold        %le\n",up->fThreshold);
	fprintf(fFile,"weights\n");
	for (I=0;I<fParm.fInNodes;++I) fprintf(fFile,"%le\n",up->fVector[I]);
	fprintf(fFile,"\n");
    }
}

void TPerceptron::ReadText() 
{
    fscanf(fFile,"\nPerceptron ID %i\n",&fParm.fPerceptronId);
    fscanf(fFile,"innodes     %i\n",&fParm.fInNodes);
    fscanf(fFile,"outnodes    %i\n",&fParm.fOutNodes);
    fscanf(fFile,"learn_step  %le\n",&fParm.fLearnStep);
    fscanf(fFile,"transfer_id %i\n",(int *)&fParm.fTransferId);
    AllocNet();
    PerceptronUnit* up;
    int I;
    for(up=fU;up<fUbound;++up) {
	fscanf(fFile,"\n");
	fscanf(fFile,"unit number      %i\n",&up->fID);
	fscanf(fFile,"threshold        %le\n",&up->fThreshold);
	fscanf(fFile,"weights\n");
	for (I=0;I<fParm.fInNodes;++I) fscanf(fFile,"%le\n",&up->fVector[I]);
	fscanf(fFile,"\n");
    }
}

Double_t* TPerceptron::Recall(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    int I;
    PerceptronUnit* up;
    Double_t* o = fOut;
    Double_t* ds = fDiffSrc;
    if (Transfer==0) Errorf((char *)"(TPerceptron) undefined transferfunction");
    for(up=fU;up<fUbound;++up) {
	Double_t* v = up->fVector;
	Double_t* i = fIn;
	Double_t  sum = 0.0;
	for (I=0;I<fParm.fInNodes;++I) sum += *i++ * *v++;
	sum -= up->fThreshold;
	Transfer(sum,o++,ds++);
    }

    return o;
}

Double_t TPerceptron::Train(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    int I;
    PerceptronUnit* up;
    Double_t* ds;
    
    // modify weights
    ds = fDiffSrc;
    for(up=fU;up<fUbound;++up) {
	Double_t* i = fIn;
	Double_t* v = up->fVector;
	Double_t* m = up->fDelta;
	for (I=0;I<fParm.fInNodes;++I) {
	    Double_t delta = *i++ * *ds;
	    *v++ += (delta + (*m * fParm.fMu)) * fParm.fLearnStep;
	    *m++ = delta;
	}
	up->fThreshold -= *ds * fParm.fLearnStep;
	++ds;
    }
    
    // propagate derivation backward if previous perceptron exists
    if (fDiffDst!=0) {
	Double_t diff;
	Double_t* dd = fDiffDst;
	for (I=0;I<fParm.fInNodes;++I) {
	    diff = 0.0;
	    ds = fDiffSrc;
	    for(up=fU;up<fUbound;++up) diff += up->fVector[I] * *ds++;
	    *dd++ *= diff;
	}
    }
    
    return 0.0;
}


