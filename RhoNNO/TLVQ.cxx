// TLVQ
//
// Implementation of the LEARNING-VECTOR-QUANTISATION (LVQ)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "RhoNNO/TLVQ.h"
#include "RhoNNO/VNeuralNetPlotter.h"

ClassImp(TLVQ)

#include <cfloat>
using namespace std;

TLVQ::TLVQ(int innodes,int cells,double winStep,string netFile)
: VUnsupervisedNet("LVQ",innodes,cells,netFile) 
{
    fXB.fCells   = cells;
    fXB.fWinStep = winStep;
    fU           = 0;
    AllocNet();
    InitNet();
}

// copy constructor
TLVQ::TLVQ(const TLVQ& lvq,string netFile)
: VUnsupervisedNet("LVQ",lvq.fParm.fInNodes,lvq.fParm.fOutNodes,netFile) 
{
    fXB = lvq.fXB;
    fU  = 0;
    AllocNet();
    InitNet();
    CopyData(lvq);
}

TLVQ::~TLVQ() 
{
    int I;
    if (fFilename!="") if (fShouldSave) Save();
    if (fU!=0) {
        TNeuralNetCell* up = fU;
        for (I=0;I<fParm.fOutNodes;++I) {
            delete[] up->fVector;
            delete[] up->fDiff;
            ++up;
        }
        delete[] fU;
    }
}


void TLVQ::ReadBinary(void) 
{
    fU = 0;
    fread(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    AllocNet();
    TNeuralNetCell* up;
    for(up=fU;up<fUbound;++up) {
        fread(up->fVector,sizeof(double),fParm.fInNodes,fFile);
        freadvar(up->fID);
    }
}


void TLVQ::ReadText(void) 
{
    fU = 0;
    fscanf(fFile,"\
           win_step     %le\n\
           cells        %i\n",
           &fXB.fWinStep,
           &fXB.fCells);
    AllocNet();
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        fscanf(fFile,"\n");
        fscanf(fFile,"TNeuralNetCell number      %i\n",&up->fID);
        fscanf(fFile,"class            %i\n",&up->fClass);
        fscanf(fFile,"vectors ");
        for (I=0;I<fParm.fInNodes;++I) fscanf(fFile,"%le ",&up->fVector[I]);
        fscanf(fFile,"\n");
    }
}

void TLVQ::WriteBinary(void) 
{
    TNeuralNetCell* up;
    fwrite(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    for(up=fU;up<fUbound;++up) {
        fwrite(up->fVector,sizeof(double),fParm.fInNodes,fFile);
        fwritevar(up->fID);
    }
}

void TLVQ::WriteText(void) 
{
    fprintf(fFile,"\
            win_step     %le\n\
            cells        %i\n",
            fXB.fWinStep,
            fXB.fCells);
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        fprintf(fFile,"\n");
        fprintf(fFile,"TNeuralNetCell number      %i\n",up->fID);
        fprintf(fFile,"class            %i\n",up->fClass);
        fprintf(fFile,"vectors ");
        for (I=0;I<fParm.fInNodes;++I) fprintf(fFile,"%le ",up->fVector[I]);
        fprintf(fFile,"\n");
    }
}

void TLVQ::AllocNet(void) 
{
    fU = new TNeuralNetCell[fParm.fOutNodes];  TestPointer(fU);
    fUbound = &fU[fXB.fCells];
    int I;
    TNeuralNetCell* up = fU;
    for (I=0;I<fParm.fOutNodes;++I) {
        up->fVector = new double[fParm.fInNodes];    TestPointer(up->fVector);
        up->fDiff = new double[fParm.fInNodes];    TestPointer(up->fDiff);
        up->fID=I;
        ++up;
    }
}

void TLVQ::InitNet(void) 
{
    TNeuralNetCell* up;
    int J;
    for(up=fU;up<fUbound;++up) {
        for (J=0;J<fParm.fInNodes;++J) up->fVector[J]=Random();
        up->fClass=0;
    }
}

void TLVQ::CopyData(const TLVQ& lvq) 
{
    TNeuralNetCell* thisup = fU;
    TNeuralNetCell* fromup = lvq.fU;
    int I;
    
    //check integrity
    if (fParm.fInNodes    !=lvq.fParm.fInNodes)  Errorf((char *)"cannot copy data; innodes not identical");
    if (fParm.fOutNodes   !=lvq.fParm.fOutNodes) Errorf((char *)"cannot copy data; outnodes not identical");
    
    fXB = lvq.fXB;
    fUbound = &fU[fXB.fCells];
    for (I=0;I<fXB.fCells;++I) {
        memcpy(thisup->fVector,fromup->fVector,sizeof(double)*fParm.fInNodes);
        thisup->fClass = fromup->fClass;
        ++thisup;
        ++fromup;
    }
}

int  TLVQ::GetWinnerCell(NNO_INTYPE* in) 
{
    int I,J;
    double s_dist;
    double min_s_dist = DBL_MAX;
    TNeuralNetCell* up;
    J=0;
    for(up=fU;up<fUbound;++up){
        double* v = up->fVector;
        double* d = up->fDiff;
        NNO_INTYPE* i = in;
        s_dist = 0.0;
        for (I=0;I<fParm.fInNodes;++I) { *d =*i++ - *v++; s_dist+=*d * *d; ++d;}
        fOut[J++] = s_dist;
        if (s_dist<min_s_dist) {
            min_s_dist=s_dist;
            fUwin=up;
        }
    }
    
    if (fPlotter) fPlotter->AddTestSample(min_s_dist);
    
    return fUwin->fID;
}

double  TLVQ::Train(NNO_INTYPE* in,NNO_OUTTYPE*) 
{
    int J;
    Recall(in);  //make output of all cells and neurons find the winner
    double* vwin = fUwin->fVector;
    double* dwin = fUwin->fDiff;
    for (J=0;J<fParm.fInNodes;++J) *vwin++ += *dwin++ * fXB.fWinStep;
    fShouldSave = kTRUE;
    return fUwin->fID;
}

void TLVQ::Implant(int nr,int c,NNO_INTYPE* in) 
{
    int I;
    fU[nr].fClass = c;
    double* v=fU[nr].fVector;
    for (I=0;I<fParm.fInNodes;++I) *v++ = *in++;
    fShouldSave = kTRUE;
}
