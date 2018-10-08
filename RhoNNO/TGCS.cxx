// TGCS
//
// Implementation of the GROWING-CELL-STRUCTURE (GCS)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "TGCS.h"
#include "VNeuralNetPlotter.h"

ClassImp(TGCS)

#include <cfloat>
using namespace std;

//cell states
#define REMOVE 1

TGCS::TGCS (int innodes,int cells,int maxCells,double winStep,double neiStep,
            double aWinCount,int connectors,
            Long_t insertStep,Long_t deleteStep,string netFile)
: VUnsupervisedNet("GCS",innodes,maxCells,netFile)
{
    fXB.fCells       = cells;
    fXB.fMinCells    = cells;
    fXB.fWinStep     = winStep;
    fXB.fNeiStep     = neiStep;
    fXB.fWinCount    = aWinCount;
    fXB.fConnectors  = connectors;
    fXB.fInsertStep  = insertStep;
    fXB.fDeleteStep  = deleteStep;
    fXB.fInsertCount = 0;
    fXB.fDeleteCount = 0;
    fU               = 0;
    AllocNet();
    InitNet();
}

// copy constructor
TGCS::TGCS(const TGCS& gcs,string netFile)
: VUnsupervisedNet("GCS",gcs.fParm.fInNodes,gcs.fParm.fOutNodes,netFile) 
{
    fXB = gcs.fXB;
    fU  = 0;
    AllocNet();
    InitNet();
    CopyData(gcs);
}

TGCS::~TGCS() 
{
    Deviation();
    if (fFilename!="") if (fShouldSave) Save();
    if (fU != 0) {
        TNeuralNetCell* up;
        for(up=fU;up<fUbound;++up) {
            delete[] up->fVector;
            delete[] up->fDiff;
            delete[] up->fC;
        }
        delete[] fU;
    }
}

void TGCS::ReadBinary(void) 
{
    fU = 0;
    fread(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    AllocNet();
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::ReadUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
        freadvar(up->fClass);
        for (I=0;I<up->fNc;++I) up->fC[I].fPtr=&fU[up->fC[I].fID];
    }
}

void TGCS::ReadText(void) 
{
    fU = 0;
    fscanf(fFile,"win_step    %le\n",&fXB.fWinStep);
    fscanf(fFile,"new_step    %le\n",&fXB.fNeiStep);
    fscanf(fFile,"a_win_count %le\n",&fXB.fWinCount);
    fscanf(fFile,"cells       %i\n",&fXB.fCells);
    fscanf(fFile,"min_cells   %i\n",&fXB.fMinCells);
    fscanf(fFile,"connectors  %i\n",&fXB.fConnectors);
    fscanf(fFile,"insert_step %li\n",&fXB.fInsertStep);
    fscanf(fFile,"delete_step %li\n",&fXB.fDeleteStep);
    fscanf(fFile,"ins_count   %li\n",&fXB.fInsertCount);
    fscanf(fFile,"del_count   %li\n",&fXB.fDeleteCount);
    
    AllocNet();
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::ReadUnitText(fFile,(TNeuralNetCell*)up,&fParm);
        fscanf(fFile,"\nclass            %i\n",&up->fClass);
        for (I=0;I<up->fNc;++I)
            up->fC[I].fPtr = &fU[up->fC[I].fID];
    }
}

void TGCS::WriteBinary(void) 
{
    TNeuralNetCell* up;
    fwrite(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::WriteUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
        fwritevar(up->fClass);
    }
}

void TGCS::WriteText(void) 
{
    fprintf(fFile,"win_step    %le\n",fXB.fWinStep);
    fprintf(fFile,"nei_step    %le\n",fXB.fNeiStep);
    fprintf(fFile,"a_win_count %le\n",fXB.fWinCount);
    fprintf(fFile,"cells       %i\n",fXB.fCells);
    fprintf(fFile,"min_cells   %i\n",fXB.fMinCells);
    fprintf(fFile,"connectors  %i\n",fXB.fConnectors);
    fprintf(fFile,"insert_step %li\n",fXB.fInsertStep);
    fprintf(fFile,"delete_step %li\n",fXB.fDeleteStep);
    fprintf(fFile,"ins_count   %li\n",fXB.fInsertCount);
    fprintf(fFile,"del_count   %li\n",fXB.fDeleteCount);
    
    TNeuralNetCell* up;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::WriteUnitText(fFile,(TNeuralNetCell*)up,&fParm);
        fprintf(fFile,"\nclass            %i\n",up->fClass);
    }
}

void TGCS::AllocNet(void) 
{
    fU = new TNeuralNetCell[fParm.fOutNodes]; TestPointer(fU);
    fUbound = &fU[fXB.fCells];
    int I;
    TNeuralNetCell* up = fU;
    for (I=0;I<fParm.fOutNodes;++I) {
        up->fVector = new double[fParm.fInNodes]; TestPointer(up->fVector);
        up->fDiff   = new double[fParm.fInNodes]; TestPointer(up->fDiff);
        up->fC      = new connector[fXB.fConnectors]; TestPointer(up->fC);
        up->fNc     = 0;
        up->fCount  = 0;
        up->fID     = I;
        up->fState  = 0;
        ++up;
    }
}

void TGCS::InitNet(void) 
{
    TNeuralNetCell* up;
    int I,J;
    for(up=fU;up<fUbound;++up) {
        for (J=0;J<fParm.fInNodes;++J) up->fVector[J] = Random();
        up->fNc = 0;
        for (J=0;J<fXB.fCells;++J)
            if (&fU[J]!=up) up->fC[up->fNc++].fPtr = &fU[J];
    }
    for (I=0;I<fParm.fOutNodes;++I) fU[I].fClass = 0;
}

void TGCS::CopyData(const TGCS& gcs) 
{
    TNeuralNetCell* thisup = fU;
    TNeuralNetCell* fromup = gcs.fU;
    int I,J;
    
    //check integrity
    if (fParm.fInNodes    !=gcs.fParm.fInNodes)     Errorf((char *)"cannot copy data: fInNodes not identical");
    if (fParm.fOutNodes   !=gcs.fParm.fOutNodes)    Errorf((char *)"cannot copy data: fOutNodes not identical");
    if (fXB.fConnectors<gcs.fXB.fConnectors)  Errorf((char *)"cannot copy data: max fConnectors too low");
    
    fXB = gcs.fXB;
    fUbound=&fU[fXB.fCells];
    for (I=0;I<fXB.fCells;++I) {
        memcpy(thisup->fVector,fromup->fVector,sizeof(double)*fParm.fInNodes);
        thisup->fNc = fromup->fNc;
        for (J=0;J<thisup->fNc;++J) thisup->fC[J].fPtr = &fU[((TNeuralNetCell*)(fromup->fC[J].fPtr))->fID];
        thisup->fCount = fromup->fCount;
        thisup->fClass = fromup->fClass;
        ++thisup;
        ++fromup;
    }
}

int TGCS::GetWinnerCell(NNO_INTYPE* in) 
{
    int I,J;
    //get distance of all cells: performed in the learnstep before
    //compute deviation for all cells: dito
    
    //make output of all cells and get the winner
    double s_dist;
    double min_s_dist = DBL_MAX;
    TNeuralNetCell* up;
    J=0;
    for(up=fU;up<fUbound;++up) {
        double* v = up->fVector;
        double* d = up->fDiff;
        NNO_INTYPE* i = in;
        
        s_dist = 0.0;
        for (I=0;I<fParm.fInNodes;++I) {
            *d = *i++ - *v++;
            s_dist += *d * *d;
            ++d;
        }
        
        fOut[J++] = s_dist;
        
        if (s_dist<min_s_dist) {
            min_s_dist = s_dist;
            fUwin = up;
        }
    }
    
    if (fPlotter) fPlotter->AddTestSample(min_s_dist);
    
    return fUwin->fID;
}

double TGCS::Train(NNO_INTYPE* in,NNO_OUTTYPE*) 
{
    int I,J;
    Recall(in);  //make output of all cells and neurons find the winner
    TNeuralNetCell* unei;
    TNeuralNetCell* up;
    double* vwin = fUwin->fVector;
    double* dwin = fUwin->fDiff;
    for (J=0;J<fParm.fInNodes;++J) *vwin++ += *dwin++ * fXB.fWinStep;
    for (I=0;I<fUwin->fNc;++I) {
        unei = (TNeuralNetCell*)fUwin->fC[I].fPtr;
        double* v = unei->fVector;
        double* d = unei->fDiff;
        for (J=0;J<fParm.fInNodes;++J) *v++ += *d++ * fXB.fNeiStep;
    }
    
    ++fUwin->fCount;
    double m_win_count = 1.0 - fXB.fWinCount;
    for(up=fU;up<fUbound;++up) {
        up->fCount *= m_win_count;
    }
    
    if (fXB.fInsertStep>0)
        if (++fXB.fInsertCount==fXB.fInsertStep)
        {
            Insert();
            fXB.fInsertCount = 0;
        }
    
    if (fXB.fDeleteStep>0)
        if (++fXB.fDeleteCount==fXB.fDeleteStep) {
            Prune();
            fXB.fDeleteCount = 0;
        }
    
    fShouldSave = true;
    
    return fUwin->fID;
}

void TGCS::Deviation(void) 
{
    TNeuralNetCell* up;
    for(up=fU;up<fUbound;++up) TNeuralNetCell::GetSDev((TNeuralNetCell*)up,&fParm);
}


int TGCS::Insert(void) 
{
    int I,J;
    TNeuralNetCell* up = nullptr;
    TNeuralNetCell* umax = nullptr;
    TNeuralNetCell* udst = nullptr;
    TNeuralNetCell* unew = nullptr;
    
    if (fXB.fCells==fParm.fOutNodes) return 0; //break if there are no cells availiable
    //find cell with most win_count
    double win_count = -1;
    for(up=fU;up<fUbound;++up)
        if (up->fCount>win_count) {
            win_count = up->fCount;
            umax = up;
        }
    
    //create new cell
    unew = fUbound++;
    ++fXB.fCells;
    unew->fNc = 0;
    unew->fState = 0;
    
    //find cell with highest distance to MaxCount
    double s_dist;
    double max_s_dist = -1;
    TNeuralNetCell* unei;
    for (I=0;I<umax->fNc;++I) {
        s_dist = 0.0;
        unei=(TNeuralNetCell*)umax->fC[I].fPtr;
        for (J=0;J<fParm.fInNodes;++J) {
            double diff = umax->fVector[J] - unei->fVector[J];
            s_dist += diff*diff;
        }
        
        if (s_dist>max_s_dist) {
            max_s_dist = s_dist;
            udst = unei;
        }
    }
    
    TNeuralNetCell::ConnectNew((TNeuralNetCell*)unew,(TNeuralNetCell*)umax,(TNeuralNetCell*)udst,(TNeuralNetCellParameters*)&fXB);
    //    TNeuralNetCell::InitVector((TNeuralNetCell*)unew,&fParm);
    TNeuralNetCell::InitVector((TNeuralNetCell*)unew,(TNeuralNetCell*)umax,(TNeuralNetCell*)udst,&fParm);
    TNeuralNetCell::InitCount ((TNeuralNetCell*)unew);
    
    for(up=fU;up<fUbound;++up) TNeuralNetCell::CheckConnections((TNeuralNetCell*)up);
    return 1;  //insertion was successful
}

void TGCS::Remove(TNeuralNetCell* urem) 
{
    int I;
    TNeuralNetCell* ulast = --fUbound;
    --fXB.fCells;
    while(urem->fNc)  TNeuralNetCell::Disconnect((TNeuralNetCell*)urem->fC[0].fPtr,(TNeuralNetCell*)urem); //disconnect all neighbours of urem from urem    //copy ulast into urem
    memmove(urem->fVector,ulast->fVector,fParm.fInNodes*sizeof(double));
    memmove(urem->fC,ulast->fC,fXB.fConnectors*sizeof(connector));
    urem->fNc = ulast->fNc;
    urem->fCount = ulast->fCount;
    urem->fState = ulast->fState;
    while(ulast->fNc) TNeuralNetCell::Disconnect((TNeuralNetCell*)ulast->fC[0].fPtr,(TNeuralNetCell*)ulast); //disconnect all neighbours of ulast from ulast
    //Connect all former neigbours of ulast to urem
    for (I=0;I<urem->fNc;++I) {
        TNeuralNetCell* unei = (TNeuralNetCell*)urem->fC[I].fPtr;
        unei->fC[unei->fNc++].fPtr = urem;
    }
}

int TGCS::Prune(void) 
{
    int I,J;
    TNeuralNetCell* up =  nullptr;
    TNeuralNetCell* umax =  nullptr;;
    TNeuralNetCell* unei =  nullptr;;
    
    //find cell with maximal fChi2
    Deviation();
    double max_s_dev = -1;
    for(up=fU;up<fUbound;++up) {
        if (up->fChi2>max_s_dev) {
            max_s_dev = up->fChi2;
            umax=up;
        }
    }
    
    int RemoveCount = 1;
    umax->fState |= REMOVE;
    for (I=0;I<umax->fNc;++I) {
        unei = (TNeuralNetCell*)umax->fC[I].fPtr;
        int xnei = 0;
        for (J=0;J<unei->fNc;++J) {
            int II;
            TNeuralNetCell* uneinei = (TNeuralNetCell*)unei->fC[J].fPtr;
            if (uneinei!=umax) {
                for (II=0;II<umax->fNc;++II)
                    if (uneinei==(TNeuralNetCell*)umax->fC[II].fPtr) break;
                if (II==umax->fNc) ++xnei;
            }
        }
        
        //if unei has no neigbour outside, set REMOVE-flag
        if (xnei==0) {
            unei->fState |= REMOVE;
            ++RemoveCount;
        }
    }
    
    if (fXB.fCells-RemoveCount<fXB.fMinCells) {
        for(up=fU;up<fUbound;++up) up->fState&=~REMOVE; //clear REMOVE-flag on all cells if network would be shrinking to zero
    }
    else {
        up=fU;
        while(up<fUbound) {
            if (up->fState&REMOVE)
                Remove(up);
            else ++up;
        }
    }
    
    return 1;
}
