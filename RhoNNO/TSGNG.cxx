// TSGNG
//
// Implementation of the SUPERVISED-GROWING-NEURAL-GAS (SGNG)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TSGNG.h"
#include "VNeuralNetPlotter.h"

#include <cfloat>
#include <string>
using namespace std;

ClassImp(TSGNG)

TSGNG::TSGNG(int innodes,
             int outnodes,
             int maxCells,
             double winStep,
             double neiStep,
             double neuStep,
             double aErrCount,
             double aEdgeCount,
             double bSDev,
             double minCount,
             int  connectors,
             long insertStep,
             long deleteStep,
             string netFile)
: VSupervisedNet("SGNG",innodes,outnodes,netFile)
{
    fXB.fWinStep     = winStep;
    fXB.fNeiStep     = neiStep;
    fXB.fNeuStep     = neuStep;
    fXB.fErrCount    = aErrCount;
    fXB.fEdgeCount   = aEdgeCount;
    fXB.fNeiCount    = bSDev;
    fXB.fMinCount    = minCount;
    fXB.fCells	     = 2;
    fXB.fMaxCells    = maxCells;
    fXB.fConnectors  = connectors;
    fXB.fInsertStep  = insertStep;
    fXB.fDeleteStep  = deleteStep;
    fXB.fInsertCount = 0;
    fXB.fDeleteCount = 0;
    
    fXB.fMainEdgeCount	= 1;
    fXB.fMainErrCount	= 1;
    
    fU = 0;
    AllocNet();
    InitNet();
}

// copy constructor
TSGNG::TSGNG(const TSGNG& sgng,string netFile)
: VSupervisedNet("SGNG",sgng.fParm.fInNodes,sgng.fParm.fOutNodes,netFile) 
{
    fXB = sgng.fXB;
    fU  = 0;
    AllocNet();
    InitNet();
    CopyData(sgng);
}

TSGNG::~TSGNG() 
{
    //transform main_err_count, main_edge_count
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        up->fCount*=fXB.fMainErrCount;
        for (I=0;I<up->fNc;++I) up->fAge[I]*=fXB.fMainEdgeCount;
    }
    fXB.fMainErrCount=1;
    fXB.fMainEdgeCount=1;
    if (fFilename!="") if (fShouldSave) Save();
    if (fU!=0) {
        up = fU;
        for (I=0;I<fXB.fMaxCells;++I) {
            delete[] up->fVector;
            delete[] up->fWeight;
            delete[] up->fDiff;
            delete[] up->fC;
            delete[] up->fAge;
            ++up;
        }
        delete[] fU;
    }
}

void TSGNG::ReadBinary(void) 
{
    fU = 0;
    fread(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    AllocNet();
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::ReadUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
        fread(up->fAge,sizeof(double),up->fNc,fFile);
        fread(up->fWeight,sizeof(double),fParm.fOutNodes,fFile);
        for (I=0;I<up->fNc;++I) up->fC[I].fPtr = &fU[up->fC[I].fID];
    }
}

void TSGNG::ReadText(void) 
{
    fU = 0;
    fscanf(fFile,"win_step     %le\n",&fXB.fWinStep);
    fscanf(fFile,"nei_step     %le\n",&fXB.fNeiStep);
    fscanf(fFile,"neu_step     %le\n",&fXB.fNeuStep);
    fscanf(fFile,"a_err_count  %le\n",&fXB.fErrCount);
    fscanf(fFile,"a_edge_count %le\n",&fXB.fEdgeCount);
    fscanf(fFile,"b_s_dev      %le\n",&fXB.fNeiCount);
    fscanf(fFile,"min_count    %le\n",&fXB.fMinCount);
    fscanf(fFile,"cells        %i\n",&fXB.fCells);
    fscanf(fFile,"max_cells    %i\n",&fXB.fMaxCells);
    fscanf(fFile,"connectors   %i\n",&fXB.fConnectors);
    fscanf(fFile,"insert_step  %li\n",&fXB.fInsertStep);
    fscanf(fFile,"delete_step  %li\n",&fXB.fDeleteStep);
    fscanf(fFile,"ins_count    %li\n",&fXB.fInsertCount);
    fscanf(fFile,"del_count    %li\n",&fXB.fDeleteCount);
    fscanf(fFile,"\n");
    fscanf(fFile,"main_err_count  %le\n",&fXB.fMainErrCount);
    fscanf(fFile,"main_edge_count %le\n",&fXB.fMainEdgeCount);
    
    AllocNet();
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::ReadUnitText(fFile,(TNeuralNetCell*)up,&fParm);
        fscanf(fFile,"\nedge count ");
        for (I=0;I<up->fNc;++I) fscanf(fFile,"%le ",&up->fAge[I]);
        fscanf(fFile,"\nweights ");
        for (I=0;I<fParm.fOutNodes;++I) fscanf(fFile,"%le ",&up->fWeight[I]);
        fscanf(fFile,"\n");
        for (I=0;I<up->fNc;++I) up->fC[I].fPtr=&fU[up->fC[I].fID];
    }
}

void TSGNG::WriteBinary(void) 
{
    TNeuralNetCell* up;
    fwrite(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::WriteUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
        fwrite(up->fAge,sizeof(double),up->fNc,fFile);
        fwrite(up->fWeight,sizeof(double),fParm.fOutNodes,fFile);
    }
}

void TSGNG::WriteText(void) 
{
    fprintf(fFile,"win_step     %le\n",fXB.fWinStep);
    fprintf(fFile,"nei_step     %le\n",fXB.fNeiStep);
    fprintf(fFile,"neu_step     %le\n",fXB.fNeuStep);
    fprintf(fFile,"a_err_count  %le\n",fXB.fErrCount);
    fprintf(fFile,"a_edge_count %le\n",fXB.fEdgeCount);
    fprintf(fFile,"b_s_dev      %le\n",fXB.fNeiCount);
    fprintf(fFile,"min_count    %le\n",fXB.fMinCount);
    fprintf(fFile,"cells        %i\n",fXB.fCells);
    fprintf(fFile,"max_cells    %i\n",fXB.fMaxCells);
    fprintf(fFile,"connectors   %i\n",fXB.fConnectors);
    fprintf(fFile,"insert_step  %li\n",fXB.fInsertStep);
    fprintf(fFile,"delete_step  %li\n",fXB.fDeleteStep);
    fprintf(fFile,"ins_count    %li\n",fXB.fInsertCount);
    fprintf(fFile,"del_count    %li\n",fXB.fDeleteCount);
    fprintf(fFile,"\n");
    fprintf(fFile,"main_err_count  %le\n",fXB.fMainErrCount);
    fprintf(fFile,"main_edge_count %le\n",fXB.fMainEdgeCount);
    
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::WriteUnitText(fFile,(TNeuralNetCell*)up,&fParm);
        fprintf(fFile,"\nedge count ");
        for (I=0;I<up->fNc;++I) fprintf(fFile,"%le ",up->fAge[I]);
        fprintf(fFile,"\nweights ");
        for (I=0;I<fParm.fOutNodes;++I) fprintf(fFile,"%le ",up->fWeight[I]);
        fprintf(fFile,"\n");
    }
}

void TSGNG::AllocNet(void) 
{
    int I;
    fU = new TNeuralNetCell[fXB.fMaxCells];  TestPointer(fU);
    fUbound = &fU[fXB.fCells];
    TNeuralNetCell* up = fU;
    for (I=0;I<fXB.fMaxCells;++I) {
        up->fVector = new double[fParm.fInNodes];    TestPointer(up->fVector);
        up->fWeight = new double[fParm.fOutNodes];    TestPointer(up->fWeight);
        up->fDiff   = new double[fParm.fInNodes];    TestPointer(up->fDiff);
        up->fC	    = new connector[fXB.fConnectors]; TestPointer(up->fC);
        up->fAge    = new double[fXB.fConnectors];    TestPointer(up->fAge);
        up->fNc     = 0;
        up->fChi2   = 0;
        up->fCount  = 0;
        up->fOut    = 0;
        up->fID     = I;
        ++up;
    }
}

void TSGNG::InitNet(void) 
{
    TNeuralNetCell* up;
    int I,J;
    for(up=fU;up<fUbound;++up) {
        for (J=0;J<fParm.fInNodes;++J) up->fVector[J]=Random();
        for (J=0;J<fParm.fOutNodes;++J) up->fWeight[J]=Random();
        up->fNc=0;
    }
    Connect(&fU[0],&fU[1]);
    for(up=fU;up<fUbound;++up) {
        double s_dist;
        up->fChi2 = 0;
        for (J=0;J<up->fNc;++J) {
            TNeuralNetCell* unei = (TNeuralNetCell*)up->fC[J].fPtr;
            s_dist = 0;
            for (I=0;I<fParm.fInNodes;++I) {double diff=up->fVector[I]-unei->fVector[I]; s_dist+=diff*diff;}
            up->fChi2 += s_dist;
        }
        up->fChi2 /= up->fNc;
    }
}

void TSGNG::CopyData(const TSGNG& sgng) 
{
    TNeuralNetCell* thisup = fU;
    TNeuralNetCell* fromup = sgng.fU;
    int I,J;
    
    //check integrity
    if (fParm.fInNodes    !=sgng.fParm.fInNodes)     Errorf((char *)"cannot copy data; innodes not identical");
    if (fParm.fOutNodes   !=sgng.fParm.fOutNodes)    Errorf((char *)"cannot copy data; outnodes not identical");
    if (fXB.fConnectors!=sgng.fXB.fConnectors) Errorf((char *)"cannot copy data; max connectors not identical");
    if (fXB.fMaxCells  !=sgng.fXB.fMaxCells )  Errorf((char *)"cannot copy data; max_cells not identical");
    
    fXB=sgng.fXB;
    fUbound=&fU[fXB.fCells];
    for (I=0;I<fXB.fCells;++I) {
        memcpy(thisup->fVector,fromup->fVector,sizeof(double)*fParm.fInNodes);
        memcpy(thisup->fWeight,fromup->fWeight,sizeof(double)*fParm.fOutNodes);
        thisup->fNc=fromup->fNc;
        for (J=0;J<thisup->fNc;++J) {
            thisup->fC[J].fPtr = &fU[((TNeuralNetCell*)(fromup->fC[J].fPtr))->fID];
            thisup->fAge[J] = fromup->fAge[J];
        }
        thisup->fCount = fromup->fCount;
        ++thisup;
        ++fromup;
    }
}

//Disconnect only if both TNeuralNetCells have more than one connectors
int  TSGNG::CondDisconnect(TNeuralNetCell* up1,TNeuralNetCell* up2) 
{
    int I;
    if ((up1->fNc==1) || (up2->fNc==1)) return 0;
    for (I=0;I<up1->fNc;++I) if ((TNeuralNetCell*)up1->fC[I].fPtr==up2) break;
    if (I<up1->fNc) up1->fC[I]=up1->fC[--up1->fNc]; else return 0;
    for (I=0;I<up2->fNc;++I) if ((TNeuralNetCell*)up2->fC[I].fPtr==up1) break;
    if (I<up2->fNc) up2->fC[I]=up2->fC[--up2->fNc];
    return 1;
}

void TSGNG::Connect(TNeuralNetCell* up1,TNeuralNetCell* up2) 
{
    if ((up1->fNc==fXB.fConnectors)||(up2->fNc==fXB.fConnectors)) return;
    up1->fAge[up1->fNc]=1.0/fXB.fMainEdgeCount;
    up1->fC[up1->fNc++].fPtr=up2;
    up2->fAge[up2->fNc]=1.0/fXB.fMainEdgeCount;
    up2->fC[up2->fNc++].fPtr=up1;
}

void TSGNG::UpdateConnector(TNeuralNetCell* up1,TNeuralNetCell* up2) 
{
    int I,J;
    for (I=0;I<up1->fNc;++I) if ((TNeuralNetCell*)up1->fC[I].fPtr==up2) break;
    if (I==up1->fNc) { Connect(up1,up2); return; }
    up1->fAge[I]+=fXB.fEdgeCount;
    for (J=0;J<up2->fNc;++J) if ((TNeuralNetCell*)up2->fC[J].fPtr==up1) break;
    up2->fAge[J]=up1->fAge[I];
}

double* TSGNG::Recall(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    int I,J;
    double s_dist;
    double sum_out = 0;
    fMinDistSquare1 = DBL_MAX;
    fMinDistSquare2 = DBL_MAX;
    TNeuralNetCell* up;
    
    for(up=fU;up<fUbound;++up){
        double* v = up->fVector;
        double* d = up->fDiff;
        NNO_INTYPE* i = in;
        s_dist = 0;
        for (I=0;I<fParm.fInNodes;++I) { *d =*i++ - *v++; s_dist+=*d * *d; ++d; }
        up->fOut = exp(-s_dist/up->fChi2);
        sum_out += up->fOut;
        if (s_dist<fMinDistSquare1) {
            fMinDistSquare2 = fMinDistSquare1;
            fUwin2=fUwin1;
            fMinDistSquare1=s_dist;
            fUwin1=up;
        }
        else if (s_dist<fMinDistSquare2) {
            fMinDistSquare2 = s_dist;
            fUwin2=up;
        }
    }
    
    for (J=0;J<fParm.fOutNodes;++J) fOut[J] = 0;
    
    if (sum_out>0) {
        for(up=fU;up<fUbound;++up) {
            up->fOut /= sum_out;
            for (J=0;J<fParm.fOutNodes;++J)
                fOut[J] += up->fOut * up->fWeight[J];
        }
    }
    
    if (fPlotter) {
        bool good = true;
        if (out!=0) good = out[0]>fParm.fThreshold;
        fPlotter->AddTestSample(fOut[0],good);
    }
    
    return fOut;
}


double TSGNG::Train(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    int I,J;
    double S_Err;
    
    double s_dist;
    double sum_out = 0;
    fMinDistSquare1 = DBL_MAX;
    fMinDistSquare2 = DBL_MAX;
    TNeuralNetCell* up;
    
    // Recall
    for(up=fU;up<fUbound;++up){
        double* v = up->fVector;
        double* d = up->fDiff;
        NNO_INTYPE* i = in;
        s_dist = 0;
        for (I=0;I<fParm.fInNodes;++I) { *d =*i++ - *v++; s_dist+=*d * *d; ++d; }
        up->fOut = exp(-s_dist/up->fChi2);
        sum_out += up->fOut;
        if (s_dist<fMinDistSquare1) {
            fMinDistSquare2 = fMinDistSquare1;
            fUwin2=fUwin1;
            fMinDistSquare1=s_dist;
            fUwin1=up;
        }
        else if (s_dist<fMinDistSquare2) {
            fMinDistSquare2 = s_dist;
            fUwin2=up;
        }
    }
    
    for (J=0;J<fParm.fOutNodes;++J) fOut[J] = 0;
    
    if (sum_out>0) {
        for(up=fU;up<fUbound;++up) {
            up->fOut /= sum_out;
            for (J=0;J<fParm.fOutNodes;++J)
                fOut[J] += up->fOut * up->fWeight[J];
        }
    }
    
    // Weights update
    UpdateConnector(fUwin1,fUwin2); //update edge_count of connector; if fUwin1,fUwin2 are not connected, connect them
    TNeuralNetCell* unei;
    double* vwin = fUwin1->fVector;
    double* dwin = fUwin1->fDiff;
    for (J=0;J<fParm.fInNodes;++J) *vwin++ += *dwin++ * fXB.fWinStep;
    fUwin1->fChi2 = 0;
    for (I=0;I<fUwin1->fNc;++I) {
        unei = (TNeuralNetCell*)fUwin1->fC[I].fPtr;
        double* v = unei->fVector;
        double* d = unei->fDiff;
        s_dist = 0.0;
        double diff;
        vwin = fUwin1->fVector;
        for (J=0;J<fParm.fInNodes;++J) {*v += *d++ * fXB.fNeiStep; diff=*vwin++ - *v++; s_dist+=diff*diff;}
        fUwin1->fChi2 += s_dist;
        unei->fChi2 += (s_dist-unei->fChi2)*fXB.fNeiCount;
    }
    
    fUwin1->fChi2 /= fUwin1->fNc;
    
    S_Err=0;
    for (J=0;J<fParm.fOutNodes;++J) {
        double diff = out[J] - fOut[J];
        S_Err+=diff*diff;
        for(up=fU;up<fUbound;++up) up->fWeight[J] += diff * up->fOut * fXB.fNeuStep;
    }
    
    fUwin1->fCount += S_Err;
    fXB.fMainErrCount *=(1.0-fXB.fErrCount);
    fXB.fMainEdgeCount*=(1.0-fXB.fEdgeCount);
    
    if (fXB.fInsertStep>0)
        if (fXB.fInsertCount++==fXB.fInsertStep) {
            Insert();
            fXB.fInsertCount=0;
        }
    
    if (fXB.fDeleteStep>0)
        if (fXB.fDeleteCount++==fXB.fDeleteStep) {
            Prune();
            fXB.fDeleteCount=0;
        }
    
    fShouldSave = true;
    
    if (fPlotter) fPlotter->AddTrainSample(out[0],out[0]>fParm.fThreshold);
    
    return S_Err;
}

int TSGNG::Insert(void) 
{
    int I,J;
    TNeuralNetCell* up = nullptr;
    TNeuralNetCell* umax1 = nullptr;
    TNeuralNetCell* umax2 = nullptr;
    TNeuralNetCell* unew = nullptr;
    if (fXB.fCells==fXB.fMaxCells) return 0; //break if there are no cells availiable
    
    //find cell with highest err_count
    double err_count = -1;
    for(up=fU;up<fUbound;++up) if (up->fCount>err_count) {err_count=up->fCount; umax1=up;}
    
    //create new cell
    unew = fUbound++;
    ++fXB.fCells;
    unew->fNc = 0;
    
    if (umax1->fNc==0) { Warningf(stdout,(char *)"CORRUPT NETWORK INTEGRITY! isolated cell found, please call developer"); return 0; }
    
    //find neighbour with highest err_count
    err_count=-1;
    for (I=0;I<umax1->fNc;++I) if (((TNeuralNetCell*)umax1->fC[I].fPtr)->fCount>err_count) { err_count=((TNeuralNetCell*)umax1->fC[I].fPtr)->fCount; umax2=(TNeuralNetCell*)umax1->fC[I].fPtr; }
    
    //connect new cell with common neighbours of umax1 and umax2
    for (I=0;I<umax1->fNc;++I)
        for (J=0;J<umax2->fNc;++J)
            if (umax1->fC[I].fPtr==umax2->fC[J].fPtr)
                Connect((TNeuralNetCell*)umax1->fC[I].fPtr,unew);
    
    Connect(unew,umax1);     //connect umax1 and unew
    Connect(unew,umax2);     //connect umax2 and unew
    TNeuralNetCell::Disconnect((TNeuralNetCell*)umax1,(TNeuralNetCell*)umax2);  //disconnect umax1 and umax2
    
    TNeuralNetCell::InitVector((TNeuralNetCell*)unew,(TNeuralNetCell*)umax1,(TNeuralNetCell*)umax2,&fParm);
    TNeuralNetCell::InitCount ((TNeuralNetCell*)unew);
    TNeuralNetCell::InitWgt   ((TNeuralNetCell*)unew,&fParm);
    TNeuralNetCell::GetSDev   ((TNeuralNetCell*)unew,&fParm);
    
    //transfom main_err_count, main_edge_count
    for (up=fU;up<fUbound;++up) {
        up->fCount*=fXB.fMainErrCount;
        for (I=0;I<up->fNc;++I) up->fAge[I] *= fXB.fMainEdgeCount;
    }
    fXB.fMainErrCount=1;
    fXB.fMainEdgeCount=1;
    
    for(up=fU;up<fUbound;++up) TNeuralNetCell::CheckConnections((TNeuralNetCell*)up);
    
    return 1;  //insertion was successful
}

void TSGNG::Prune(void) 
{ // remove edges
    int I;
    TNeuralNetCell* up;
    
    //transfom main_err_count, main_edge_count
    for(up=fU;up<fUbound;++up) {
        up->fCount*=fXB.fMainErrCount;
        for (I=0;I<up->fNc;++I) up->fAge[I] *= fXB.fMainEdgeCount;
    }
    fXB.fMainErrCount  = 1;
    fXB.fMainEdgeCount = 1;
    
    //remove all edges with edge_count<min_count
    for(up=fU;up<fUbound;++up) {
        I=0;
        while (I<up->fNc) {
            if (up->fAge[I]<fXB.fMinCount) {
                if (!CondDisconnect(up,(TNeuralNetCell*)up->fC[I].fPtr)) ++I;
            } else ++I;
        }
    }
}

