// TGNGTracker
//
// Implementation of the GROWING-NEURAL-GAS Trackfinder
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Heidelberg University, 1995

#include "TNtuple.h"

#include "TGNGTracker.h"
#include "VNeuralNetPlotter.h"
#include "Graph.h"

#include <cmath>
#include <cfloat>
#include <iostream>
using namespace std;

ClassImp(TGNGTracker)

TGNGTracker::TGNGTracker(int innodes,int maxCells,double winStep,double neiStep,
                         double aWinCount,double aEdgeCount,double minCount,int connectors,
                         Long_t insertStep,Long_t deleteStep,const char* netFile)
: VUnsupervisedNet("GNG",innodes,maxCells,netFile)
{
    fXB.fCells        = 2;
    fXB.fWinStep      = winStep;
    fXB.fNeiStep      = neiStep;
    fXB.fWinCount     = aWinCount;
    fXB.fEdgeCount    = aEdgeCount;
    fXB.fMinCount     = minCount;
    fXB.fConnectors   = connectors;
    fXB.fInsertStep   = insertStep;
    fXB.fDeleteStep   = deleteStep;
    fXB.fInsertCount  = 0;
    fXB.fDeleteCount  = 0;
    
    fXB.fMainEdgeCount = 1;
    fXB.fMainWinCount  = 1;
    
    fU = 0;
    AllocNet();
    InitNet();
}

// copy constructor
TGNGTracker::TGNGTracker(const TGNGTracker& gng,const char* netFile)
: VUnsupervisedNet("GNG",gng.fParm.fInNodes,gng.fParm.fOutNodes,netFile)
{
    fXB = gng.fXB;
    fU  = 0;
    AllocNet();
    InitNet();
    CopyData(gng);
}


TGNGTracker::~TGNGTracker()
{
    Deviation();
    //transfom fMainWinCount, fMainEdgeCount
    TNeuralNetCell* up;
    Long_t I;
    for(up=fU;up<fUbound;++up) {
        up->fCount*=fXB.fMainWinCount;
        for (I=0;I<up->fNc;++I) up->fAge[I]*=fXB.fMainEdgeCount;
    }
    fXB.fMainWinCount=1;
    fXB.fMainEdgeCount=1;
    if (fFilename!="") if (fShouldSave) Save();
    if (fU!=0) {
        up = fU;
        for (I=0;I<fParm.fOutNodes;++I) {
            delete[] up->fVector;
            delete[] up->fDiff;
            delete[] up->fC;
            delete[] up->fAge;
            ++up;
        }
        delete[] fU;
    }
}

void TGNGTracker::ReadBinary(void)
{
    fU = 0;
    fread(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    AllocNet();
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::ReadUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
        fread(up->fAge,sizeof(double),fXB.fConnectors,fFile);
        freadvar(up->fClass);
        for (I=0;I<up->fNc;++I) up->fC[I].fPtr = &fU[up->fC[I].fID];
    }
}

void TGNGTracker::ReadText(void)
{
    fU = 0;
    fscanf(fFile,"win_step     %le\n",&fXB.fWinStep);
    fscanf(fFile,"nei_step     %le\n",&fXB.fNeiStep);
    fscanf(fFile,"a_win_count  %le\n",&fXB.fWinCount);
    fscanf(fFile,"a_edge_count %le\n",&fXB.fEdgeCount);
    fscanf(fFile,"min_count    %le\n",&fXB.fMinCount);
    fscanf(fFile,"cells        %i\n",&fXB.fCells);
    fscanf(fFile,"connectors   %i\n",&fXB.fConnectors);
    fscanf(fFile,"insert_step  %li\n",&fXB.fInsertStep);
    fscanf(fFile,"delete_step  %li\n",&fXB.fDeleteStep);
    fscanf(fFile,"ins_count    %li\n",&fXB.fInsertCount);
    fscanf(fFile,"del_count    %li\n",&fXB.fDeleteCount);
    fscanf(fFile,"\n");
    fscanf(fFile,"main_win_count  %le\n",&fXB.fMainWinCount);
    fscanf(fFile,"main_edge_count %le\n",&fXB.fMainEdgeCount);
    
    AllocNet();
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::ReadUnitText(fFile,(TNeuralNetCell*)up,&fParm);
        fscanf(fFile,"\nedge count ");
        for (I=0;I<up->fNc;++I) fscanf(fFile,"%le ",&up->fAge[I]);
        fscanf(fFile,"\nclass            %i\n",&up->fClass);
        for (I=0;I<up->fNc;++I) up->fC[I].fPtr=&fU[up->fC[I].fID];
    }
}

void TGNGTracker::WriteBinary(void)
{
    TNeuralNetCell* up;
    fwrite(&fXB,sizeof(TNeuralNetCellParameters),1,fFile);
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::WriteUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
        fwrite(up->fAge,sizeof(double),fXB.fConnectors,fFile);
        fwritevar(up->fClass);
    }
}

void TGNGTracker::WriteText(void)
{
    fprintf(fFile,"win_step     %le\n",fXB.fWinStep);
    fprintf(fFile,"nei_step     %le\n",fXB.fNeiStep);
    fprintf(fFile,"a_win_count  %le\n",fXB.fWinCount);
    fprintf(fFile,"a_edge_count %le\n",fXB.fEdgeCount);
    fprintf(fFile,"min_count    %le\n",fXB.fMinCount);
    fprintf(fFile,"cells        %i\n",fXB.fCells);
    fprintf(fFile,"connectors   %i\n",fXB.fConnectors);
    fprintf(fFile,"insert_step  %li\n",fXB.fInsertStep);
    fprintf(fFile,"delete_step  %li\n",fXB.fDeleteStep);
    fprintf(fFile,"ins_count    %li\n",fXB.fInsertCount);
    fprintf(fFile,"del_count    %li\n",fXB.fDeleteCount);
    fprintf(fFile,"\n");
    fprintf(fFile,"main_win_count  %le\n",fXB.fMainWinCount);
    fprintf(fFile,"main_edge_count %le\n",fXB.fMainEdgeCount);
    
    TNeuralNetCell* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::WriteUnitText(fFile,(TNeuralNetCell*)up,&fParm);
        fprintf(fFile,"\nedge count ");
        for (I=0;I<up->fNc;++I) fprintf(fFile,"%le ",up->fAge[I]);
        fprintf(fFile,"\nclass            %i\n",up->fClass);
    }
}

void TGNGTracker::AllocNet(void)
{
    fU = new TNeuralNetCell[fParm.fOutNodes]; TestPointer(fU);
    fUbound=&fU[fXB.fCells];
    int I;
    TNeuralNetCell* up = fU;
    for (I=0;I<fParm.fOutNodes;++I) {
        up->fVector    = new double[fParm.fInNodes];  TestPointer(up->fVector);
        up->fDiff      = new double[fParm.fInNodes];  TestPointer(up->fDiff);
        up->fC         = new connector[fXB.fConnectors]; TestPointer(up->fC);
        up->fAge = new double[fXB.fConnectors];  TestPointer(up->fAge);
        up->fNc = 0;
        up->fCount = 0;
        up->fID = I;
        ++up;
    }
}

void TGNGTracker::InitNet(void)
{
    TNeuralNetCell* up;
    int J;
    for(up=fU;up<fUbound;++up) {
        for (J=0;J<fParm.fInNodes;++J) up->fVector[J]=Random();
        up->fNc = 0;
    }
    
    for (J=0;J<fParm.fOutNodes;++J) fU[J].fClass=0;
    
    Connect(&fU[0],&fU[1]);
}

void TGNGTracker::CopyData(const TGNGTracker& gng)
{
    TNeuralNetCell* thisup = fU;
    TNeuralNetCell* fromup = gng.fU;
    int I,J;
    
    //check integrity
    if (fParm.fInNodes    !=gng.fParm.fInNodes)     Errorf((char *)"cannot copy data; innodes not identical");
    if (fParm.fOutNodes   !=gng.fParm.fOutNodes)    Errorf((char *)"cannot copy data; outnodes not identical");
    if (fXB.fConnectors!=gng.fXB.fConnectors) Errorf((char *)"cannot copy data; max connectors not identical");
    
    fXB = gng.fXB;
    fUbound = &fU[fXB.fCells];
    for (I=0;I<fXB.fCells;++I) {
        memcpy(thisup->fVector,fromup->fVector,sizeof(double)*fParm.fInNodes);
        thisup->fNc = fromup->fNc;
        for (J=0;J<thisup->fNc;++J) {
            thisup->fC[J].fPtr=&fU[((TNeuralNetCell*)(fromup->fC[J].fPtr))->fID];
            thisup->fAge[J]=fromup->fAge[J];
        }
        thisup->fCount = fromup->fCount;
        thisup->fClass = fromup->fClass;
        ++thisup;
        ++fromup;
    }
}

//Disconnect only if both units have more than one connectors
int  TGNGTracker::CondDisconnect(TNeuralNetCell* up1,TNeuralNetCell* up2)
{
    int I;
    if ((up1->fNc==1) || (up2->fNc==1)) return 0;
    for (I=0;I<up1->fNc;++I) if ((TNeuralNetCell*)up1->fC[I].fPtr==up2) break;
    if (I<up1->fNc) up1->fC[I]=up1->fC[--up1->fNc]; else return 0;
    for (I=0;I<up2->fNc;++I) if ((TNeuralNetCell*)up2->fC[I].fPtr==up1) break;
    if (I<up2->fNc) up2->fC[I]=up2->fC[--up2->fNc];
    return 1;
}

void TGNGTracker::Connect(TNeuralNetCell* up1,TNeuralNetCell* up2)
{
    if ((up1->fNc==fXB.fConnectors)||(up2->fNc==fXB.fConnectors)) return;
    up1->fAge[up1->fNc]=1.0/fXB.fMainEdgeCount;
    up1->fC[up1->fNc++].fPtr = up2;
    up2->fAge[up2->fNc]=1.0/fXB.fMainEdgeCount;
    up2->fC[up2->fNc++].fPtr = up1;
}

void TGNGTracker::UpdateConnector(TNeuralNetCell* up1,TNeuralNetCell* up2)
{
    int I,J;
    for (I=0;I<up1->fNc;++I) if ((TNeuralNetCell*)up1->fC[I].fPtr==up2) break;
    if (I==up1->fNc) { Connect(up1,up2); return; }
    up1->fAge[I] += fXB.fEdgeCount;
    for (J=0;J<up2->fNc;++J) if ((TNeuralNetCell*)up2->fC[J].fPtr==up1) break;
    up2->fAge[J] = up1->fAge[I];
}

int  TGNGTracker::GetWinnerCell(NNO_INTYPE* in)
{
    int I,J;
    double s_dist;
    fMinDistSquare1 = DBL_MAX;
    fMinDistSquare2 = DBL_MAX;
    TNeuralNetCell* up;
    J = 0;
    for(up=fU;up<fUbound;++up){
        double* v = up->fVector;
        double* d = up->fDiff;
        NNO_INTYPE* i = in;
        s_dist = 0.0;
        for (I=0;I<fParm.fInNodes;++I) {
            *d =*i++ - *v++;
            s_dist+=*d * *d;
            ++d;
        }
        
        fOut[J++] = s_dist;
        
        if (s_dist<fMinDistSquare1) {
            fMinDistSquare2 = fMinDistSquare1;
            fUwin2 = fUwin1;
            fMinDistSquare1 = s_dist;
            fUwin1 = up;
        }
        else if (s_dist<fMinDistSquare2) {
            fMinDistSquare2 = s_dist;
            fUwin2=up;
        }
        
    }
    
    if (fPlotter) fPlotter->AddTestSample(fMinDistSquare1);
    
    return fUwin1->fID;
}

double  TGNGTracker::Train(NNO_INTYPE* in,NNO_OUTTYPE*)
{
    int I,J;
    Inference(in);  //make output of all cells and neurons find the winners
    TNeuralNetCell* unei;
    double* vwin = fUwin1->fVector;
    double* dwin = fUwin1->fDiff;
    for (J=0;J<fParm.fInNodes;++J) *vwin++ += *dwin++ * fXB.fWinStep;
    for (I=0;I<fUwin1->fNc;++I) {
        unei = (TNeuralNetCell*)fUwin1->fC[I].fPtr;
        double* v = unei->fVector;
        double* d = unei->fDiff;
        for (J=0;J<fParm.fInNodes;++J) *v++ += *d++ * fXB.fNeiStep;
    }
    UpdateConnector(fUwin1,fUwin2); //update edge_count of connector; if Uwin1,Uwin2 are not connected, connect them
    fUwin1->fCount += fMinDistSquare1;
    fXB.fMainWinCount  *= (1.0 - fXB.fWinCount);
    fXB.fMainEdgeCount *= (1.0 - fXB.fEdgeCount);
    
    if (fXB.fInsertStep>0)
        if (++fXB.fInsertCount==fXB.fInsertStep) {
            Insert();
            fXB.fInsertCount = 0;
        }
    
    if (fXB.fDeleteStep>0)
        if (++fXB.fDeleteCount==fXB.fDeleteStep) {
            Prune();
            fXB.fDeleteCount = 0;
        }
    
    fShouldSave = true;
    
    return fUwin1->fID;
}

void TGNGTracker::Deviation(void)
{
    TNeuralNetCell* up;
    for(up=fU;up<fUbound;++up) TNeuralNetCell::GetSDev((TNeuralNetCell*)up,&fParm);
}

int TGNGTracker::Insert(void)
{
    int I,J;
    TNeuralNetCell* up=0;
    TNeuralNetCell* umax1=0;
    TNeuralNetCell* umax2=0;
    TNeuralNetCell* unew=0;
    
    if (fXB.fCells==fParm.fOutNodes) return 0; //break if there are no cells availiable
    //find cell with highest win_count
    double win_count=-1;
    for(up=fU;up<fUbound;++up)
        if (up->fCount>win_count) {
            win_count = up->fCount;
            umax1 = up;
        }
    
    //create new cell
    unew = fUbound++;
    ++fXB.fCells;
    unew->fNc = 0;
    
    if (umax1->fNc==0) {
        Warningf(stdout,(char *)"CORRUPT NETWORK INTEGRITY! isolated cell found, please call developer");
        return 0;
    }
    
    //find neigbour to MaxCount with highest win_count
    win_count=-1;
    for (I=0;I<umax1->fNc;++I)
        if (((TNeuralNetCell*)umax1->fC[I].fPtr)->fCount>win_count) {
            win_count = ((TNeuralNetCell*)umax1->fC[I].fPtr)->fCount;
            umax2=(TNeuralNetCell*)umax1->fC[I].fPtr;
        }
    
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
    
    //transfom fMainWinCount, fMainEdgeCount
    for(up=fU;up<fUbound;++up) {
        up->fCount *= fXB.fMainWinCount;
        for (I=0;I<up->fNc;++I) up->fAge[I] *= fXB.fMainEdgeCount;
    }
    
    fXB.fMainWinCount  = 1;
    fXB.fMainEdgeCount = 1;
    
    for(up=fU;up<fUbound;++up) TNeuralNetCell::CheckConnections((TNeuralNetCell*)up);
    
    // Modification 1: Look for next data point and move cells to there
    
    if (fTuple == 0) {
        cout << "No ntuple set for training" << endl;
        return 0; // No training data available, leave cell where it is
    }
    
    double ux[3];
    ux[0] = unew->fVector[0];
    ux[1] = unew->fVector[1];
    ux[2] = unew->fVector[2];
    
    int imin = 0;
    double dmin = DBL_MAX;
    for (int i=0;i<fTuple->GetEntries();i++) {
        double dx[3];
        fTuple->GetEvent(i,1);
        float *x=fTuple->GetArgs();
        dx[0] = ux[0]-x[0];
        dx[1] = ux[1]-x[1];
        dx[2] = ux[2]-x[2];
        double d = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        if (d < dmin) {
            imin = i;
            dmin = d;
        }
    }
    
    fTuple->GetEvent(imin,1);
    float *x=fTuple->GetArgs();
    unew->fVector[0] = x[0];
    unew->fVector[1] = x[1];
    unew->fVector[2] = x[2];
    
    return 1;  //insertion was successful
}

void TGNGTracker::Prune(void)
{ // remove edges
    int I;
    TNeuralNetCell* up;
    
    //transfom fMainWinCount, fMainEdgeCount
    for(up=fU;up<fUbound;++up) {
        up->fCount*=fXB.fMainWinCount;
        for (I=0;I<up->fNc;++I) up->fAge[I]*=fXB.fMainEdgeCount;
    }
    fXB.fMainWinCount  = 1;
    fXB.fMainEdgeCount = 1;
    
    //remove all edges with edge_count<fMinCount
    for(up=fU;up<fUbound;++up) {
        I = 0;
        while (I<up->fNc) {
            // Modification 2: Favor straight lines
            const TNeuralNetCell *c = GetCell(I);
            unsigned int numberConnections = c->GetNumberOfConnections();
            if (up!=c) { // Make sure we have different cells
                const double *x0 = up->GetVector();
                const double *xi = c->GetVector();
                double x1[3]; // vector to neighbour
                x1[0] = xi[0] - x0[0];
                x1[1] = xi[1] - x0[1];
                x1[2] = xi[2] - x0[2];
                double rx1 = sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]); // Disctance between cells
                for (long j=0;j<numberConnections;j++) {
                    const TNeuralNetCell *u = c->GetConnectedCell(j);
                    const double *xj = u->GetVector();
                    double x2[3];
                    x2[0] = xj[0] - x0[0];
                    x2[1] = xj[1] - x0[1];
                    x2[2] = xj[2] - x0[2];
                    double rx2 = sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);
                    double cosine = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
                    double penalty = 1. - cosine/(rx1*rx2); // punish small opening angles
                    up->fAge[I] *= penalty;
                    up->fAge[j] *= penalty;
                }
            }
            
            if (up->fAge[I]<fXB.fMinCount)
                if (CondDisconnect(up,(TNeuralNetCell*)up->fC[I].fPtr)) --I; // Remove the connection
            
            ++I;
        }
        
    }
    
}

int TGNGTracker::GetCluster(void) {

    Graph g(GetNumberOfCells());

    TNeuralNetCell* up;
    for(up=fU;up<fUbound;++up) {
        const TNeuralNetCell c = *up;
 //       cout << "Cell " << c.GetID() << " Connections:" << c.GetNumberOfConnections() << endl;
        for (int I=0; I<c.GetNumberOfConnections(); I++) {
            const TNeuralNetCell *cc = c.GetConnectedCell(I);
//            cout << cc->GetID() << " ";
            g.addEdge(c.GetID(),cc->GetID());
        }
//        cout << endl;
    }

    cout << "Following are strongly connected components in "
    "given graph \n";
    int n = g.printSCCs();

    return n;
}
