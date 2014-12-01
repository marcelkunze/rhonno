//////////////////////////////////////////////////////////////////////////
//									//
// TSGCS								//
//									//
// Implementation of the SUPERVISED-GROWING-CELL-STRUCTURE (SGCS)	//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "RhoNNO/TSGCS.h"
#include "RhoNNO/VNeuralNetPlotter.h"

ClassImp(TSGCS)

//cell states
#define REMOVE 1

TSGCS::TSGCS(Int_t innodes,
	     Int_t cells,
	     Int_t outnodes,
	     Int_t maxCells,
	     Double_t winStep,
	     Double_t neiStep,
	     Double_t neuStep,
	     Double_t aErrCount,
	     Double_t bSDev,
	     Int_t connectors,
	     Long_t insertStep,
	     Long_t deleteStep,
	     const char* netFile)
	     : VSupervisedNet("SGCS",innodes,outnodes,netFile) 
{
    fXB.fWinStep	= winStep;
    fXB.fNeiStep	= neiStep;
    fXB.fNeuStep	= neuStep;
    fXB.fErrCount	= aErrCount;
    fXB.fNeiCount	= bSDev;
    fXB.fCells	= cells;
    fXB.fMinCells	= cells;
    fXB.fMaxCells	= maxCells;
    fXB.fConnectors	= connectors;
    fXB.fInsertStep	= insertStep;
    fXB.fDeleteStep	= deleteStep;
    fXB.fInsertCount	= 0;
    fXB.fDeleteCount	= 0;
    fU = 0;
    AllocNet();
    InitNet();
}

// copy constructor
TSGCS::TSGCS(const TSGCS& sgcs,const char* netFile)
: VSupervisedNet("SGCS",sgcs.fParm.fInNodes,sgcs.fParm.fOutNodes,netFile) 
{
    fXB = sgcs.fXB;
    fU  = 0;
    AllocNet();
    InitNet();
    CopyData(sgcs);
}

TSGCS::~TSGCS() 
{
    if (fFilename!="") if (fShouldSave) Save();
    if (fU!=0) {
	TNeuralNetCell* up = fU;
	Int_t I;
	for (I=0;I<fXB.fMaxCells;++I) {
	    delete[] up->fVector;
	    delete[] up->fWeight;
	    delete[] up->fDiff;
	    delete[] up->fC;
	    ++up;
	}
	delete[] fU;
    }
}

void TSGCS::ReadBinary(void) 
{
    fU = 0;
    fread(&fXB,sizeof(TNeuralNetCell),1,fFile);
    AllocNet();
    TNeuralNetCell* up;
    Int_t I;
    for(up=fU;up<fUbound;++up) {
        TNeuralNetCell::ReadUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
	fread(up->fWeight,sizeof(Double_t),fParm.fOutNodes,fFile);
	for (I=0;I<up->fNc;++I) up->fC[I].fPtr=&fU[up->fC[I].fID];
    }
}

void TSGCS::ReadText(void) 
{
    fU = 0;
    fscanf(fFile,"win_step    %le\n",&fXB.fWinStep);
    fscanf(fFile,"nei_step    %le\n",&fXB.fNeiStep);
    fscanf(fFile,"neu_step    %le\n",&fXB.fNeuStep);
    fscanf(fFile,"a_err_count %le\n",&fXB.fErrCount);
    fscanf(fFile,"b_s_dev     %le\n",&fXB.fNeiCount);
    fscanf(fFile,"cells       %i\n",&fXB.fCells);
    fscanf(fFile,"min_cells   %i\n",&fXB.fMinCells);
    fscanf(fFile,"max_cells   %i\n",&fXB.fMaxCells);
    fscanf(fFile,"connectors  %i\n",&fXB.fConnectors);
    fscanf(fFile,"insert_step %li\n",&fXB.fInsertStep);
    fscanf(fFile,"delete_step %li\n",&fXB.fDeleteStep);
    fscanf(fFile,"ins_count   %li\n",&fXB.fInsertCount);
    fscanf(fFile,"del_count   %li\n",&fXB.fDeleteCount);
    
    AllocNet();
    TNeuralNetCell* up;
    Int_t I;
    for(up=fU;up<fUbound;++up) {
	TNeuralNetCell::ReadUnitText(fFile,(TNeuralNetCell*)up,&fParm);
	fscanf(fFile,"\nweights\n");
	for (I=0;I<fParm.fOutNodes;++I) fscanf(fFile,"%le",&up->fWeight[I]);
	fscanf(fFile,"\n");
	for (I=0;I<up->fNc;++I) up->fC[I].fPtr=&fU[up->fC[I].fID];
    }
}

void TSGCS::WriteBinary(void) 
{
    TNeuralNetCell* up;
    fwrite(&fXB,sizeof(TNeuralNetCell),1,fFile);
    for(up=fU;up<fUbound;++up) {
	TNeuralNetCell::WriteUnitBinary(fFile,(TNeuralNetCell*)up,&fParm);
	fwrite(up->fWeight,sizeof(Double_t),fParm.fOutNodes,fFile);
    }
}

void TSGCS::WriteText(void) 
{
    fprintf(fFile,"win_step    %le\n",fXB.fWinStep);
    fprintf(fFile,"nei_step    %le\n",fXB.fNeiStep);
    fprintf(fFile,"neu_step    %le\n",fXB.fNeuStep);
    fprintf(fFile,"a_err_count %le\n",fXB.fErrCount);
    fprintf(fFile,"b_s_dev     %le\n",fXB.fNeiCount);
    fprintf(fFile,"cells       %i\n",fXB.fCells);
    fprintf(fFile,"min_cells   %i\n",fXB.fMinCells);
    fprintf(fFile,"max_cells   %i\n",fXB.fMaxCells);
    fprintf(fFile,"connectors  %i\n",fXB.fConnectors);
    fprintf(fFile,"insert_step %li\n",fXB.fInsertStep);
    fprintf(fFile,"delete_step %li\n",fXB.fDeleteStep);
    fprintf(fFile,"ins_count   %li\n",fXB.fInsertCount);
    fprintf(fFile,"del_count   %li\n",fXB.fDeleteCount);
    
    TNeuralNetCell* up;
    Int_t I;
    for(up=fU;up<fUbound;++up) {
	TNeuralNetCell::WriteUnitText(fFile,(TNeuralNetCell*)up,&fParm);
	fprintf(fFile,"\nweights\n");
	for (I=0;I<fParm.fOutNodes;++I) fprintf(fFile,"%le\n",up->fWeight[I]);
	fprintf(fFile,"\n");
    }
}

void TSGCS::AllocNet(void) 
{
    Int_t I;
    fU = new TNeuralNetCell[fXB.fMaxCells]; TestPointer(fU);
    fUbound = &fU[fXB.fCells];
    TNeuralNetCell* up = fU;
    for (I=0;I<fXB.fMaxCells;++I) {
	up->fVector	= new Double_t[fParm.fInNodes];  TestPointer(up->fVector);
	up->fWeight	= new Double_t[fParm.fOutNodes]; TestPointer(up->fWeight);
	up->fDiff	= new Double_t[fParm.fInNodes];  TestPointer(up->fDiff);
	up->fC		= new connector[fXB.fConnectors]; TestPointer(up->fC);
	up->fNc		= 0;
	up->fChi2	= 0;
	up->fCount	= 0;
	up->fOut	= 0;
	up->fID		= I;
	up->fState	= 0;
	++up;
    }
}

void TSGCS::InitNet(void) 
{
    TNeuralNetCell* up;
    Int_t I,J;
    for(up=fU;up<fUbound;++up) {
	for (J=0;J<fParm.fInNodes;++J) up->fVector[J]=Random();
	for (J=0;J<fParm.fOutNodes;++J) up->fWeight[J]=Random();
	up->fNc=0;
	for (J=0;J<fXB.fCells;++J) 
	    if (&fU[J]!=up) up->fC[up->fNc++].fPtr = &fU[J];
    }
    for(up=fU;up<fUbound;++up) {
	Double_t s_dist;
	up->fChi2 = 0.0;
	for (J=0;J<up->fNc;++J) {
	    TNeuralNetCell* unei=(TNeuralNetCell*)up->fC[J].fPtr;
	    s_dist = 0.0;
	    for (I=0;I<fParm.fInNodes;++I) { 
		Double_t diff = up->fVector[I]-unei->fVector[I]; 
		s_dist += diff*diff;
	    }

	    up->fChi2 += s_dist;
	}

	up->fChi2 /= up->fNc;
    }
}

void TSGCS::CopyData(const TSGCS& sgcs) 
{
    TNeuralNetCell* thisup = fU;
    TNeuralNetCell* fromup = sgcs.fU;
    Int_t I,J;
    
    //check integrity
    if (fParm.fInNodes    !=sgcs.fParm.fInNodes)     Errorf((char *)"cannot copy data; innodes not identical");
    if (fParm.fOutNodes   !=sgcs.fParm.fOutNodes)    Errorf((char *)"cannot copy data; outnodes not identical");
    if (fXB.fConnectors< sgcs.fXB.fConnectors) Errorf((char *)"cannot copy data; max connectors too low");
    if (fXB.fMaxCells  < sgcs.fXB.fCells )     Errorf((char *)"cannot copy data; not enough cells availiable");
    
    fXB = sgcs.fXB;
    fUbound = &fU[fXB.fCells];
    for (I=0;I<fXB.fCells;++I) {
	memcpy(thisup->fVector,fromup->fVector,sizeof(Double_t)*fParm.fInNodes);
	memcpy(thisup->fWeight,fromup->fWeight,sizeof(Double_t)*fParm.fOutNodes);
	thisup->fNc = fromup->fNc;
	for (J=0;J<thisup->fNc;++J) 
	    thisup->fC[J].fPtr = &fU[((TNeuralNetCell*)(fromup->fC[J].fPtr))->fID];
	thisup->fCount = fromup->fCount;
	++thisup;
	++fromup;
    }
}

Double_t* TSGCS::Recall(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    Int_t I,J;
    //get distance of all cells: performed in the learnstep before
    //compute deviation for all cells: dito
    
    //make output of all cells and get the winner
    Double_t s_dist;
    Double_t sum_out = 0.0;
    Double_t min_s_dist = DBL_MAX;
    TNeuralNetCell* up;
    
    for(up=fU;up<fUbound;++up){
	Double_t* v = up->fVector;
	Double_t* d = up->fDiff;
	NNO_INTYPE* i = in;
	s_dist = 0.0;
	for (I=0;I<fParm.fInNodes;++I) { 
	    *d = *i++ - *v++; 
	    s_dist += *d * *d; 
	    ++d;
	}
	up->fOut = exp(-s_dist/up->fChi2);
	sum_out += up->fOut;
	if (s_dist<min_s_dist) {
	    min_s_dist = s_dist; 
	    fUwin = up;
	}
    }
    for (J=0;J<fParm.fOutNodes;++J) fOut[J] = 0.0;
    if (sum_out>0) {
	for(up=fU;up<fUbound;++up) {
	    up->fOut /= sum_out;
	    for (J=0;J<fParm.fOutNodes;++J) 
		fOut[J] += up->fOut * up->fWeight[J];
	}
    }

    if (fPlotter) {
	Bool_t good = kTRUE;
	if (out!=0) good = out[0]>fParm.fThreshold;
	fPlotter->AddTestSample(fOut[0],good);
    }

    return fOut;
}


Double_t TSGCS::Train(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    Int_t I,J;
    Double_t S_Err;

    // Recall
    Double_t s_dist;
    Double_t sum_out = 0.0;
    Double_t min_s_dist = DBL_MAX;
    TNeuralNetCell* up;
    
    for(up=fU;up<fUbound;++up){
	Double_t* v = up->fVector;
	Double_t* d = up->fDiff;
	NNO_INTYPE* i = in;
	s_dist = 0.0;
	for (I=0;I<fParm.fInNodes;++I) { 
	    *d = *i++ - *v++; 
	    s_dist += *d * *d; 
	    ++d;
	}
	up->fOut = exp(-s_dist/up->fChi2);
	sum_out += up->fOut;
	if (s_dist<min_s_dist) {
	    min_s_dist = s_dist; 
	    fUwin = up;
	}
    }
    for (J=0;J<fParm.fOutNodes;++J) fOut[J] = 0.0;
    if (sum_out>0) {
	for(up=fU;up<fUbound;++up) {
	    up->fOut /= sum_out;
	    for (J=0;J<fParm.fOutNodes;++J) 
		fOut[J] += up->fOut * up->fWeight[J];
	}
    }

    // Weights update
    TNeuralNetCell* unei;
    Double_t* vwin = fUwin->fVector;
    Double_t* dwin = fUwin->fDiff;
    for (J=0;J<fParm.fInNodes;++J) 
	*vwin++ += *dwin++ * fXB.fWinStep;

    fUwin->fChi2 = 0.0;
    for (I=0;I<fUwin->fNc;++I) {
	unei = (TNeuralNetCell*)fUwin->fC[I].fPtr;
	Double_t* v = unei->fVector;
	Double_t* d = unei->fDiff;
	Double_t s_dist = 0.0;
	Double_t diff;
	vwin = fUwin->fVector;
	for (J=0;J<fParm.fInNodes;++J) {
	    *v += *d++ * fXB.fNeiStep; 
	    diff = *vwin++ - *v++; 
	    s_dist += diff*diff;
	}
	fUwin->fChi2 += s_dist;
	unei->fChi2  += (s_dist-unei->fChi2) * fXB.fNeiCount;
    }
    
    fUwin->fChi2 /= fUwin->fNc;
    
    S_Err = 0.0;
    for (J=0;J<fParm.fOutNodes;++J) {
	Double_t diff = out[J] - fOut[J];
	S_Err += diff * diff;
	for(up=fU;up<fUbound;++up) 
	    up->fWeight[J] += diff*up->fOut * fXB.fNeuStep;
    }
    fUwin->fCount += S_Err;
    Double_t m_err_count = 1.0 - fXB.fErrCount;
    for(up=fU;up<fUbound;++up) { 
	up->fCount *= m_err_count; 
    }
    
    if (fXB.fInsertStep>0) 
	if (fXB.fInsertCount++ == fXB.fInsertStep) {
	    Insert(); 
	    fXB.fInsertCount=0;
	}
	
    if (fXB.fDeleteStep>0) 
	if (fXB.fDeleteCount++ == fXB.fDeleteStep) {
	    Prune(); 
	    fXB.fDeleteCount=0;
	}
	
    fShouldSave = kTRUE;
    
    if (fPlotter) fPlotter->AddTrainSample(out[0],out[0]>fParm.fThreshold);

    return S_Err;
}

Int_t TSGCS::Insert(void) 
{
    Int_t I,J;
    TNeuralNetCell* up;
    TNeuralNetCell* umax;
    TNeuralNetCell* udst;
    TNeuralNetCell* unew;
    
    if (fXB.fCells == fXB.fMaxCells) return 0; //break if there are no cells availiable
    
    //find cell with highest err_count
    Double_t err_count = -1;
    for(up=fU;up<fUbound;++up) 
	if (up->fCount>err_count) {
	    err_count = up->fCount; 
	    umax = up;
	}
	
	//create new cell
	unew = fUbound++;
	++fXB.fCells;
	unew->fNc = 0;
	unew->fState = 0;
	
	//find cell with highest distance to MaxCount
	Double_t s_dist;
	Double_t max_s_dist=-1;
	TNeuralNetCell* unei;
	for (I=0;I<umax->fNc;++I) {
	    s_dist=0;
	    unei=(TNeuralNetCell*)umax->fC[I].fPtr;
	    for (J=0;J<fParm.fInNodes;++J) {
		Double_t diff = umax->fVector[J] - unei->fVector[J]; 
		s_dist += diff*diff;
	    }
	    if (s_dist>max_s_dist) {
		max_s_dist = s_dist; 
		udst = unei;
	    }
	}
	
	TNeuralNetCell::ConnectNew((TNeuralNetCell*)unew,(TNeuralNetCell*)umax,(TNeuralNetCell*)udst,(TNeuralNetCellParameters*)&fXB);
	//TNeuralNetCell::InitVector((TNeuralNetCell*)unew,&B);
	TNeuralNetCell::InitVector((TNeuralNetCell*)unew,(TNeuralNetCell*)umax,(TNeuralNetCell*)udst,&fParm);
	TNeuralNetCell::InitCount ((TNeuralNetCell*)unew);
	TNeuralNetCell::InitWgt   ((TNeuralNetCell*)unew,&fParm);
	TNeuralNetCell::GetSDev   ((TNeuralNetCell*)unew,&fParm);
	
	for(up=fU;up<fUbound;++up) TNeuralNetCell::CheckConnections((TNeuralNetCell*)up);
	return 1;  //insertion was successful
}


void TSGCS::Remove(TNeuralNetCell* urem) 
{
    Int_t I;
    TNeuralNetCell* ulast = --fUbound;
    --fXB.fCells;
    while(urem->fNc)  TNeuralNetCell::Disconnect((TNeuralNetCell*)urem->fC[0].fPtr,(TNeuralNetCell*)urem); //disconnect all neighbours of urem from urem
    //copy ulast into urem
    memcpy(urem->fVector,ulast->fVector,fParm.fInNodes*sizeof(Double_t));
    memcpy(urem->fWeight,ulast->fWeight,fParm.fOutNodes*sizeof(Double_t));
    memcpy(urem->fC,ulast->fC,fXB.fConnectors*sizeof(connector));
    urem->fNc = ulast->fNc;
    urem->fChi2 = ulast->fChi2;
    urem->fCount = ulast->fCount;
    urem->fState = ulast->fState;
    while(ulast->fNc) TNeuralNetCell::Disconnect((TNeuralNetCell*)ulast->fC[0].fPtr,(TNeuralNetCell*)ulast); //disconnect all neighbours of ulast from ulast
    //Connect all former neigbours of ulast to urem
    for (I=0;I<urem->fNc;++I) {
	TNeuralNetCell* unei = (TNeuralNetCell*)urem->fC[I].fPtr;
	unei->fC[unei->fNc++].fPtr = urem;
    }
}

Int_t TSGCS::Prune(void) 
{
    Int_t I,J;
    TNeuralNetCell* up;
    TNeuralNetCell* umax;
    TNeuralNetCell* unei;
    
    //find cell with maximal s_dev
    Double_t max_s_dev=-1;
    for(up=fU;up<fUbound;++up) {
	TNeuralNetCell::GetSDev((TNeuralNetCell*)up,&fParm);
	if (up->fChi2>max_s_dev) { 
	    max_s_dev = up->fChi2; 
	    umax = up; 
	}
    }

    Int_t RemoveCount=1;
    umax->fState |= REMOVE;
    for (I=0;I<umax->fNc;++I) {
	unei = (TNeuralNetCell*)umax->fC[I].fPtr;
	Int_t xnei = 0;
	for (J=0;J<unei->fNc;++J) {
	    Int_t I;
	    TNeuralNetCell* uneinei = (TNeuralNetCell*)unei->fC[J].fPtr;
	    if (uneinei!=umax) {
		for (I=0;I<umax->fNc;++I) 
		    if (uneinei==(TNeuralNetCell*)umax->fC[I].fPtr) break;
		    if (I==umax->fNc) ++xnei;
	    }
	}

	//if unei has no neigbour outside, set REMOVE-flag
	if (xnei==0) {unei->fState|=REMOVE; ++RemoveCount; }
    }

    if (fXB.fCells-RemoveCount<fXB.fMinCells) {
	for(up=fU;up<fUbound;++up) up->fState&=~REMOVE; //clear REMOVE-flag on all cells if network would be shrinking to zero
    } else {
	up = fU;
	while(up<fUbound) { 
	    if (up->fState&REMOVE) 
		Remove(up); 
	    else 
		++up; 
	}
    }

    return 1;
}

