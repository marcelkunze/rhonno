//////////////////////////////////////////////////////////////////////////
//									//
// Routines for Connected-Cell-Networks					//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/TNeuralNetCell.h"
#include "RhoNNO/VNeuralNet.h"

ClassImp(TNeuralNetCellParameters)

TNeuralNetCellParameters::TNeuralNetCellParameters() 
{
    fWinStep = 0;
    fNeiStep = 0;
    fNeuStep = 0;
    fCells = 0;
    fConnectors = 0;
    fWinCount = 0;
    fMinCells = 0;
    fMaxCells = 0;
    fInsertStep = 0;
    fDeleteStep = 0;
    fInsertCount = 0;
    fDeleteCount = 0;
    fEdgeCount = 0;
    fErrCount = 0;
    fNeiCount = 0;
    fMinCount = 0;
    fMainWinCount = 0;
    fMainErrCount = 0;
    fMainEdgeCount = 0;
}

ClassImp(TNeuralNetCell)

TNeuralNetCell::TNeuralNetCell()
{
    fVector = 0;
    fID = 0;
    fCount = 0;
    fDiff  = 0;
    fState  = 0;
    fC  = 0;
    fNc = 0;
    fChi2 = 0;
    fWeight = 0;
    fOut = 0;
    fAge = 0;
    fClass = 0;
}

void TNeuralNetCell::Disconnect(TNeuralNetCell* up1,TNeuralNetCell* up2) 
{
     int I;
     for (I=0;I<up1->fNc;++I) if (up1->fC[I].fPtr==up2) break;
     if (I<up1->fNc) up1->fC[I]=up1->fC[--up1->fNc];
     for (I=0;I<up2->fNc;++I) if (up2->fC[I].fPtr==up1) break;
     if (I<up2->fNc) up2->fC[I]=up2->fC[--up2->fNc];
}

void TNeuralNetCell::Connect(TNeuralNetCell* up1,TNeuralNetCell* up2,TNeuralNetCellParameters* XB) 
{
     if ((up1->fNc==XB->fConnectors)||(up2->fNc==XB->fConnectors)) return;
     up1->fC[up1->fNc++].fPtr=up2;
     up2->fC[up2->fNc++].fPtr=up1;
}

//connect new cell with common neighbours of unit1 and unit2
void TNeuralNetCell::ConnectNew(TNeuralNetCell* unew,TNeuralNetCell* unit1,TNeuralNetCell* unit2,TNeuralNetCellParameters* XB) 
{
    int I,J;
    for (I=0;I<unit1->fNc;++I) 
	for (J=0;J<unit2->fNc;++J)
	if (unit1->fC[I].fPtr==unit2->fC[J].fPtr)
	   Connect((TNeuralNetCell*)unit1->fC[I].fPtr,unew,XB);
    Connect(unew,unit1,XB);
    Connect(unew,unit2,XB);
    Disconnect(unit1,unit2);
}

//this procedure places unew beween all neighbours in the centre of gravity
void TNeuralNetCell::InitVector(TNeuralNetCell* unew,TNeuralNetParameters* B) 
{
     int I,J;
     for (I=0;I<B->fInNodes;++I) unew->fVector[I]=0;
     for (J=0;J<unew->fNc;++J) {
	 TNeuralNetCell* unei = (TNeuralNetCell*) unew->fC[J].fPtr;
	 for (I=0;I<B->fInNodes;++I) 
	     unew->fVector[I] += unei->fVector[I];
     }
     for (I=0;I<B->fInNodes;++I) unew->fVector[I] /= unew->fNc;
}

//this procedure places unew beween U1 and U2
void TNeuralNetCell::InitVector(TNeuralNetCell* unew,TNeuralNetCell* U1,TNeuralNetCell* U2,TNeuralNetParameters* B) 
{
     int I;
     for (I=0;I<B->fInNodes;++I) 
	 unew->fVector[I] = (U1->fVector[I] + U2->fVector[I]) * 0.5;
}

void TNeuralNetCell::InitWgt(TNeuralNetCell* unew,TNeuralNetParameters* B) 
{
     int I,J;
     Double_t fract=1.0 / (1.0 + unew->fNc);
     for (I=0;I<B->fOutNodes;++I) unew->fWeight[I]=0;
     for (J=0;J<unew->fNc;++J) {
	TNeuralNetCell* unei = (TNeuralNetCell*) unew->fC[J].fPtr;
	for (I=0;I<B->fOutNodes;++I) {
	    Double_t WgtFract = unei->fWeight[I] * fract;
	    unew->fWeight[I] += WgtFract;
	    unei->fWeight[I] -= WgtFract;
	}
     }
}

void TNeuralNetCell::GetSDev(TNeuralNetCell* unit,TNeuralNetParameters* B) 
{
     int I,J;
     unit->fChi2=0;
     for (I=0;I<unit->fNc;++I) {
	TNeuralNetCell* unei=(TNeuralNetCell*)unit->fC[I].fPtr;
	Double_t sdist = 0.0;
	for (J=0;J<B->fInNodes;++J) { 
	  Double_t d = unei->fVector[J] - unit->fVector[J]; 
	  sdist += d * d; 
	}
	unit->fChi2 += sdist;
     }
     unit->fChi2 /= unit->fNc;
}

void TNeuralNetCell::InitSDev(TNeuralNetCell* unew,TNeuralNetParameters* B) 
{
     int I;
     GetSDev(unew,B);
     for (I=0;I<unew->fNc;++I) GetSDev((TNeuralNetCell*)unew->fC[I].fPtr,B);
}

void TNeuralNetCell::InitCount(TNeuralNetCell* unew) 
{
    int I;
    Double_t fract = 1.0 / (Double_t) (unew->fNc+1);
    unew->fCount=0;
    for (I=0;I<unew->fNc;++I) {
	TNeuralNetCell* unei=(TNeuralNetCell*)unew->fC[I].fPtr;
	Double_t errFract=((TNeuralNetCell*)unew->fC[I].fPtr)->fCount * fract;
	unei->fCount -= errFract;
	unew->fCount += errFract;
    }
}


void TNeuralNetCell::CheckConnections(TNeuralNetCell* unit) 
{
     int I,J;
     for (I=0;I<unit->fNc;++I) {
         TNeuralNetCell* unei=(TNeuralNetCell*)unit->fC[I].fPtr;
	 for (J=0;J<unei->fNc;++J) {if ((TNeuralNetCell*)unei->fC[J].fPtr==unit) break; }
	 if (J==unei->fNc) fprintf(stdout,"one-way-connection from cell %i to cell %i",unit->fID,unei->fID);
     }
}


//               -------------  FILE I/O  ----------------

const char *f_id          = "\nunit number    %i\n";
const char *f_fChi2       = "square deviation %le\n";
const char *f_count       = "count            %le\n";
const char *f_vector      = "vector\n";
const char *f_connections = "\nconnections %i\n";
const char *f_con_list    = "connected with ";


void TNeuralNetCell::WriteUnitText(FILE* file,TNeuralNetCell* unit,TNeuralNetParameters* B) 
{
     int I;
     fprintf(file,f_id,unit->fID);
     fprintf(file,f_fChi2,unit->fChi2);
     fprintf(file,f_count,unit->fCount);
     fprintf(file,f_vector);
     for (I=0;I<B->fInNodes;++I) fprintf(file,"%le\n",unit->fVector[I]);
     fprintf(file,f_connections,unit->fNc);
     fprintf(file,f_con_list);
     for (I=0;I<unit->fNc;++I) fprintf(file,"%i ",((TNeuralNetCell*)(unit->fC[I].fPtr))->fID);
}

void TNeuralNetCell::ReadUnitText(FILE* file,TNeuralNetCell* unit,TNeuralNetParameters* B) 
{
     int I;
     fscanf(file,f_id,&unit->fID);
     fscanf(file,f_fChi2,&unit->fChi2);
     fscanf(file,f_count,&unit->fCount);
     fscanf(file,f_vector);
     for (I=0;I<B->fInNodes;++I) fscanf(file,"%le",&unit->fVector[I]);
     fscanf(file,f_connections,&unit->fNc);
     fscanf(file,f_con_list);
     for (I=0;I<unit->fNc;++I) fscanf(file,"%i",&unit->fC[I].fID);
}

void TNeuralNetCell::WriteUnitBinary(FILE* file,TNeuralNetCell* unit,TNeuralNetParameters* B) 
{
     fwrite(&unit->fID,sizeof(unit->fID),1,file);
     fwrite(&unit->fChi2,sizeof(unit->fChi2),1,file);
     fwrite(&unit->fCount,sizeof(unit->fCount),1,file);
     fwrite(unit->fVector,sizeof(Double_t),B->fInNodes,file);
     fwrite(&unit->fNc,sizeof(unit->fNc),1,file);
     int I;
     for (I=0;I<unit->fNc;++I) {
	 Int_t i = ((TNeuralNetCell*)(unit->fC[I].fPtr))->fID;
	 fwrite(&i,sizeof(i),1,file);
     }
}

void TNeuralNetCell::ReadUnitBinary(FILE* file,TNeuralNetCell* unit,TNeuralNetParameters* B) 
{
     fread(&unit->fID,sizeof(unit->fID),1,file);
     fread(&unit->fChi2,sizeof(unit->fChi2),1,file);
     fread(&unit->fCount,sizeof(unit->fCount),1,file);
     fread(unit->fVector,sizeof(Double_t),B->fInNodes,file);
     fread(&unit->fNc,sizeof(unit->fNc),1,file);
     fread(unit->fC,sizeof(connector),unit->fNc,file);
}


