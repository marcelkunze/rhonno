//////////////////////////////////////////////////////////////////////////
//									//
// VUnsupervisedNet							//
//									//
// Base classes for unsupervised learning				//
// Abstract base class of all unsupervised networks			//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VUnsupervisedNet.h"

ClassImp(VUnsupervisedNet)

Long_t  VUnsupervisedNet::TrainEpoch(FILE* file) 
{
    if (file==0) return -1;
    NNO_INTYPE* Buf = new NNO_INTYPE[fParm.fInNodes];
    TestPointer(Buf);
    Long_t records = 0;
    rewind(file);
    while(!feof(file)) {
	if (fread(Buf,sizeof(NNO_INTYPE),fParm.fInNodes,file)==fParm.fInNodes) {
	    Train(Buf);
	    ++records;
	}
    }
    delete[] Buf;
    return records;
}
