#ifndef TNNK_H
#define TNNK_H
//////////////////////////////////////////////////////////////////////////
//									//
// TNNK									//
//									//
// Interface to J.P Ernenwein's neural network kernel TNNKernel		//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// M.Kunze, Bochum University						//
// (C) Copyright 2001, Ruhr-University Bochum.				//
//									//
//////////////////////////////////////////////////////////////////////////

#include "RhoNNO/VSupervisedNet.h"

class TNNKernel;

class TNNK : public VSupervisedNet {
private:
    TNNKernel *fKernel;

    void AllocNet(void);
    void InitNet(void);
    void WriteText(void);
    void WriteBinary(void);
    void ReadText(void);
    void ReadBinary(void);
public:
    TNNK() {};
    TNNK(Double_t learn,Double_t fse,Double_t mu,Int_t innodes,Text_t *hidnodes,Int_t outnodes,const char* netFile);
    TNNK(const char* netFile);
    virtual ~TNNK();                          //destructor of network  (File will be saved)
    Double_t Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    Double_t* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    
    ClassDef(TNNK,1)	    // Interface to TNNKernel
};

//////////////////////////////////////////////////////////////////
//
//  Neural Network classes :
//  TNNKernel
//                             J.P. Ernenwein (rnenwein@in2p3.fr)
//////////////////////////////////////////////////////////////////

#include <math.h>
#include "TNamed.h"
#include "TRandom.h"
#include "TTree.h"

class TNNTree;
class TNNControlE;

class TNNKernel : public TNamed 
{
    
private:
    Int_t  fNHiddL;  // number of hidden layers
    Float_t **fValues;     //! array of activations 
    Double_t **fErrors;    //! array of errors
    Double_t **fBiases;    //! array of biases
    Int_t *fNUnits;    //! array of units numbers
    Double_t ***fW;    //! array of weights
    
    Int_t  fNTrainEvents;  // number of events for training 
    Int_t  fNValidEvents;  // number of events for validation
    TNNTree *fValidTree;   // !validation tree 
    Double_t fLearnParam;   // learning parameter
    Float_t fLowerInitWeight;  // minimum weight for initialisation
    Float_t fUpperInitWeight;  // maximum weight for initialisation
    Float_t **fArrayOut;  //! Internal Array with teaching  values for output units
    Float_t *fTeach; //! array of teaching values for outputs
    Float_t **fArrayIn; //! Internal Array with values of input units
    Int_t *fEventsList; //! array of events numbers, for mixing
    Int_t fNTrainCycles; // Number of training cycles done
    Double_t fUseBiases;  // flag for use of biases or not (1=use, 0=no use)
    TRandom fRandom; // Random object used in initialisation and mixing
    Int_t fNWeights; // number of weights in neural network
    Double_t fMu; // backpropagation momentum parameter
    Double_t fFlatSE; // Flat Spot elimination paramater
    Double_t ***fDW;  //! array of delta weights used by backprop momentum
    Double_t **fDB;   //! array of delta biases  used by backprop momentum
    
    void GetArrayEvt(Int_t iEvent)
    {
	Int_t l;      
	for(l=0;l<fNUnits[0];l++)fValues[0][l]=fArrayIn[iEvent][l];
	for(l=0;l<fNUnits[fNHiddL+1];l++)fTeach[l]=fArrayOut[iEvent][l];
    };
    void LearnBackward();   // gradient retropropagation (updates of biases and weights)
    void Forward(); // do a simple forward propagation
    Double_t Error();// compute the error between forward propagation and teaching
    Double_t ErrorO();// compute the error between forward propagation and teaching
    void FreeVW();  
    void AllocateVW(Int_t nInput, Text_t *hidden, Int_t nOutput);  
    void SetHidden(Text_t *ttext); 
    Float_t Alea();  
    void DeleteArray();  
    
protected: 
    virtual Double_t Sigmoide(Double_t x)
    {  
	if(x> 10.) return 0.99999; // probability MUST be < 1
	if(x<-10.) return 0.;
	return (1./(1.+exp(-x)));
    };
    virtual Double_t SigPrim(Double_t x){return (x*(1.-x));};   
    
public:
    TNNKernel();
    TNNKernel(Text_t *name, Int_t nInput=5, Text_t *hidden=(char *)"6:7:8", Int_t nOutput=4); 
    virtual ~TNNKernel(); // destructor
    virtual void SetKernel(Int_t nInput, Text_t *hidden, Int_t nOutput);
    virtual void SetLearnParam(Double_t learnParam=0.2,Double_t fse=0.,Double_t mu=0.);
    virtual void SetInitParam(Float_t lowerInitWeight=-1., Float_t upperInitWeight=1.);
    virtual void Init();   // init biases and weights
    virtual void PrintS(); // print structure of network
    virtual void Mix();    // mix the events before learning
    virtual Double_t TrainOneCycle();  // one loop on internal events = one cycle
    virtual void ResetCycles(){fNTrainCycles=0;};
    virtual void Export(Text_t *fileName=(char *)"exportNN.dat");
    virtual void Import(Text_t *fileName=(char *)"exportNN.dat");  
    virtual void SetUseBiases(Bool_t trueForUse=1){fUseBiases=(Double_t)trueForUse;};
    virtual void SetRandomSeed(UInt_t seed=0){fRandom.SetSeed(seed);};
    virtual UInt_t GetRandomSeed(){return fRandom.GetSeed();};    
    virtual Bool_t IsTrained(){return fNTrainCycles;};
    virtual Int_t GetNTrainCycles(){return fNTrainCycles;};
    virtual Int_t GetNTrainEvents(){return fNTrainEvents;};
    virtual Int_t GetNValidEvents(){return fNValidEvents;};
    virtual void SetArraySize(Int_t s=0);
    virtual void Fill(Int_t iev=0)
    {
	Int_t i;
	for(i=0;i<fNUnits[0];i++)fArrayIn[iev][i]=fValues[0][i];
	for(i=0;i<fNUnits[fNHiddL+1];i++)fArrayOut[iev][i]=fTeach[i];
    }
    virtual Float_t* GetInputAdr(){return fValues[0];};
    virtual void     SetInput(Float_t v,Int_t i){fValues[0][i]=v;};
    virtual Int_t    GetNInput(){return fNUnits[0];};
    virtual Int_t    GetNOutput(){return fNUnits[fNHiddL+1];};
    virtual Float_t  GetOutput(Int_t unit=0){return fValues[fNHiddL+1][unit];};
    virtual Float_t* GetOutputAdr(){return fValues[fNHiddL+1];};
    virtual Float_t* GetTeachAdr(){return fTeach;};
    virtual void     SetTeach(Float_t v,Int_t i){fTeach[i]=v;};
    virtual Double_t GoThrough(){Forward();return ErrorO();};
    virtual Float_t  GetSumO()
    {
	Int_t i; Float_t s=0.;
	for(i=0;i<fNUnits[fNHiddL+1];i++)s+=fValues[fNHiddL+1][i];
	return s;
    };
    //virtual void SetTrainTree(TNNTree *t);  
    //virtual void SetValidTree(TNNTree *t);  
    //virtual Double_t Valid();
    //virtual void TrainNCycles(TNNControlE *conte, Int_t period=5, Int_t nCycles=10);
    virtual Int_t GetNWeights()
    {
	if(!fNUnits)return 0;
	Int_t n=0;
	for(Int_t i=0;i<fNHiddL+1;i++)
	{
	    n+=fNUnits[i]*fNUnits[i+1];
	}
	return n;
    };
    
    ClassDef(TNNKernel,1)

    friend class TNNK;
	
};
#endif
