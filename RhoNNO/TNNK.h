#ifndef TNNK_H
#define TNNK_H
// TNNK
//
// Interface to J.P Ernenwein's neural network kernel TNNKernel
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include "VSupervisedNet.h"

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
    TNNK(double learn,double fse,double mu,int innodes,Text_t *hidnodes,int outnodes,std::string netFile);
    TNNK(std::string netFile);
    virtual ~TNNK();                          //destructor of network  (File will be saved)
    double Train(NNO_INTYPE* in,NNO_OUTTYPE* out);
    double* Recall(NNO_INTYPE* in,NNO_OUTTYPE* out);
    
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
    int  fNHiddL;  // number of hidden layers
    float **fValues;     //! array of activations 
    double **fErrors;    //! array of errors
    double **fBiases;    //! array of biases
    int *fNUnits;    //! array of units numbers
    double ***fW;    //! array of weights
    
    int  fNTrainEvents;  // number of events for training 
    int  fNValidEvents;  // number of events for validation
    TNNTree *fValidTree;   // !validation tree 
    double fLearnParam;   // learning parameter
    float fLowerInitWeight;  // minimum weight for initialisation
    float fUpperInitWeight;  // maximum weight for initialisation
    float **fArrayOut;  //! Internal Array with teaching  values for output units
    float *fTeach; //! array of teaching values for outputs
    float **fArrayIn; //! Internal Array with values of input units
    int *fEventsList; //! array of events numbers, for mixing
    int fNTrainCycles; // Number of training cycles done
    double fUseBiases;  // flag for use of biases or not (1=use, 0=no use)
    TRandom fRandom; // Random object used in initialisation and mixing
    int fNWeights; // number of weights in neural network
    double fMu; // backpropagation momentum parameter
    double fFlatSE; // Flat Spot elimination paramater
    double ***fDW;  //! array of delta weights used by backprop momentum
    double **fDB;   //! array of delta biases  used by backprop momentum
    
    void GetArrayEvt(int iEvent)
    {
	int l;      
	for(l=0;l<fNUnits[0];l++)fValues[0][l]=fArrayIn[iEvent][l];
	for(l=0;l<fNUnits[fNHiddL+1];l++)fTeach[l]=fArrayOut[iEvent][l];
    };
    void LearnBackward();   // gradient retropropagation (updates of biases and weights)
    void Forward(); // do a simple forward propagation
    double Error1();// compute the error between forward propagation and teaching
    double ErrorO();// compute the error between forward propagation and teaching
    void FreeVW();  
    void AllocateVW(int nInput, Text_t *hidden, int nOutput);  
    void SetHidden(Text_t *ttext); 
    float Alea();  
    void DeleteArray();  
    
protected: 
    virtual double Sigmoide(double x)
    {  
	if(x> 10.) return 0.99999; // probability MUST be < 1
	if(x<-10.) return 0.;
	return (1./(1.+exp(-x)));
    };
    virtual double SigPrim(double x){return (x*(1.-x));};   
    
public:
    TNNKernel();
    TNNKernel(Text_t *name, int nInput=5, Text_t *hidden=(char *)"6:7:8", int nOutput=4); 
    virtual ~TNNKernel(); // destructor
    virtual void SetKernel(int nInput, Text_t *hidden, int nOutput);
    virtual void SetLearnParam(double learnParam=0.2,double fse=0.,double mu=0.);
    virtual void SetInitParam(float lowerInitWeight=-1., float upperInitWeight=1.);
    virtual void Init();   // init biases and weights
    virtual void PrintS(); // print structure of network
    virtual void Mix();    // mix the events before learning
    virtual double TrainOneCycle();  // one loop on internal events = one cycle
    virtual void ResetCycles(){fNTrainCycles=0;};
    virtual void Export(Text_t *fileName=(char *)"exportNN.dat");
    virtual void Import(Text_t *fileName=(char *)"exportNN.dat");  
    virtual void SetUseBiases(bool trueForUse=1){fUseBiases=(double)trueForUse;};
    virtual void SetRandomSeed(unsigned int seed=0){fRandom.SetSeed(seed);};
    virtual unsigned int GetRandomSeed(){return fRandom.GetSeed();};
    virtual bool IsTrained(){return fNTrainCycles;};
    virtual int GetNTrainCycles(){return fNTrainCycles;};
    virtual int GetNTrainEvents(){return fNTrainEvents;};
    virtual int GetNValidEvents(){return fNValidEvents;};
    virtual void SetArraySize(int s=0);
    virtual void Fill(int iev=0)
    {
	int i;
	for(i=0;i<fNUnits[0];i++)fArrayIn[iev][i]=fValues[0][i];
	for(i=0;i<fNUnits[fNHiddL+1];i++)fArrayOut[iev][i]=fTeach[i];
    }
    virtual float* GetInputAdr(){return fValues[0];};
    virtual void     SetInput(float v,int i){fValues[0][i]=v;};
    virtual int    GetNInput(){return fNUnits[0];};
    virtual int    GetNOutput(){return fNUnits[fNHiddL+1];};
    virtual float  GetOutput(int unit=0){return fValues[fNHiddL+1][unit];};
    virtual float* GetOutputAdr(){return fValues[fNHiddL+1];};
    virtual float* GetTeachAdr(){return fTeach;};
    virtual void     SetTeach(float v,int i){fTeach[i]=v;};
    virtual double GoThrough(){Forward();return ErrorO();};
    virtual float  GetSumO()
    {
	int i; float s=0.;
	for(i=0;i<fNUnits[fNHiddL+1];i++)s+=fValues[fNHiddL+1][i];
	return s;
    };
    //virtual void SetTrainTree(TNNTree *t);  
    //virtual void SetValidTree(TNNTree *t);  
    //virtual double Valid();
    //virtual void TrainNCycles(TNNControlE *conte, int period=5, int nCycles=10);
    virtual int GetNWeights()
    {
	if(!fNUnits)return 0;
	int n=0;
	for(int i=0;i<fNHiddL+1;i++)
	{
	    n+=fNUnits[i]*fNUnits[i+1];
	}
	return n;
    };
    
    ClassDef(TNNKernel,1)

    friend class TNNK;
	
};
#endif
