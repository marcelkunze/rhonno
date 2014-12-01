#ifndef TDataServe_H
#define TDataServe_H
//////////////////////////////////////////////////////////////////////////
//									//
// TDataServe								//
//									//
// A simple datastack for use with the Neural Network Objects (NNO)	//
// It supports: - separate storing of train and testvectors		//
//		- mixing of the trainvectors				//
//		- deleting of bad vectors				//
//		- reading data from TNtuple				//
//									//
// Usage: - Call the normal constructors				//
//	  - Fill the TDataServe objects via				//
//		- Putvec function					//
//		- TNtuple_DataRead function				//
//		- TNtuple_XDataRead function				//
//	  - Call Init function						//
//									//
// Author List:								//
// H.Schmücker, Bochum University					//
// M.Kunze, Bochum University, Feb. 01					//
// Copyright (C) 1999-2001, Ruhr-University Bochum.			//
//									//
//									//
//////////////////////////////////////////////////////////////////////////

#include <assert.h>
#include "TNamed.h"

class TNtuple;

class TDataServe : public TNamed {
    
private:  		
    UInt_t	fInvecLen;	//length of inputvector
    UInt_t	fOutvecLen;	//length of outputvectors
    UInt_t	fMaxvecs;	
    Float_t**	fInvecAr;	//array of inputputvectors
    Float_t*	fInvecMean;	//mean values
    Float_t*	fInvecScale;	//input scale
    Float_t**	fOutvecAr;	//array of outputvectors
    Float_t*	fOutvecMean;	//mean values
    Float_t*	fOutvecScale;	//output scale
    
    UInt_t	fNumTrnvecs;	//Number of trainvectors
    UInt_t	fNumTstvecs;	//Number of testvevtors
    UInt_t	fNumvecs;	//Number of all vectors
    UInt_t* 	fIndexTrn;	//index array
    UInt_t*	fIndexTst;	//index array
    Bool_t	fData_OK;
    Bool_t	fBalance;	//!Balance samples
    
    void	Mix(UInt_t*, const UInt_t);
    void	BlastAr2DD(Float_t**& ar,const UInt_t oldlen1,const UInt_t newlen1,const UInt_t len2);
    void	DelAr2DD(Float_t**& ar,const UInt_t len);
    
public:
    TDataServe();
    TDataServe(Text_t* name,Text_t* title,const UInt_t,const UInt_t);			
    virtual ~TDataServe();
    void	MixTrn();
    void	Init(const UInt_t);				
    void	DataRead(const char*, const Float_t*, const UInt_t=0, const UInt_t=0);	
    UInt_t	GetNumAllvecs() const { return fNumTrnvecs+fNumTstvecs; }
    UInt_t	GetNumTrnvecs() const { return fNumTrnvecs; }
    UInt_t	GetNumTstvecs() const { return fNumTstvecs; }
    UInt_t	GetNumvecs() const { return fNumvecs; } 
    void	FillTNtuple(TNtuple&) const;
    void	Reset();

    Bool_t	TTreeDataRead(const char* file,const char *tree,const char* in, const char* out, const char *cut="");

    void	TNtupleDataRead(TNtuple&,const UInt_t=0,UInt_t=0);
    
    void	TNtupleXDataRead(TNtuple&,const UInt_t,const UInt_t,const UInt_t*,
			    const UInt_t*,const UInt_t=0,UInt_t=0);
    void	SetInvecElem(const UInt_t, const UInt_t, const Float_t);
    void	SetOutvecElem(const UInt_t, const UInt_t, const Float_t);
    void	Putvec(const Float_t*, const Float_t*);
    void	Deletevec(const UInt_t);
    UInt_t	GetInvecLen() const {return fInvecLen;}
    UInt_t	GetOutvecLen() const {return fOutvecLen;}
    Float_t*	GetOutputMean();
    Float_t*	GetOutputScale();
    Float_t*	GetInputMean();
    Float_t*	GetInputScale();

    
    Float_t* GetInvec(const UInt_t i) const 
    {
	assert(i<fNumvecs);
	return fInvecAr[i];
    }
    
    Float_t* GetOutvec(const UInt_t i) const 
    {
	assert(i<fNumvecs);
	return fOutvecAr[i];
    }  
    
    Float_t* GetInvecTrn(const UInt_t i) const 
    {
	assert(fData_OK && i<fNumTrnvecs);
	return fInvecAr[fIndexTrn[i]];	
    }

    Float_t* GetInvecTst(const UInt_t i) const 
    {
	assert(fData_OK && i<fNumTstvecs);
	return fInvecAr[fIndexTst[i]];
    }

    Float_t* GetOutvecTrn(const UInt_t i) const 
    {
	assert(fData_OK && i<fNumTrnvecs);
	return fOutvecAr[fIndexTrn[i]];
    }

    Float_t* GetOutvecTst(const UInt_t i) const {
	assert(fData_OK && i<fNumTstvecs);
	return fOutvecAr[fIndexTst[i]];
    }
    
    ClassDef(TDataServe,2)	//Database for network training	  
};
#endif
