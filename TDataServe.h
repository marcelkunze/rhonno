#ifndef TDataServe_H
#define TDataServe_H
// TDataServe
//
// A simple datastack for use with the Neural Network Objects (NNO)
// It supports: - separate storing of train and testvectors
//        - mixing of the trainvectors
//        - deleting of bad vectors
//        - reading data from TNtuple
//
// Usage: - Call the normal constructors
//      - Fill the TDataServe objects via
//        - Putvec function
//        - TNtuple_DataRead function
//        - TNtuple_XDataRead function
//      - Call Init function
//
// Author List:
// H.Schmuecker, Bochum University
// M.Kunze, Bochum University, Feb. 01

#include <assert.h>
#include "TNamed.h"

#include <string>

class TNtuple;

class TDataServe : public TNamed {
    
private:  		
    unsigned int	fInvecLen;	//length of inputvector
    unsigned int	fOutvecLen;	//length of outputvectors
    unsigned int	fMaxvecs;
    float**	fInvecAr;	//array of inputputvectors
    float*	fInvecMean;	//mean values
    float*	fInvecScale;//input scale
    float**	fOutvecAr;	//array of outputvectors
    float*	fOutvecMean;//mean values
    float*	fOutvecScale;//output scale
    
    unsigned int	fNumTrnvecs;//Number of trainvectors
    unsigned int	fNumTstvecs;//Number of testvevtors
    unsigned int	fNumvecs;	//Number of all vectors
    unsigned int* fIndexTrn;	//index array
    unsigned int*	fIndexTst;	//index array
    bool	fData_OK;
    bool	fBalance;	//!Balance samples
    
    void	Mix(unsigned int*, const unsigned int);
    void	BlastAr2DD(float**& ar,const unsigned int oldlen1,const unsigned int newlen1,const unsigned int len2);
    void	DelAr2DD(float**& ar,const unsigned int len);
    
public:
    TDataServe();
    TDataServe(std::string name,std::string title,const unsigned int,const unsigned int);
    virtual ~TDataServe();
    void	MixTrn();
    void	Init(const unsigned int);
    void	DataRead(std::string name, const float*, const unsigned int=0, const unsigned int=0);
    unsigned int	GetNumAllvecs() const { return fNumTrnvecs+fNumTstvecs; }
    unsigned int	GetNumTrnvecs() const { return fNumTrnvecs; }
    unsigned int	GetNumTstvecs() const { return fNumTstvecs; }
    unsigned int	GetNumvecs() const { return fNumvecs; }
    void	FillTNtuple(TNtuple&) const;
    void	Reset();
    
    bool	TTreeDataRead(std::string file,std::string tree,std::string in,std::string out,std::string cut="");
    
    void	TNtupleDataRead(TNtuple&,const unsigned int=0,unsigned int=0);
    
    void	TNtupleXDataRead(TNtuple&,const unsigned int,const unsigned int,const unsigned int*,
                             const unsigned int*,const unsigned int=0,unsigned int=0);
    void	SetInvecElem(const unsigned int, const unsigned int, const float);
    void	SetOutvecElem(const unsigned int, const unsigned int, const float);
    void	Putvec(const float*, const float*);
    void	Deletevec(const unsigned int);
    unsigned int	GetInvecLen() const {return fInvecLen;}
    unsigned int	GetOutvecLen() const {return fOutvecLen;}
    float*	GetOutputMean();
    float*	GetOutputScale();
    float*	GetInputMean();
    float*	GetInputScale();
    
    
    float* GetInvec(const unsigned int i) const
    {
        assert(i<fNumvecs);
        return fInvecAr[i];
    }
    
    float* GetOutvec(const unsigned int i) const
    {
        assert(i<fNumvecs);
        return fOutvecAr[i];
    }
    
    float* GetInvecTrn(const unsigned int i) const
    {
        assert(fData_OK && i<fNumTrnvecs);
        return fInvecAr[fIndexTrn[i]];
    }
    
    float* GetInvecTst(const unsigned int i) const
    {
        assert(fData_OK && i<fNumTstvecs);
        return fInvecAr[fIndexTst[i]];
    }
    
    float* GetOutvecTrn(const unsigned int i) const
    {
        assert(fData_OK && i<fNumTrnvecs);
        return fOutvecAr[fIndexTrn[i]];
    }
    
    float* GetOutvecTst(const unsigned int i) const {
        assert(fData_OK && i<fNumTstvecs);
        return fOutvecAr[fIndexTst[i]];
    }
    
    ClassDef(TDataServe,2)	//Database for network training
};
#endif
