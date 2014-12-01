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

#include <cmath>
#include <cstdlib>
#include "TFile.h"
#include "TNtuple.h"
#include "TTreeFormula.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"

#include "RhoNNO/TDataServe.h"

#include <iostream>
#include <fstream>
#include <strstream>
using namespace std;

#ifdef WIN32
//DllImport TRandom *gRandom;
#endif

ClassImp(TDataServe)

// This is using the old streams

void TDataServe::Streamer(TBuffer &b)
{
    // Streamer
    UInt_t i;
    
    if (b.IsReading()){
	Version_t v=b.ReadVersion();
	TNamed::Streamer(b);
	b >> fInvecLen;
	b >> fOutvecLen;
	b >> fNumTrnvecs;
	b >> fNumTstvecs;
	b >> fNumvecs;
	b >> fData_OK;
	fMaxvecs=fNumvecs;
	if (fNumTstvecs) { fIndexTst=new UInt_t [fNumTstvecs]; assert(fIndexTst!=0); }
	if (fNumTrnvecs) { fIndexTrn=new UInt_t [fNumTrnvecs]; assert(fIndexTrn!=0); }
	if (fNumTrnvecs) b.ReadFastArray(fIndexTrn, fNumTrnvecs);
	if (fNumTstvecs) b.ReadFastArray(fIndexTst, fNumTstvecs);    
	BlastAr2DD(fInvecAr, 0, fNumvecs, fInvecLen);
	BlastAr2DD(fOutvecAr, 0, fNumvecs, fOutvecLen);
	for (i=0; i<fNumvecs; i++){
	    b.ReadFastArray(fInvecAr[i], fInvecLen);
	    b.ReadFastArray(fOutvecAr[i], fOutvecLen);
	}
    }
    else{
	b.WriteVersion(TDataServe::IsA());
	TNamed::Streamer(b);
	b << fInvecLen;
	b << fOutvecLen;
	b << fNumTrnvecs;
	b << fNumTstvecs;
	b << fNumvecs;
	b << fData_OK;
	b.WriteFastArray(fIndexTrn, fNumTrnvecs);
	b.WriteFastArray(fIndexTst, fNumTstvecs);
	for (i=0; i<fNumvecs; i++){
	    b.WriteFastArray(fInvecAr[i], fInvecLen);
	    b.WriteFastArray(fOutvecAr[i], fOutvecLen);
	}
    }
}

TBuffer &operator>>(TBuffer &buf, TDataServe *&obj)
{
   obj = (TDataServe *) buf.ReadObject(TDataServe::Class());
   return buf;
}

TDataServe::TDataServe() :
TNamed(),
fInvecLen(0),
fOutvecLen(0),
fInvecScale(0),
fOutvecScale(0),
fMaxvecs(0),
fNumvecs(0),
fNumTrnvecs(0),
fNumTstvecs(0),
fData_OK(kFALSE),
fInvecAr(0),
fInvecMean(0),
fOutvecAr(0),
fOutvecMean(0),
fIndexTrn(0),
fIndexTst(0)
{}

TDataServe::TDataServe(Text_t* name,Text_t* title,const UInt_t in,const  UInt_t out):
fInvecLen(in),
fOutvecLen(out),
fInvecScale(0),
fOutvecScale(0),
TNamed(name,title),
fMaxvecs(0),
fNumvecs(0),
fNumTrnvecs(0),
fNumTstvecs(0),
fData_OK(kFALSE),
fInvecAr(0),
fInvecMean(0),
fOutvecAr(0),
fOutvecMean(0),
fIndexTrn(0),
fIndexTst(0)
{
    // normal constructor
    // Parameters: name, title, inputvectorlength, outputvectorlength
    
}

void TDataServe::FillTNtuple(TNtuple& tup) const
{
    // fill tup with input and outputvectors
    // e.g. TNtuple tup("tup","in out","i1:i2:i3:i4... :o1:o2:o3");  
    Float_t* helpar=new Float_t [fInvecLen+fOutvecLen];
    assert(helpar!=0);
    UInt_t i,j;
    for (i=0; i<fNumvecs; i++){
	for (j=0; j<fInvecLen; j++) helpar[j]=fInvecAr[i][j];
	for (j=0; j<fOutvecLen; j++) helpar[j+fInvecLen]=fOutvecAr[i][j]; 
	tup.Fill(helpar);
    }
    delete [] helpar;
}

TDataServe::~TDataServe()
{
    // destructor
    if (fIndexTrn!=0) { delete [] fIndexTrn; fIndexTrn=0; }
    if (fIndexTst!=0) { delete [] fIndexTst; fIndexTst=0; }
    if (fInvecMean!=0) { delete [] fInvecMean; fInvecMean=0; }
    if (fOutvecMean!=0) { delete [] fOutvecMean; fOutvecMean=0; }
    if (fInvecScale!=0) { delete [] fInvecScale; fInvecScale=0; }
    if (fOutvecScale!=0) { delete [] fOutvecScale; fOutvecScale=0; }
    DelAr2DD(fInvecAr, fNumvecs);
    DelAr2DD(fOutvecAr, fNumvecs);
}

void TDataServe::Reset()
{
    // clears TDataServe object and sets allocated memory free 
    if (fIndexTrn!=0) { delete [] fIndexTrn; fIndexTrn=0; }
    if (fIndexTst!=0) { delete [] fIndexTst; fIndexTst=0; }
    if (fInvecMean!=0) { delete [] fInvecMean; fInvecMean=0; }
    if (fOutvecMean!=0) { delete [] fOutvecMean; fOutvecMean=0; }
    DelAr2DD(fInvecAr, fNumvecs);
    DelAr2DD(fOutvecAr, fNumvecs);
    fNumvecs=0;
    fData_OK=kFALSE;
    fNumTrnvecs=0;
    fNumTstvecs=0;
    fMaxvecs=0;
}

void TDataServe::DelAr2DD(Float_t**& ar,const UInt_t len)
{
    if (ar==0) return;
    for (UInt_t i=0; i<len; i++) delete [] ar[i];
    delete [] ar;
}

void TDataServe::BlastAr2DD(Float_t**& ar,
			    const UInt_t oldlen1,
			    const UInt_t newlen1,
			    const UInt_t len2)
{
    UInt_t i;
    Float_t** temp = ar;
    assert(newlen1>=oldlen1);
    ar = new Float_t* [newlen1];
    assert(ar!=0);
    for (i=0; i<oldlen1; i++) ar[i] = temp[i];
    delete [] temp;
    for (i=oldlen1; i<newlen1; i++){
	ar[i] = new Float_t [len2];
	assert(ar[i]!=0);
    }
}

void TDataServe::Mix(UInt_t* vec, const UInt_t len)
{
    UInt_t i;
    UInt_t temp;
    UInt_t x1=0,x2=0;
    
    for (i=0; i<len; i++){
	x1 = (UInt_t) (gRandom->Rndm(x2)*len);
	x2 = (UInt_t) (gRandom->Rndm(x1)*len);
	temp = vec[x1];
	vec[x1] = vec[x2];
	vec[x2] = temp;
    }
}

void TDataServe::MixTrn()
{
    //mixes trainvectors of TDataServe object.
    assert(fData_OK);
    Mix(fIndexTrn, fNumTrnvecs);
}

void TDataServe::Init(const UInt_t tst)
{
    // Init(number of testvectors)
    // Call Init after every change of data.
    UInt_t i;
    assert(tst<=fNumvecs);
    UInt_t* helpvec = new UInt_t [fNumvecs];
    assert(helpvec!=0);
    if (fIndexTrn!=0) { delete [] fIndexTrn; fIndexTrn=0; }
    if (fIndexTst!=0) { delete [] fIndexTst; fIndexTst=0; }
    if (fInvecMean!=0) { delete [] fInvecMean; fInvecMean=0; }
    if (fOutvecMean!=0) { delete [] fOutvecMean; fOutvecMean=0; }
    if (fInvecMean) { fInvecMean=new Float_t [fInvecLen]; assert(fInvecMean!=0); }
    if (fOutvecMean) { fOutvecMean=new Float_t [fOutvecLen]; assert(fOutvecMean!=0); }
    fNumTstvecs=tst;
    fNumTrnvecs=fNumvecs-tst;
    if (fNumTstvecs) { fIndexTst=new UInt_t [fNumTstvecs]; assert(fIndexTst!=0); }
    if (fNumTrnvecs) { fIndexTrn=new UInt_t [fNumTrnvecs]; assert(fIndexTrn!=0); }
    gRandom->SetSeed(); // Randomize the numbers
    for (i=0; i<fNumvecs; i++) helpvec[i]=i;
    Mix(helpvec, fNumvecs);
    for (i=0; i<fNumTstvecs; i++) fIndexTst[i]=helpvec[fNumvecs-i-1];
    for (i=0; i<fNumTrnvecs; i++) fIndexTrn[i]=helpvec[i];
    fData_OK=kTRUE;
    delete helpvec;
}

void TDataServe::DataRead(const char* name,
			     const Float_t* out,
			     const UInt_t start,
			     const UInt_t laenge)
{
    // forget it
    UInt_t len,i,j;
    ifstream f(name, ios::in);
    if (!f){
	cerr << "CAN'T OPEN FILE: " << name;
	return;
    }
    cerr << "reading file: " << name << endl;
    f.seekg(0L, ios::end);
    len=f.tellg();				 
    f.seekg((long)(start*4*fInvecLen), ios::beg);
    len-=f.tellg();
    len = (UInt_t) 0.25 * len;
    if (laenge>0){
	assert(laenge*fInvecLen<=len); 
	len=laenge*fInvecLen;
    }
    Float_t* fvec=new Float_t [len];
    assert(fvec!=0);
    f.read((char*)(fvec), len*4);
    if (f.fail()){
	delete [] fvec;
	cerr << "READ ERROR IN FILE: " << name;
	exit(5);
    }
    len/=fInvecLen;
    BlastAr2DD(fInvecAr, fNumvecs, fNumvecs+len, fInvecLen);
    BlastAr2DD(fOutvecAr, fNumvecs, fNumvecs+len, fOutvecLen);
    for (i=0; i<len; i++){
	for (j=0; j<fInvecLen; j++) fInvecAr[i+fNumvecs][j]=fvec[i*fInvecLen+j];
	for (j=0; j<fOutvecLen; j++) fOutvecAr[i+fNumvecs][j]=out[j];
    }
    fMaxvecs=fNumvecs+=len;
    delete [] fvec;
    fData_OK=kFALSE;
}

void TDataServe::TNtupleXDataRead(TNtuple& tup,
				   const UInt_t inlen,
				   const UInt_t outlen,			
				   const UInt_t* inselect,
				   const UInt_t* outselect,
				   const UInt_t start,
				   UInt_t len)
{
    // fill TDataServe from TNtuple
    //
    // inlen: length of inselect, i.e. number of selected Ntuple columns 
    // inselect: array containing the indexnumbers of selected TNtuple columns
    //	       (first column: 0), if the number of selected columns is smaller than 
    //	       the length of inputvectors, the skipped elements of inputvectors 
    //	       are filled with zero.	
    // start=0: startevent of TNtuple
    // len=0: Number of events   if len=0, tup is read from 
    //	    start position to last row of Ntuple.    
    UInt_t i,j;
    Float_t* val;
    if (len==0) len=(UInt_t)tup.GetEntries()-start;
    assert(start+len<=tup.GetEntries() && inlen<=fInvecLen && outlen<=fOutvecLen);
    BlastAr2DD(fInvecAr, fNumvecs, fNumvecs+len, fInvecLen);
    BlastAr2DD(fOutvecAr, fNumvecs, fNumvecs+len, fOutvecLen);
    j=tup.GetNvar();
    for (i=0; i<inlen; i++) assert(inselect[i]<j && inselect[i]>0);
    for (i=0; i<outlen; i++) assert(outselect[i]<j && outselect[i]>0);
    for (i=0; i<len; i++){
	tup.GetEvent(start+i,1);
	val=tup.GetArgs();
	for (j=0; j<inlen; j++) fInvecAr[i+fNumvecs][j]=val[inselect[j]];
	for (j=inlen; j<fInvecLen; j++) fInvecAr[i+fNumvecs][j]=0.0;
	for (j=0; j<outlen; j++) fOutvecAr[i+fNumvecs][j]=val[outselect[j]];
	for (j=outlen; j<fOutvecLen; j++) fOutvecAr[i+fNumvecs][j]=0.0;
    }
    fMaxvecs=fNumvecs+=len;
    fData_OK=kFALSE;
}

void TDataServe::TNtupleDataRead(TNtuple& tup,
				  const UInt_t start,
				  UInt_t len)
{
    // fill TDataServe from TNtuple
    // tup("tup","in out","i1:i2:i3:i4... :o1:o2:o3");
    // start=0: startevent of TNtuple
    // len=0: Number of events   if len=0, TNuple is read, from
    //	    startposition to the last row of Ntuple.   
    UInt_t i;
    UInt_t* iv=new UInt_t [fInvecLen];
    UInt_t* ov=new UInt_t [fOutvecLen];
    for (i=0; i<fInvecLen; i++) iv[i]=i;
    for (i=0; i<fOutvecLen; i++) ov[i]=fInvecLen+i;
    TNtupleXDataRead(tup, fInvecLen, fOutvecLen, iv, ov, start, len);
    delete [] iv;
    delete [] ov;
}

void TDataServe::SetInvecElem(const UInt_t ind1,
			      const UInt_t ind2,
			      const Float_t val)
{
    // set single inputvector element 
    // ind1: number of vector
    // ind2: number of element
    assert(ind1<fNumvecs && ind2<fInvecLen);
    fInvecAr[ind1][ind2]=val;
}

void TDataServe::SetOutvecElem(const UInt_t ind1,
			       const UInt_t ind2,
			       const Float_t val)
{
    // set single inputvector element 
    // ind1: number of vector
    // ind2; number of element
    assert(ind1<fNumvecs && ind2<fOutvecLen);
    fOutvecAr[ind1][ind2]=val;
}

void TDataServe::Putvec(const Float_t* invec, const Float_t* outvec)
{
    //Add new input and outputvector
    UInt_t i;
    if (fNumvecs+1>fMaxvecs)
    {
	fMaxvecs+=100;
	Float_t** temp = fInvecAr;
	fInvecAr = new Float_t* [fMaxvecs];
	assert(fInvecAr!=0);
	for (i=0; i<fNumvecs; i++) fInvecAr[i] = temp[i];
	delete [] temp;
	temp=fOutvecAr;
	fOutvecAr = new Float_t* [fMaxvecs];
	assert(fOutvecAr!=0);
	for (i=0; i<fNumvecs; i++) fOutvecAr[i] = temp[i];
	delete [] temp;
    }
    fInvecAr[fNumvecs]=new Float_t [fInvecLen];
    assert(fInvecAr[fNumvecs]!=0);
    fOutvecAr[fNumvecs]=new Float_t [fOutvecLen];
    assert(fOutvecAr[fNumvecs]!=0);
    fInvecMean=new Float_t [fInvecLen];
    assert(fInvecMean!=0);
    fOutvecMean=new Float_t [fOutvecLen];
    assert(fOutvecMean!=0);
    for (i=0; i<fInvecLen; i++) fInvecAr[fNumvecs][i]=invec[i];
    for (i=0; i<fOutvecLen; i++) fOutvecAr[fNumvecs][i]=outvec[i];
    fNumvecs++;
    fData_OK=kFALSE;    
}

void TDataServe::Deletevec(const UInt_t ind)
{
    // Deletes input and outputvector 
    assert(ind<fNumvecs);
    delete [] fInvecAr[ind];
    delete [] fOutvecAr[ind];
    fInvecAr[ind]=fInvecAr[fNumvecs];
    fOutvecAr[ind]=fOutvecAr[fNumvecs];
    fNumvecs--;
    fData_OK=kFALSE;
}


Bool_t TDataServe::TTreeDataRead(const char* file,const char *tree,const char* in, const char* out, const char *c)
{
    TFile f(file);
    TTree *t = (TTree*) f.Get(tree);
    if (t==0) {
	cerr << "TDataServe::TTreeDataRead: Tree "<< tree << " not found in " << file << endl;
	return kFALSE;
    }

    // Set up the cut
    TString cut(c);
    TTreeFormula *cutForm = 0;
    if (cut!="") cutForm = new TTreeFormula("Cut",cut,t);

    // Determine input branches
    TString input(in);
    input.ReplaceAll(":"," ");
    input += " ";
    TTreeFormula *inForm[100];
    istrstream inStream((char *) input.Data());
    TString inName;
    UInt_t nInputs = 0;
    while (inStream >> inName && nInputs<fInvecLen){
	inForm[nInputs++] = new TTreeFormula("Input",inName.Data(),t);
    }

    // Determine output branches
    TString output(out);
    output.ReplaceAll(":"," ");
    output += " ";

    Bool_t tag = (output=="-1" || output=="0" || output=="1");
    Int_t value = 0;
    if (output=="-1") value = -1;
    if (output=="1")  value =  1;

    TTreeFormula *outForm[100];
    UInt_t nOutputs = 0;
    if (!tag) {
	istrstream outStream((char *) output.Data());
	TString outName;
	while (outStream >> outName && nOutputs<fOutvecLen){
	    outForm[nOutputs++] = new TTreeFormula("Output",outName.Data(),t);
	}
    }
    else
	nOutputs = 1;

    Float_t *inValues = new Float_t[fInvecLen];
    Float_t *outValues = new Float_t[fOutvecLen];

    Stat_t nEntries = t->GetEntries();
    for (int n=0;n<nEntries;n++) {

	t->LoadTree(n);

	if (cutForm!=0 && cutForm->EvalInstance(0)!=0.0) continue; // No good sample

	for (UInt_t i=0;i<nInputs;i++) 
	    inValues[i] = (Float_t) inForm[i]->EvalInstance(0);

	for (UInt_t o=0;o<nOutputs;o++)
	    if (!tag)
		outValues[o] = (Float_t) outForm[o]->EvalInstance(0);
	    else
		outValues[o] = (Float_t) value;

	Putvec(inValues,outValues);
    }

    delete [] inValues;
    delete [] outValues;

    //cout << "TDataServe::TTreeDataRead: Accumulated " << fNumvecs << " vectors " << endl;

    return kTRUE;
}

Float_t* TDataServe::GetInputMean()
{
    UInt_t i,j;
    if (fInvecMean==0) fInvecMean = new Float_t[fInvecLen];
    for (i=0; i<fInvecLen; i++) fInvecMean[i] = 0.0;
    if (fNumvecs<=0) return fInvecMean;

    for (i=0; i<fNumvecs; i++){
	for (j=0; j<fInvecLen; j++) fInvecMean[j]+=fInvecAr[i][j];
    }

    for (i=0; i<fInvecLen; i++) fInvecMean[i] /= fNumvecs;

    return fInvecMean;
}

Float_t* TDataServe::GetOutputMean()
{
    UInt_t i,j;
    if (fOutvecMean==0) fOutvecMean = new Float_t[fOutvecLen];
    for (i=0; i<fOutvecLen; i++) fOutvecMean[i] = 0.0;
    if (fNumvecs<=0) return fOutvecMean;

    for (i=0; i<fNumvecs; i++){
	for (j=0; j<fOutvecLen; j++) fOutvecMean[j]+=fOutvecAr[i][j];
    }

    for (i=0; i<fOutvecLen; i++) fOutvecMean[i] /= fNumvecs;

    return fOutvecMean;
}

Float_t* TDataServe::GetInputScale()
{
    UInt_t i,j;
    if (fInvecScale==0) fInvecScale = new Float_t[fInvecLen];
    for (i=0; i<fInvecLen; i++) fInvecScale[i] = 1.0;
    if (fNumvecs<=0) return fInvecScale;

    TCanvas *c1 = new TCanvas("c1","Input parameter space",200,100,700,900);
    TPad *pmain = new TPad("pmain","pmain",0,0,1,.8);
    pmain->SetFillColor(20);
    pmain->Draw();
    pmain->Divide(3,3);
    for (j=0; j<fInvecLen; j++) {
	Float_t max=FLT_MIN,min=FLT_MAX;
	for (i=0; i<fNumvecs; i++){
	    Float_t vec = fInvecAr[i][j];
	    if (vec>max) max = vec;
	    else if (vec<min) min = vec;
	}
	//if (min==0) min = -1.;
	//if (max==0) max = 1.;
	pmain->cd(j%9+1);
	TF1 fit("g","gaus");
	fit.SetLineColor(kRed);
	TString name("Input ");
	name += j;
	TH1F *hist = new TH1F(name,name,400,2*min,2*max);
	hist->SetFillColor(45);
	for (i=0; i<fNumvecs; i++) hist->Fill(fInvecAr[i][j]);
	hist->Fit("g");
	hist->Draw();
	c1->Update();
	Double_t par[3];
	fit.GetParameters(par);
	fInvecScale[j] = 1. / fabs(par[1]); // Position of peak
	
    }

    return fInvecScale;
}

Float_t* TDataServe::GetOutputScale()
{
    UInt_t i,j;
    if (fOutvecScale==0) fOutvecScale = new Float_t[fOutvecLen];
    for (i=0; i<fOutvecLen; i++) fOutvecScale[i] = 1.0;
    if (fNumvecs<=0) return fOutvecScale;

    for (j=0; j<fOutvecLen; j++) {
	Float_t max=FLT_MIN,min=FLT_MAX;
	for (i=0; i<fNumvecs; i++){
	    Float_t vec = fInvecAr[i][j];
	    if (vec>100 || vec<-100) continue;
	    if (vec>max) max = vec;
	    else if (vec<min) min = vec;
	}
	fOutvecScale[j] = 2.0 / (max - min);
    }

    return fOutvecScale;
}

