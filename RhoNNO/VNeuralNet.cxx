//////////////////////////////////////////////////////////////////////////
//									//
// VNeuralNet								//
//									//
// Base classes for unsupervised and supervised networks		//
// Partof the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// Johannes Steffens, Bochum University					//
// M.Kunze, Bochum University						//
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.	//
//									//
//////////////////////////////////////////////////////////////////////////

static const char* NNO_VERSION="1.3ROOT";

#include "TRandom.h"

#include "RhoNNO/VNeuralNet.h"
#include "RhoNNO/VNeuralNetPlotter.h"
#include "RhoNNO/TDataServe.h"

#ifdef WIN32
//DllImport TRandom *gRandom;
#endif

ClassImp(TNeuralNetParameters)

#include <iostream>
#include <cstdlib>
using namespace std;

TNeuralNetParameters::TNeuralNetParameters() 
{
    for (int i=0;i<9;i++) fNetId[i] = 0;
    fLayers = 0;
    fInScale = 1.0;
    fInNodes = 0;
    fOutNodes = 0;
    fLearnStep = 0.01;
    fTransferId = TNeuralNetParameters::TR_FERMI;
    fPerceptronId = 0;
    fThreshold = 0.0;
    fMu = 0.0;
    fFse = 0.0;
}

ClassImp(VNeuralNet)

VNeuralNet::VNeuralNet() 
  : TNamed("NNO","NNO"), fParm(), fPlotter(0), fOwnPlotter(kFALSE), fBalance(kFALSE)
{
    fShouldSave = kFALSE; 
    fFile  = 0;
    fOut = 0;
}

VNeuralNet::VNeuralNet(const char* netID,Int_t innodes,Int_t outnodes,const char* netFile)
 : TNamed(netID,netID), fParm(), fPlotter(0), fOwnPlotter(kFALSE), fBalance(kFALSE) 
{
#ifndef NNORAND
    gRandom->SetSeed(); // Randomize the numbers
#endif
    fFilename = netFile;
    strncpy(fParm.fNetId,netID,9);
    fShouldSave  = kTRUE;
    fFiletype    = FILE_TEXT;
    fParm.fInNodes  = innodes;
    fParm.fOutNodes = outnodes;
    fFile   = 0;
    fOut    = 0;
    if (outnodes>0) {
	fOut = new Double_t[fParm.fOutNodes];
	TestPointer(fOut);
    }
}

VNeuralNet::VNeuralNet(const char* netFile)
 : TNamed(netFile,netFile), fParm(), fPlotter(0), fOwnPlotter(kFALSE), fBalance(kFALSE) 
{
    fFilename   = netFile;
    fShouldSave = kFALSE;
    fFile	= 0;
    fOut        = 0;
}

VNeuralNet::~VNeuralNet() 
{ 
    if (fOut!=0) { delete[] fOut; fOut = 0; }
    if (fOwnPlotter) { delete fPlotter; fPlotter = 0; }
}	

void VNeuralNet::Save() 
{
    fFile = fopen(fFilename,"wb");
    WriteNet();
}

void VNeuralNet::Save(char* file) 
{
    fFile = fopen(file,"wb");
    WriteNet();
}

void VNeuralNet::WriteNet() 
{
    char ftype[16];
    if (fFiletype==FILE_BINARY) strcpy(ftype,"binary"); else strcpy(ftype,"text");
    if (fFile==0) { cerr << "VNeuralNet::WriteNet:: Could not open for writing " << fFilename << endl; return; }
    fprintf(fFile,"C++  NEURAL NETWORK OBJECTS   VERSION %s\n(C) Copyright Johannes Steffens\nFiletype %s\n",NNO_VERSION,ftype);
    if (fFiletype==FILE_BINARY) WriteNetBinary(); else WriteNetText();
    if (fFiletype==FILE_BINARY) WriteBinary();     else WriteText();
    fclose(fFile);
}

void VNeuralNet::WriteNetText() 
{
    fprintf(fFile,"\nnetwork id  %s\n",fParm.fNetId);
    fprintf(fFile,"innodes     %i\n",fParm.fInNodes);
    fprintf(fFile,"outnodes    %i\n",fParm.fOutNodes);
}

void VNeuralNet::WriteNetBinary() 
{
    fwrite(&fParm,sizeof(TNeuralNetParameters),1,fFile);
}

void VNeuralNet::ReadNet(const char* netID) 
{
    fFile = fopen(fFilename,"rb");
    if (fFile==0) Errorf((char *)"file %s not found",(char *)fFilename.Data());
    char ftype[16];
    char Version[16];
    fscanf(fFile,"C++  NEURAL NETWORK OBJECTS   VERSION %s\n(C) Copyright Johannes Steffens\nFiletype %s\n",Version,ftype);
    if      (!strcmp(ftype,"binary")) fFiletype = FILE_BINARY;
    else if (!strcmp(ftype,"text"))   fFiletype = FILE_TEXT;
    else Errorf((char *)"illegal fileformat: %s",(char *)fFilename.Data());

    if (fFiletype==FILE_BINARY) 
	ReadNetBinary(); 
    else 
	ReadNetText();

    fParm.fNetId[4]=0;
    if (strcmp(netID,fParm.fNetId)) {
	fclose(fFile);
	Errorf((char *)"file %s  (incompatible network)\nnetwork ID is %s and should be %s",fFilename.Data(),fParm.fNetId,netID);
    }

    if (strcmp(Version,NNO_VERSION)) {
	fclose(fFile);
	Errorf((char *)"illegal NNO version number of file %s\nversion number is %s and should be %s",fFilename.Data(),Version,(char *)NNO_VERSION);
    }

    if (fFiletype==FILE_BINARY)
	ReadBinary();
    else
	ReadText();

    fOut = new Double_t[fParm.fOutNodes];

    TestPointer(fOut);

    fclose(fFile);
}

void VNeuralNet::ReadNetText() 
{
     fscanf(fFile,"\nnetwork id  %s\n",fParm.fNetId);
     fscanf(fFile,"innodes     %i\n",&fParm.fInNodes);
     fscanf(fFile,"outnodes    %i\n",&fParm.fOutNodes);
}

void VNeuralNet::ReadNetBinary() 
{
    fread(&fParm,sizeof(TNeuralNetParameters),1,fFile);
}

void VNeuralNet::Errorf(char* format,...) 
{
     va_list ap;
     va_start(ap, format);
     char MainFormat[256];
     sprintf(MainFormat,"NNO ERROR: %s\n",format);
     vfprintf(stderr,MainFormat,ap);
     exit(1);
}

void VNeuralNet::Warningf(FILE* f,char* format,...) 
{
     va_list ap;
     va_start(ap, format);
     char MainFormat[256];
     sprintf(MainFormat,"NNO WARNING: %s\n",format);
     vfprintf(f,MainFormat,ap);
}

void VNeuralNet::Messagef(FILE* f,char* format,...) 
{
     va_list ap;
     va_start(ap, format);
     char MainFormat[256];
     sprintf(MainFormat,"NNO INFO: %s\n",format);
     vfprintf(f,MainFormat,ap);
}

double VNeuralNet::Random(void) 
{
#ifdef NNORAND
//  Machine independent random number generator.
//  Produces uniformly-distributed floating points between 0 and 1.
//  Identical sequence on all machines of >= 32 bits.
//  Universal version (Fred james 1985).
//  Return numbers in the range -0.5..0.5 (MK)

   const float kCONS = 4.6566128730774E-10;
   const int kMASK31 = 2147483647;
   static unsigned int fSeed = 65539;
 
   fSeed *= 69069;
      // keep only lower 31 bits
   fSeed &= kMASK31;
      // Set lower 8 bits to zero to assure exact float
   int jy = (fSeed/256)*256;
   double random = kCONS*jy;
   return 1.0 - 2.0*random;
#else
    return 1.0 - 2.0*gRandom->Rndm();
#endif
}

void VNeuralNet::TestPointer(void* Ptr) 
{
     if (Ptr==0) Errorf((char *)"not enough memory");
}

void VNeuralNet::SetupPlots(VNeuralNetPlotter *plotter)
{
    if (plotter==0) {
	cout << "Instantiating plotter for " << GetName() << endl;
	if (fOwnPlotter) delete fPlotter;
	fPlotter = new TSimpleNeuralNetPlotter(GetName());
	fOwnPlotter = kTRUE;
    }

    fPlotter->Initialize();
}

void VNeuralNet::FillPlots(Double_t trn, Double_t tst)
{
    if (fPlotter==0) return;
    fPlotter->AddTrainSample(trn,kTRUE);
    fPlotter->AddTestSample(tst,kTRUE);
}

void VNeuralNet::ShowPlots()
{
    if (fPlotter==0) return;
    fPlotter->ShowPlots();
}

Double_t VNeuralNet::TrainEpoch(TDataServe *server, Int_t nEpoch)
{
    double       error;			// squared error collector
    unsigned int classError;		// classification Error
    unsigned int n;			// number of samples

    const Int_t samples = server->GetNumTrnvecs();
    const Int_t tests   = server->GetNumTstvecs();

    for (Int_t epo=0; epo<nEpoch; epo++){
	error = 0.0;
	classError = 0;
	n = 0;
	
	server->MixTrn(); // Shuffle the dataset

	for (Int_t i=0; i<samples; i++){

	    Int_t trnind = i;
	    if (fBalance) trnind = BalancedTrnIndex(server);

	    Float_t *inv  = (Float_t *) server->GetInvecTrn(trnind);
	    Float_t *outv = (Float_t *) server->GetOutvecTrn(trnind);
	    
	    error += Train(inv,outv);
	    n++;
	}

	classError = (UInt_t) TestEpoch(server);
	double percentage = 100. * classError;
	if (tests>0) percentage /= tests; else percentage = 0;
	
	// print training info
	cout << GetNetID() << ": Epoch " << epo << 
	    ", samples " << n << 
	    ", Error " << error << 
	    ", classError " << classError << 
	    " (" << (int)percentage << "%)" << endl;

	// Fill the plots (for a random recall)
	if (fPlotter!=0) {
	    fPlotter->AddTrainGraph(error);
	    fPlotter->AddTestGraph(classError);
	    fPlotter->ShowPlots();
	    fPlotter->Reset();
	}
    }

    if (fShouldSave) Save(); // Store the net

    return error;
}

// Check the network performance

Double_t VNeuralNet::TestEpoch(TDataServe *server)
{
    unsigned int classError = 0;    // classification Error
    const Int_t samples = server->GetNumTstvecs();
    TNeuralNetParameters &parm = GetParameters();

    for (Int_t i=0; i<samples; i++){

	Int_t tstind = i;
	if (fBalance) tstind = BalancedTstIndex(server);

	Float_t *inv  = server->GetInvecTst(tstind);
	Float_t *outv = server->GetOutvecTst(tstind);
	
	// compare network recall with server
	Recall(inv,outv);

	for (int i=0;i<parm.fOutNodes;++i) {
	    Double_t answer = GetOutput()[i];
	    if ((answer>parm.fThreshold && outv[i]<=parm.fThreshold) ||
		(answer<=parm.fThreshold && outv[i]>parm.fThreshold) )
		++classError; // classification ok ?
	}
	
    }

    return classError;
}

Double_t  VNeuralNet::Test(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    Recall(in,out);
    int I;
    Double_t* o = fOut;
    Double_t diff,totalError = 0.0;
    for (I=0;I<fParm.fOutNodes;++I) {
	diff = *out++ - *o++; 
	totalError += diff * diff;
    }
    return totalError;
}

Double_t VNeuralNet::TrainEpoch(const char *file, Int_t nEpoch)
{
    FILE* ftrn=fopen(file,"rb");
    if (ftrn==0) {
	cerr << "Training file does not exist:" << file << endl;
	return 0.0;
    }
    
    Int_t epoch=0;		// epoch counter
    double error;		// squared error collector
    UInt_t classError;		// classification Error

    NNO_INTYPE   *in  = new NNO_INTYPE[fParm.fInNodes];	    // inputvector
    NNO_OUTTYPE  *out = new NNO_OUTTYPE[fParm.fOutNodes];   // outputvector
    
    // begin of training
    do {
	error      = 0.0;
	classError = 0;
	
	while ( fread(in,sizeof(NNO_INTYPE),fParm.fInNodes,ftrn) ) {  // read inputvector
	    fread(out,sizeof(NNO_OUTTYPE),fParm.fOutNodes,ftrn);      // read outputvector
	    Double_t output = out[0];

	    error += Train(in,out);                // perform learnstep
	    
	    // compare network output with 'Out'
	    Double_t *net = GetOutput();
	    Double_t answer = net[0];
	    for (int I=0;I<fParm.fOutNodes;++I) {
		if ((net[I]>fParm.fThreshold && out[I]<=fParm.fThreshold) ||
		    (net[I]<=fParm.fThreshold && out[I]>fParm.fThreshold) )
		    ++classError; // classification ok ?
	    }

	}
	rewind(ftrn);  // epoch completed, rewind filepointer
	++epoch;
	
	// print training info
	cout << GetNetID() << ": Epoch " << epoch << ", Error " << error << ", classError " << classError << endl;

	// Fill the plots (for a random recall)
	if (fPlotter!=0) {
	    fPlotter->AddTrainGraph(error);
	    fPlotter->AddTestGraph(classError);
	    fPlotter->ShowPlots();
	    fPlotter->Reset();
	}
	
    } while (classError>0 && epoch<nEpoch);
    
    fclose(ftrn);

    delete in; delete out;

    if (fShouldSave) Save(); // Store the net
    
    return error;
}

Double_t VNeuralNet::TestEpoch(const char *file)
{
    FILE* ftst=fopen(file,"rb");
    if (ftst==0) {
	cerr << "Test file does not exist:" << file << endl;
	return 0.0;
    }

    unsigned int classError = 0;    // classification Error
    NNO_INTYPE   *in  = new NNO_INTYPE[fParm.fInNodes];	    // inputvector
    NNO_OUTTYPE  *out = new NNO_OUTTYPE[fParm.fOutNodes];   // outputvector
    
    while ( fread(in,sizeof(NNO_INTYPE),fParm.fInNodes,ftst) ) {   // read inputvector
	fread(out,sizeof(NNO_OUTTYPE),fParm.fOutNodes,ftst);		    // read outputvector
	
	// compare network recall with file
	Double_t *net = Recall(in,out);
	for (int I=0;I<fParm.fOutNodes;++I) {
	    if ((net[I]>fParm.fThreshold && out[I]<=fParm.fThreshold) ||
		(net[I]<=fParm.fThreshold && out[I]>fParm.fThreshold) )
		++classError; // classification ok ?
	}

    }
    
    fclose(ftst);

    delete in; delete out;
    
    return classError;
}

UInt_t VNeuralNet::BalancedTrnIndex(TDataServe *server)
{
    static ULong_t ngood=0, nbad=0;
    UInt_t samples = server->GetNumTrnvecs();
    UInt_t index = (UInt_t) (gRandom->Rndm()*samples);
    Float_t *outv = server->GetOutvecTrn(index);

    if (ngood<nbad) 
	while (outv[0]<=0.0) {
	  index++;
	  outv = server->GetOutvecTrn(index%samples);
	}

    if (outv[0]>0.0) ngood++; else nbad++;

    return index%samples;
}

UInt_t VNeuralNet::BalancedTstIndex(TDataServe *server)
{
    static ULong_t ngood=0, nbad=0;
    UInt_t samples = server->GetNumTstvecs();
    UInt_t index = (UInt_t) (gRandom->Rndm()*samples);
    Float_t *outv = server->GetOutvecTst(index);

    if (ngood<nbad) 
	while (outv[0]<=0.0) {
	  index++;
	  outv = server->GetOutvecTst(index%samples);
	}

    if (outv[0]>0.0) ngood++; else nbad++;

    return index%samples;
}

void VNeuralNet::SetMomentumTerm(Double_t f) 
{ 
    fParm.fMu = f;
}

void VNeuralNet::SetFlatSpotElimination(Double_t f) 
{ 
    fParm.fFse = f;
}
