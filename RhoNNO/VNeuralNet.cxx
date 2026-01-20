// VNeuralNet
//
// Base classes for unsupervised and supervised networks
// Partof the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

static const char* NNO_VERSION="2.0ROOT";

#include "VNeuralNet.h"
#include "VNeuralNetPlotter.h"
#include "TDataServe.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

ClassImp(TNeuralNetParameters)

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
: TNamed("NNO","NNO"), fBalance(false), fOwnPlotter(false), fParm(), fPlotter(0)
{
    fShouldSave = false;
    fFile  = 0;
    fOut = 0;
}

VNeuralNet::VNeuralNet(string netID,int innodes,int outnodes,string netFile)
: TNamed(netID.data(),netID.data()), fBalance(false), fOwnPlotter(false), fParm(), fPlotter(0)
{
    fFilename = netFile;
    strncpy(fParm.fNetId,netID.data(),9);
    fShouldSave  = true;
    fFiletype    = FILE_TEXT;
    fParm.fInNodes  = innodes;
    fParm.fOutNodes = outnodes;
    fFile   = 0;
    fOut    = 0;
    if (outnodes>0) {
        fOut = new double[fParm.fOutNodes];
        TestPointer(fOut);
    }
}

VNeuralNet::VNeuralNet(string netFile)
: TNamed(netFile.data(),netFile.data()), fBalance(false), fOwnPlotter(false), fParm(), fPlotter(0)
{
    fFilename   = netFile;
    fShouldSave = false;
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
    fFile = fopen(fFilename.data(),"wb");
    WriteNet();
}

void VNeuralNet::Save(string file)
{
    fFile = fopen(file.data(),"wb");
    WriteNet();
}

void VNeuralNet::WriteNet() 
{
    char ftype[16];
    if (fFiletype==FILE_BINARY) strcpy(ftype,"binary"); else strcpy(ftype,"text");
    if (fFile==0) { cerr << "VNeuralNet::WriteNet:: Could not open for writing " << fFilename << endl; return; }
    fprintf(fFile,"C++  NEURAL NETWORK OBJECTS   VERSION %s\nFiletype %s\n",NNO_VERSION,ftype);
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
    fFile = fopen(fFilename.data(),"rb");
    if (fFile==0) Errorf((char *)"file %s not found",(char *)fFilename.data());
    char ftype[16];
    char Version[16];
    fscanf(fFile,"C++  NEURAL NETWORK OBJECTS   VERSION %s\nFiletype %s\n",Version,ftype);
    if      (!strcmp(ftype,"binary")) fFiletype = FILE_BINARY;
    else if (!strcmp(ftype,"text"))   fFiletype = FILE_TEXT;
    else Errorf((char *)"illegal fileformat: %s",(char *)fFilename.data());
    
    if (fFiletype==FILE_BINARY)
        ReadNetBinary();
    else
        ReadNetText();
    
    fParm.fNetId[4]=0;
    if (strcmp(netID,fParm.fNetId)) {
        fclose(fFile);
        Errorf((char *)"file %s  (incompatible network)\nnetwork ID is %s and should be %s",fFilename.data(),fParm.fNetId,netID);
    }
    
    if (strcmp(Version,NNO_VERSION)) {
        fclose(fFile);
        Errorf((char *)"illegal NNO version number of file %s\nversion number is %s and should be %s",fFilename.data(),Version,(char *)NNO_VERSION);
    }
    
    if (fFiletype==FILE_BINARY)
        ReadBinary();
    else
        ReadText();
    
    fOut = new double[fParm.fOutNodes];
    
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
    return rand() / (RAND_MAX + 1.);
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
        fOwnPlotter = true;
    }
    
    fPlotter->Initialize();
}

void VNeuralNet::FillPlots(double trn, double tst)
{
    if (fPlotter==0) return;
    fPlotter->AddTrainSample(trn,true);
    fPlotter->AddTestSample(tst,true);
}

void VNeuralNet::ShowPlots()
{
    if (fPlotter==0) return;
    fPlotter->ShowPlots();
}

double VNeuralNet::TrainEpoch(TDataServe *server, int nEpoch)
{
    double       error = 0.0;			// squared error collector
    unsigned int classError;		// classification Error
    unsigned int n;			// number of samples
    
    const int samples = server->GetNumTrnvecs();
    const int tests   = server->GetNumTstvecs();
    
    for (int epo=0; epo<nEpoch; epo++){
        error = 0.0;
        classError = 0;
        n = 0;
        
        server->MixTrn(); // Shuffle the dataset
        
        for (int i=0; i<samples; i++){
            
            int trnind = i;
            if (fBalance) trnind = BalancedTrnIndex(server);
            
            float *inv  = (float *) server->GetInvecTrn(trnind);
            float *outv = (float *) server->GetOutvecTrn(trnind);
            
            double de = Train(inv,outv);
	    if (!isnan(de)) error += de;
            n++;
        }
        
        classError = (unsigned int) TestEpoch(server);
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

double VNeuralNet::TestEpoch(TDataServe *server)
{
    unsigned int classError = 0;    // classification Error
    const int samples = server->GetNumTstvecs();
    TNeuralNetParameters &parm = GetParameters();
    
    int i;
    
    for (i=0; i<samples; i++){
        
        int tstind = i;
        if (fBalance) tstind = BalancedTstIndex(server);
        
        float *inv  = server->GetInvecTst(tstind);
        float *outv = server->GetOutvecTst(tstind);
        
        // compare network recall with server
        Recall(inv,outv);
        
        for (int ii=0;ii<parm.fOutNodes;++ii) {
            double answer = GetOutput()[ii];
            if ((answer>parm.fThreshold && outv[ii]<=parm.fThreshold) ||
                (answer<=parm.fThreshold && outv[ii]>parm.fThreshold) )
                ++classError; // classification ok ?
        }
        
    }
    
    return classError;
}

double  VNeuralNet::Test(NNO_INTYPE* in,NNO_OUTTYPE* out) 
{
    Recall(in,out);
    int I;
    double* o = fOut;
    double diff,totalError = 0.0;
    for (I=0;I<fParm.fOutNodes;++I) {
        diff = *out++ - *o++;
        totalError += diff * diff;
    }
    return totalError;
}

double VNeuralNet::TrainEpoch(string file, int nEpoch)
{
    FILE* ftrn=fopen(file.data(),"rb");
    if (ftrn==0) {
        cerr << "Training file does not exist:" << file << endl;
        return 0.0;
    }
    
    int epoch=0;		// epoch counter
    double error;		// squared error collector
    unsigned int classError;		// classification Error
    
    NNO_INTYPE   *in  = new NNO_INTYPE[fParm.fInNodes];	    // inputvector
    NNO_OUTTYPE  *out = new NNO_OUTTYPE[fParm.fOutNodes];   // outputvector
    
    // begin of training
    do {
        error      = 0.0;
        classError = 0;
        
        while ( fread(in,sizeof(NNO_INTYPE),fParm.fInNodes,ftrn) ) {  // read inputvector
            fread(out,sizeof(NNO_OUTTYPE),fParm.fOutNodes,ftrn);      // read outputvector
            //double output = out[0];
            
            error += Train(in,out);                // perform learnstep
            
            // compare network output with 'Out'
            double *net = GetOutput();
            //double answer = net[0];
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
    
    delete[] in; delete[] out;
    
    if (fShouldSave) Save(); // Store the net
    
    return error;
}

double VNeuralNet::TestEpoch(string file)
{
    FILE* ftst=fopen(file.data(),"rb");
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
        double *net = Recall(in,out);
        for (int I=0;I<fParm.fOutNodes;++I) {
            if ((net[I]>fParm.fThreshold && out[I]<=fParm.fThreshold) ||
                (net[I]<=fParm.fThreshold && out[I]>fParm.fThreshold) )
                ++classError; // classification ok ?
        }
        
    }
    
    fclose(ftst);
    
    delete[] in; delete[] out;
    
    return classError;
}

unsigned int VNeuralNet::BalancedTrnIndex(TDataServe *server)
{
    static ULong_t ngood=0, nbad=0;
    unsigned int samples = server->GetNumTrnvecs();
    unsigned int index = (unsigned int) (rand()%samples);
    float *outv = server->GetOutvecTrn(index);
    
    if (ngood<nbad)
        while (outv[0]<=0.0) {
            index++;
            outv = server->GetOutvecTrn(index%samples);
        }
    
    if (outv[0]>0.0) ngood++; else nbad++;
    
    return index%samples;
}

unsigned int VNeuralNet::BalancedTstIndex(TDataServe *server)
{
    static ULong_t ngood=0, nbad=0;
    unsigned int samples = server->GetNumTstvecs();
    unsigned int index = (unsigned int) (rand()%samples);
    float *outv = server->GetOutvecTst(index);
    
    if (ngood<nbad)
        while (outv[0]<=0.0) {
            index++;
            outv = server->GetOutvecTst(index%samples);
        }
    
    if (outv[0]>0.0) ngood++; else nbad++;
    
    return index%samples;
}

void VNeuralNet::SetMomentumTerm(double f) 
{ 
    fParm.fMu = f;
}

void VNeuralNet::SetFlatSpotElimination(double f) 
{ 
    fParm.fFse = f;
}
