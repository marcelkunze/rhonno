// TPerceptron
//
// Implementation of the perceptron (Supervised Learning)
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include <iostream>
#include <math.h>
#include <algorithm>  // for std::max, std::min
//#include "TMath.h"
#include "TPerceptron.h"

ClassImp(TPerceptron)

// Transferfunctions - optimized versions
void TransferFermi(double in,double* out,double* deriv) 
{
    // Clamp input to avoid overflow in exp()
    in = std::max(-10.0, std::min(10.0, in));
    double exp_neg_in = exp(-in);
    double O = 1.0 / (1.0 + exp_neg_in);
    *out   = O;
    *deriv = O * (1.0 - O);
}

void TransferSigmoid(double in,double* out,double* deriv) 
{
    // Clamp input to avoid overflow in exp()
    in = std::max(-10.0, std::min(10.0, in));
    double exp_neg_in = exp(-in);
    double O = 1.0 - 2.0 / (1.0 + exp_neg_in);
    *out   = O;
    *deriv = 2.0 * O * (1.0 - O);
}

void TransferLinear(double in,double* out,double* deriv) 
{
    *out = in;
    *deriv = 1.0;
}

void TransferLinearBend(double in,double* out,double* deriv) 
{
    if (in < -1.0) {
        *out = -0.9 + in * 0.1;
        *deriv = 0.1;
    } else
        if (in > 1.0) {
            *out = 0.9 + in * 0.1;
            *deriv = 0.1;
        } else {
            *out = in;
            *deriv = 1.0;
        }
}

void TransferReLU(double in,double* out,double* deriv) 
{
    if (in > 0.0) {
        *out = in;
        *deriv = 1.0;
    } else {
        *out = 0.0;
        *deriv = 0.0;
    }
}

TPerceptron::TPerceptron(int inNodes,
                         int outNodes,
                         double learnStep,
                         TNeuralNetParameters::TRANSFER transferId,
                         int perceptronId)
{
    
    fParm.fInNodes      = inNodes;
    fParm.fOutNodes     = outNodes;
    fParm.fLearnStep    = learnStep;
    fParm.fTransferId   = transferId;
    fParm.fPerceptronId = perceptronId;
    fPrev     = 0;
    AllocNet();
    InitNet();
}

TPerceptron::TPerceptron(TPerceptron* prev,
                         int outNodes,
                         double learnStep,
                         TNeuralNetParameters::TRANSFER transferId,
                         int perceptronId)
{
    
    fParm.fInNodes      = prev->fParm.fOutNodes;
    fParm.fOutNodes     = outNodes;
    fParm.fLearnStep    = learnStep;
    fParm.fTransferId   = transferId;
    fParm.fPerceptronId = perceptronId;
    fPrev             = prev;
    AllocNet();
    InitNet();
}

TPerceptron::TPerceptron(void) 
{
    fPrev = 0;
}

TPerceptron::TPerceptron(TPerceptron* prev) 
{
    fPrev = prev;
}

void TPerceptron::AllocNet(void) 
{
    int I;
    fU       = new PerceptronUnit[fParm.fOutNodes]; TestPointer(fU);
    fOut     = new double[fParm.fOutNodes]; TestPointer(fOut);
    fDiffSrc = new double[fParm.fOutNodes]; TestPointer(fDiffSrc);
    if (fPrev == 0) {
        fIn      = new double[fParm.fInNodes]; TestPointer(fIn);
        fDiffDst = 0;
    } else {
        fIn      = fPrev->fOut;
        fDiffDst = fPrev->fDiffSrc;
    }
    fUbound = &fU[fParm.fOutNodes];
    PerceptronUnit* up = fU;
    for (I=0;I<fParm.fOutNodes;++I) {
        up->fVector = new double[fParm.fInNodes];
        TestPointer(up->fVector);
        up->fDelta = new double[fParm.fInNodes];
        TestPointer(up->fDelta);
        up->fThreshold = 0;
        up->fID = I;
        ++up;
    }
    switch (fParm.fTransferId) {
        case TNeuralNetParameters::TR_FERMI  :      Transfer=TransferFermi;      break;
        case TNeuralNetParameters::TR_SIGMOID:      Transfer=TransferSigmoid;      break;
        case TNeuralNetParameters::TR_LINEAR :      Transfer=TransferLinear;     break;
        case TNeuralNetParameters::TR_LINEAR_BEND:  Transfer=TransferLinearBend; break;
        case TNeuralNetParameters::TR_RELU   :      Transfer=TransferReLU;       break;
        default:             Transfer=0;
    }
}

void TPerceptron::InitNet(void) 
{
    PerceptronUnit* up;
    int J;
    for(up=fU;up<fUbound;++up) {
        for (J=0;J<fParm.fInNodes;++J) {
            up->fVector[J] = Random();
            up->fDelta[J]  = 0.0;
        }
        up->fThreshold = Random();
    }
}


TPerceptron ::~TPerceptron() 
{
    PerceptronUnit* up = fU;
    if (fU!=0)  {
        for(up=fU;up<fUbound;++up) {
            delete[] up->fVector;
            delete[] up->fDelta;
        }
        delete[] fU;
    }
    delete[] fOut;
    delete[] fDiffSrc;
    if (fPrev==0) delete[] fIn;
}


void TPerceptron::WriteBinary() 
{
    PerceptronUnit* up;
    fwrite(&fParm,sizeof(PerceptronBase),1,fFile);
    for(up=fU;up<fUbound;++up) {
        fwrite(up->fVector,sizeof(double),fParm.fInNodes,fFile);
        fwritevar(up->fThreshold);
        fwritevar(up->fID);
    }
}

void TPerceptron::ReadBinary() 
{
    PerceptronUnit* up;
    fread(&fParm,sizeof(PerceptronBase),1,fFile);
    AllocNet();
    for(up=fU;up<fUbound;++up) {
        fread(up->fVector,sizeof(double),fParm.fInNodes,fFile);
        freadvar(up->fThreshold);
        freadvar(up->fID);
    }
}

void TPerceptron::WriteText() 
{
    fprintf(fFile,"\nPerceptron ID %i\n",fParm.fPerceptronId);
    fprintf(fFile,"innodes     %i\n",fParm.fInNodes);
    fprintf(fFile,"outnodes    %i\n",fParm.fOutNodes);
    fprintf(fFile,"learn_step  %le\n",fParm.fLearnStep);
    fprintf(fFile,"transfer_id %i\n",fParm.fTransferId);
    PerceptronUnit* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        fprintf(fFile,"\n");
        fprintf(fFile,"unit number      %i\n",up->fID);
        fprintf(fFile,"threshold        %le\n",up->fThreshold);
        fprintf(fFile,"weights\n");
        for (I=0;I<fParm.fInNodes;++I) fprintf(fFile,"%le\n",up->fVector[I]);
        fprintf(fFile,"\n");
    }
}

void TPerceptron::ReadText() 
{
    fscanf(fFile,"\nPerceptron ID %i\n",&fParm.fPerceptronId);
    fscanf(fFile,"innodes     %i\n",&fParm.fInNodes);
    fscanf(fFile,"outnodes    %i\n",&fParm.fOutNodes);
    fscanf(fFile,"learn_step  %le\n",&fParm.fLearnStep);
    fscanf(fFile,"transfer_id %i\n",(int *)&fParm.fTransferId);
    AllocNet();
    PerceptronUnit* up;
    int I;
    for(up=fU;up<fUbound;++up) {
        fscanf(fFile,"\n");
        fscanf(fFile,"unit number      %i\n",&up->fID);
        fscanf(fFile,"threshold        %le\n",&up->fThreshold);
        fscanf(fFile,"weights\n");
        for (I=0;I<fParm.fInNodes;++I) fscanf(fFile,"%le\n",&up->fVector[I]);
        fscanf(fFile,"\n");
    }
}

double* TPerceptron::Inference(NNO_INTYPE*,NNO_OUTTYPE*) 
{
    if (Transfer==0) Errorf((char *)"(TPerceptron) undefined transferfunction");
    
    // Optimized version using loop unrolling and better memory access
    const int unroll_factor = 4;  // Unroll factor for better ILP
    PerceptronUnit* up = fU;
    double* o = fOut;
    double* ds = fDiffSrc;
    
    // Process units in chunks for better cache performance
    for(; up < fUbound - (unroll_factor - 1); up += unroll_factor, o += unroll_factor, ds += unroll_factor) {
        // Unrolled loop for better instruction-level parallelism
        for(int k = 0; k < unroll_factor; ++k) {
            PerceptronUnit* current_up = up + k;
            double* v = current_up->fVector;
            double* i = fIn;
            double sum = 0.0;
            
            // Main dot product loop - optimized for cache performance
            int I;
            for (I = 0; I < fParm.fInNodes - 3; I += 4) {
                sum += i[I] * v[I] + i[I+1] * v[I+1] + i[I+2] * v[I+2] + i[I+3] * v[I+3];
            }
            // Handle remaining elements
            for(; I < fParm.fInNodes; ++I) {
                sum += i[I] * v[I];
            }
            
            sum -= current_up->fThreshold;
            Transfer(sum, o + k, ds + k);
        }
    }
    
    // Handle remaining units
    for(; up < fUbound; ++up, ++o, ++ds) {
        double* v = up->fVector;
        double* i = fIn;
        double sum = 0.0;
        
        // Optimized dot product with loop unrolling
        int I;
        for (I = 0; I < fParm.fInNodes - 3; I += 4) {
            sum += i[I] * v[I] + i[I+1] * v[I+1] + i[I+2] * v[I+2] + i[I+3] * v[I+3];
        }
        // Handle remaining elements
        for(; I < fParm.fInNodes; ++I) {
            sum += i[I] * v[I];
        }
        
        sum -= up->fThreshold;
        Transfer(sum, o, ds);
    }
    
    return fOut + fParm.fOutNodes;
}

double TPerceptron::Train(NNO_INTYPE*,NNO_OUTTYPE*) 
{
    // modify weights - optimized version
    double* ds = fDiffSrc;
    const double learn_step = fParm.fLearnStep;
    const double mu = fParm.fMu;
    
    for(PerceptronUnit* up = fU; up < fUbound; ++up, ++ds) {
        double* i = fIn;
        double* v = up->fVector;
        double* m = up->fDelta;
        
        // Optimized weight update loop with reduced NaN checks
        int I;
        for (I = 0; I < fParm.fInNodes - 3; I += 4) {
            // Process 4 elements at once for better ILP
            double i0 = i[I], i1 = i[I+1], i2 = i[I+2], i3 = i[I+3];
            double ds_val = *ds;
            
            // Check for NaN only once per group
            if (std::isnan(ds_val)) ds_val = 1.0;
            
            double delta0 = i0 * ds_val;
            double delta1 = i1 * ds_val; 
            double delta2 = i2 * ds_val;
            double delta3 = i3 * ds_val;
            
            v[I]   += (delta0 + m[I]   * mu) * learn_step;
            v[I+1] += (delta1 + m[I+1] * mu) * learn_step;
            v[I+2] += (delta2 + m[I+2] * mu) * learn_step;
            v[I+3] += (delta3 + m[I+3] * mu) * learn_step;
            
            m[I]   = delta0;
            m[I+1] = delta1;
            m[I+2] = delta2;
            m[I+3] = delta3;
        }
        
        // Handle remaining elements
        for(; I < fParm.fInNodes; ++I) {
            double i_val = i[I];
            double ds_val = *ds;
            
            if (std::isnan(ds_val)) ds_val = 1.0;
            if (std::isnan(i_val)) i_val = 1.0;
            
            double delta = i_val * ds_val;
            if (!std::isnan(delta)) {
                v[I] += (delta + m[I] * mu) * learn_step;
                m[I] = delta;
            }
        }
        
        // Update threshold
        double delta_thresh = *ds * learn_step;
        if (!std::isnan(delta_thresh)) {
            up->fThreshold -= delta_thresh;
        }
    }
    
    // propagate derivation backward if previous perceptron exists
    if (fDiffDst != nullptr) {
        double* dd = fDiffDst;
        
        // Pre-compute ds values to avoid repeated access
        for (int I = 0; I < fParm.fInNodes; ++I) {
            double diff = 0.0;
            ds = fDiffSrc;
            
            // Optimized inner loop with better cache locality
            for(PerceptronUnit* up = fU; up < fUbound; ++up, ++ds) {
                double ds_val = *ds;
                double weight_val = up->fVector[I];
                
                if (std::isnan(ds_val)) ds_val = 1.0;
                if (std::isnan(weight_val)) weight_val = 1.0;
                
                double delta = weight_val * ds_val;
                if (!std::isnan(delta)) {
                    diff += delta;
                }
            }
            
            dd[I] *= diff;
        }
    }
    
    return 0.0;
}


