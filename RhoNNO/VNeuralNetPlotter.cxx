// VNeuralNetPlotter
// TSimpleNeuralNetPlotter
//
// Base classes for network plotters
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University
// (C) Copyright Johannes Steffens 1995, Ruhr-University Bochum.

#include <TCanvas.h>
#include <TFrame.h>
#include <TGraph.h>
#include <TH1.h>
#include <TText.h>
#include <TTree.h>
#include <TBranch.h>
#include "RhoNNO/VNeuralNetPlotter.h"

#define NPMAX 500

using namespace std;

ClassImp(VNeuralNetPlotter)

VNeuralNetPlotter::VNeuralNetPlotter(string name) : TNamed(name.data(),name.data())
{}

void VNeuralNetPlotter::DrawT(Text_t *text, float x, float y, float angle, int color)
{
    TText *tText= new TText(x,y,text);
    tText->SetNDC(kTRUE);
    tText->SetTextColor(color);
    tText->SetTextAngle(angle);
    tText->Draw();
}


ClassImp(TSimpleNeuralNetPlotter)

TSimpleNeuralNetPlotter::TSimpleNeuralNetPlotter(string name) 
: VNeuralNetPlotter(name), fCanvas(0),
fTrnCurve(0), fNtrn(0), fTstCurve(0), fNtst(0)
{
    fPlots = kFALSE;
    fXtrn = new double[NPMAX];
    fYtrn = new double[NPMAX];
    fXtst = new double[NPMAX];
    fYtst = new double[NPMAX];
}

TSimpleNeuralNetPlotter::~TSimpleNeuralNetPlotter() 
{ 
    if (fTrnHistTrue!=nullptr) { delete fTrnHistTrue; fTrnHistTrue=nullptr;}
    if (fTrnHistFalse!=nullptr) { delete fTrnHistFalse; fTrnHistFalse=nullptr;}
    if (fTstHistTrue!=nullptr) { delete fTstHistTrue; fTstHistTrue=nullptr;}
    if (fTstHistFalse!=nullptr) { delete fTstHistFalse; fTstHistFalse=nullptr;}
    if (fTrnCurve!=nullptr) { delete fTrnCurve; fTrnCurve=nullptr;}
    if (fTstCurve!=nullptr) { delete fTstCurve; fTstCurve=nullptr;}
    if (fXtrn!=nullptr) { delete [] fXtrn; fXtrn=nullptr;}
    if (fYtrn!=nullptr) { delete [] fYtrn; fYtrn=nullptr;}
    if (fXtst!=nullptr) { delete [] fXtst; fXtst=nullptr;}
    if (fYtst!=nullptr) { delete [] fYtst; fYtst=nullptr;}
    if (fCanvas!=nullptr) { fCanvas->Close(); /*delete fCanvas; fCanvas=nullptr;*/}
}

void TSimpleNeuralNetPlotter::Initialize() 
{
    if (fPlots) return; // Nothing to do
    
    // Make sure a canvas exists
    
    if (fCanvas==0) {
        fCanvas = new TCanvas("fCanvas",GetName(),0,0,800,800);
        fCanvas->SetFillColor(41);
        fCanvas->SetGridx();
        fCanvas->SetGridy();
        fCanvas->GetFrame()->SetFillColor(21);
        fCanvas->GetFrame()->SetBorderSize(12);
        fCanvas->Divide(2,2);
    }
    
    string trnLabel("Train");
    fTrnPlot = trnLabel + GetName();
    string histname = fTrnPlot + "+";
    fTrnHistTrue = new TH1D(histname.data(),histname.data(),100,-1.1,1.1);
    fTrnHistTrue->SetFillColor(kGreen);
    histname = fTrnPlot + "-";
    fTrnHistFalse = new TH1D(histname.data(),histname.data(),100,-1.1,1.1);
    fTrnHistFalse->SetFillColor(kRed);
    string tstLabel("Recall");
    fTstPlot = tstLabel + GetName();
    histname = fTstPlot + "+";
    fTstHistTrue = new TH1D(histname.data(),histname.data(),100,-1.1,1.1);
    fTstHistTrue->SetFillColor(kGreen);
    histname = fTstPlot + "-";
    fTstHistFalse = new TH1D(histname.data(),histname.data(),100,-1.1,1.1);
    fTstHistFalse->SetFillColor(kRed);
    
    fPlots = kTRUE;
}

void TSimpleNeuralNetPlotter::AddTrainSample(double trn, bool good)
{
    if (good)
        fTrnHistTrue->Fill(trn);
    else
        fTrnHistFalse->Fill(trn);
}

void TSimpleNeuralNetPlotter::AddTestSample(double tst, bool good)
{
    if (good)
        fTstHistTrue->Fill(tst);
    else
        fTstHistFalse->Fill(tst);
}

void TSimpleNeuralNetPlotter::AddTrainGraph(double trn) 
{
    int n = fNtrn%NPMAX;
    fXtrn[n] = (double) fNtrn;
    fYtrn[n] = trn;
    fNtrn++;
    if (fTrnCurve==0) {
        fTrnCurve = new TGraph;
        fTrnCurve->SetFillColor(19);
        fTrnCurve->SetLineColor(3);
        fTrnCurve->SetLineWidth(3);
        fTrnCurve->SetMarkerColor(3);
        fTrnCurve->SetMarkerStyle(20);
        string histname = fTstPlot + ": Error on training data";
        fTrnCurve->SetTitle(histname.data());
    }
}

void TSimpleNeuralNetPlotter::AddTestGraph(double tst) 
{
    int n = fNtst%NPMAX;
    fXtst[n] = (double) fNtst;
    fYtst[n] = tst;
    fNtst++;
    if (fTstCurve==0) {
        fTstCurve = new TGraph;
        fTstCurve->SetFillColor(19);
        fTstCurve->SetLineColor(5);
        fTstCurve->SetLineWidth(3);
        fTstCurve->SetMarkerColor(5);
        fTstCurve->SetMarkerStyle(20);
        string histname = fTstPlot + ": Error on test data";
        fTstCurve->SetTitle(histname.data());
    }
}

void TSimpleNeuralNetPlotter::ShowPlots() 
{
    if (fCanvas==0) return;
    
    if (fTrnHistTrue!=0) {
        fCanvas->cd(1);
        fTrnHistFalse->Draw();
        fTrnHistTrue->Draw("same");
        DrawT((char *)"Training",0.7f, 0.8f, 0.f, 3);
    }
    
    if (fTstHistTrue!=0) {
        fCanvas->cd(2);
        fTstHistFalse->Draw();
        fTstHistTrue->Draw("same");
        DrawT((char *)"Test",0.7f,0.8f,0.f,5);
    }
    
    if (fTrnCurve!=0 && fNtrn%NPMAX>0) {
        fCanvas->cd(3);
        gPad->Clear();
        fTrnCurve->DrawGraph(fNtrn%NPMAX,fXtrn,fYtrn,"ALP");
        DrawT((char *)"Training",0.7f, 0.8f, 0.f, 3);
    }
    
    if (fTstCurve!=0 && fNtst%NPMAX>0) {
        fCanvas->cd(4);
        gPad->Clear();
        fTstCurve->DrawGraph(fNtst%NPMAX,fXtst,fYtst,"ALP");
        DrawT((char *)"Test",0.7f,0.8f,0.f,5);
    }
    
    fCanvas->Modified();
    fCanvas->Update();
}

void TSimpleNeuralNetPlotter::Reset() 
{
    if (fTrnHistTrue!=0) {
        fTrnHistTrue->Reset();
    }
    
    if (fTrnHistFalse!=0) {
        fTrnHistFalse->Reset();
    }
    
    if (fTstHistTrue!=0) {
        fTstHistTrue->Reset();
    }
    
    if (fTstHistFalse!=0) {
        fTstHistFalse->Reset();
    }
}
