// VNeuralNetPlotter
// TSimpleNeuralNetPlotter
//
// Base classes for network plotters
// Part of the Neural Network Objects package (NNO)
//
// Author List:
// Johannes Steffens, Bochum University
// M.Kunze, Bochum University, 1995

#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TH1.h"
#include "TText.h"
#include "TTree.h"
#include "TBranch.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"

#include "VNeuralNetPlotter.h"

#include <algorithm>
#include <cstdio>

#define NPMAX 1000

using namespace std;

ClassImp(VNeuralNetPlotter)

VNeuralNetPlotter::VNeuralNetPlotter(string name) : TNamed(name.data(),name.data())
{}

void VNeuralNetPlotter::DrawT(Text_t *text, float x, float y, float angle, int color)
{
    TText *tText= new TText(x,y,text);
    tText->SetNDC(true);
    tText->SetTextColor(color);
    tText->SetTextAngle(angle);
    tText->Draw();
}


ClassImp(TSimpleNeuralNetPlotter)

TSimpleNeuralNetPlotter::TSimpleNeuralNetPlotter(string name) 
: VNeuralNetPlotter(name), fCanvas(0),
fTrnCurve(0), fNtrn(0), fTstCurve(0), fNtst(0)
{
    fPlots = false;
    fXtrn = new double[NPMAX];
    fYtrn = new double[NPMAX];
    fXtst = new double[NPMAX];
    fYtst = new double[NPMAX];
    fTrnHistTrue = nullptr;
    fTrnHistFalse = nullptr;
    fTstHistTrue = nullptr;
    fTstHistFalse = nullptr;
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

static TH1D* MakeScoreHist(const string& name, Color_t color)
{
    // Unique name avoids ROOT gDirectory collisions between runs/models
    string uname = name + Form("_%p", (void*)gROOT);
    // Classification scores live in [0,1]; keep small margin for linear nets
    TH1D* h = new TH1D(uname.data(), name.data(), 50, -0.05, 1.05);
    h->SetDirectory(nullptr); // not owned by gDirectory
    h->SetFillColor(color);
    h->SetLineColor(color);
    h->SetLineWidth(2);
    // Hollow-ish hatch: solid fill of two overlapping classes looks like a "box"
    h->SetFillStyle(3004);
    h->SetMarkerColor(color);
    return h;
}

void TSimpleNeuralNetPlotter::Initialize() 
{
    if (fPlots) return; // Nothing to do
    
    // Make sure a canvas exists
    
    if (fCanvas==0) {
        fCanvas = new TCanvas(Form("nno_canvas_%s", GetName()), GetName(), 0, 0, 900, 800);
        fCanvas->SetFillColor(41);
        fCanvas->SetGridx();
        fCanvas->SetGridy();
        fCanvas->GetFrame()->SetFillColor(21);
        fCanvas->GetFrame()->SetBorderSize(12);
        fCanvas->Divide(2,2);
    }
    
    string trnLabel("Train ");
    fTrnPlot = trnLabel + GetName();
    fTrnHistTrue  = MakeScoreHist(fTrnPlot + " signal", kGreen+2);
    fTrnHistFalse = MakeScoreHist(fTrnPlot + " background", kRed+1);
    fTrnHistFalse->SetFillStyle(3005);

    string tstLabel("Test ");
    fTstPlot = tstLabel + GetName();
    fTstHistTrue  = MakeScoreHist(fTstPlot + " signal", kGreen+2);
    fTstHistFalse = MakeScoreHist(fTstPlot + " background", kRed+1);
    fTstHistFalse->SetFillStyle(3005);
    
    fPlots = true;
}

void TSimpleNeuralNetPlotter::AddTrainSample(double trn, bool good)
{
    if (!fTrnHistTrue || !fTrnHistFalse) return;
    if (good)
        fTrnHistTrue->Fill(trn);
    else
        fTrnHistFalse->Fill(trn);
}

void TSimpleNeuralNetPlotter::AddTestSample(double tst, bool good)
{
    if (!fTstHistTrue || !fTstHistFalse) return;
    if (good)
        fTstHistTrue->Fill(tst);
    else
        fTstHistFalse->Fill(tst);
}

void TSimpleNeuralNetPlotter::AddTrainGraph(double trn) 
{
    int n = fNtrn%NPMAX;
    fXtrn[n] = (double) fNtrn + 1;
    fYtrn[n] = trn;
    fNtrn++;
    if (fTrnCurve==0) {
        fTrnCurve = new TGraph;
        fTrnCurve->SetFillColor(19);
        fTrnCurve->SetLineColor(3);
        fTrnCurve->SetLineWidth(1);
        fTrnCurve->SetMarkerColor(3);
        fTrnCurve->SetMarkerStyle(20);
        string histname = fTstPlot + ": Loss function"+";Loss;Epoch";
        fTrnCurve->SetTitle(histname.data());
        fTrnCurve->SetName(histname.data());
    }
}

void TSimpleNeuralNetPlotter::AddTestGraph(double tst) 
{
    int n = fNtst%NPMAX;
    fXtst[n] = (double) fNtst + 1;
    fYtst[n] = tst;
    fNtst++;
    if (fTstCurve==0) {
        fTstCurve = new TGraph;
        fTstCurve->SetFillColor(19);
        fTstCurve->SetLineColor(5);
        fTstCurve->SetLineWidth(1);
        fTstCurve->SetMarkerColor(5);
        fTstCurve->SetMarkerStyle(20);
        string histname = fTstPlot + ": Loss function"+";Loss;Epoch";
        fTstCurve->SetTitle(histname.data());
        fTstCurve->SetName(histname.data());
    }
}

static void DrawScorePair(TH1D* bg, TH1D* sig, const char* title)
{
    if (!bg || !sig) return;
    if (gPad) gPad->Clear();
    bg->SetTitle(title);
    sig->SetTitle(title);
    bg->GetXaxis()->SetTitle("network output");
    bg->GetYaxis()->SetTitle("entries");
    // Common y-scale so both classes stay visible
    const double ymax = std::max(bg->GetMaximum(), sig->GetMaximum()) * 1.15;
    bg->SetMaximum(ymax > 0 ? ymax : 1.0);
    sig->SetMaximum(ymax > 0 ? ymax : 1.0);
    // DrawCopy: pad owns a snapshot; later Reset()/extra TestEpoch fills
    // must not blank the peaks the user is looking at.
    bg->DrawCopy("hist");
    sig->DrawCopy("hist same");
}

void TSimpleNeuralNetPlotter::ShowPlots() 
{
    if (fCanvas==0) return;
    
    if (fTrnHistTrue!=0) {
        fCanvas->cd(1);
        gPad->SetLogy(0);
        DrawScorePair(fTrnHistFalse, fTrnHistTrue, "Training: output score");
        DrawT((char *)"Training",0.12f, 0.92f, 0.f, kBlack);
        DrawT((char *)"red=bg  green=signal",0.12f, 0.87f, 0.f, kBlack);
    }
    
    if (fTstHistTrue!=0) {
        fCanvas->cd(2);
        gPad->SetLogy(0);
        DrawScorePair(fTstHistFalse, fTstHistTrue, "Validation: output score");
        DrawT((char *)"Validation",0.12f,0.92f,0.f, kBlack);
        DrawT((char *)"red=bg  green=signal",0.12f, 0.87f, 0.f, kBlack);
    }
    
    if (fTrnCurve!=0 && fNtrn%NPMAX>0) {
        fCanvas->cd(3);
        gPad->Clear();
        fTrnCurve->DrawGraph(fNtrn%NPMAX,fXtrn,fYtrn,"ALP");
        DrawT((char *)"Training loss",0.2f, 0.8f, 0.f, 3);
    }
    
    if (fTstCurve!=0 && fNtst%NPMAX>0) {
        fCanvas->cd(4);
        gPad->Clear();
        fTstCurve->DrawGraph(fNtst%NPMAX,fXtst,fYtst,"ALP");
        DrawT((char *)"Validation classError",0.2f,0.8f,0.f,5);
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
