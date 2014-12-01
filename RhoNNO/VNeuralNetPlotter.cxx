//////////////////////////////////////////////////////////////////////////
//									//
// VNeuralNetPlotter							//
// TSimpleNeuralNetPlotter						//
//									//
// Base classes for network plotters					//
// Part of the Neural Network Objects package (NNO)			//
//									//
// Author List:								//
// M.Kunze, Bochum University						//
// (C) Copyright Marcel Kunze 2001, Ruhr-University Bochum.		//
//									//
//////////////////////////////////////////////////////////////////////////

#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TH1.h"
#include "TText.h"
#include "TTree.h"
#include "TBranch.h"
#include "RhoNNO/VNeuralNetPlotter.h"

#define NPMAX 500

ClassImp(VNeuralNetPlotter)

VNeuralNetPlotter::VNeuralNetPlotter(const char *name) : TNamed(name,name)
{}

void VNeuralNetPlotter::DrawT(Text_t *text, Float_t x, Float_t y, Float_t angle, Int_t color)
{
    TText *tText= new TText(x,y,text);
    tText->SetNDC(kTRUE);
    tText->SetTextColor(color);
    tText->SetTextAngle(angle);
    tText->Draw();
}


ClassImp(TSimpleNeuralNetPlotter)

TSimpleNeuralNetPlotter::TSimpleNeuralNetPlotter(const char *name) 
  : VNeuralNetPlotter(name), 
    fCanvas(0), fTrnCurve(0), fNtrn(0), fTstCurve(0), fNtst(0)
{
    fPlots = kFALSE; 
    fXtrn = new Double_t[NPMAX];
    fYtrn = new Double_t[NPMAX];
    fXtst = new Double_t[NPMAX];
    fYtst = new Double_t[NPMAX];
}

TSimpleNeuralNetPlotter::~TSimpleNeuralNetPlotter() 
{ 
    if (fCanvas!=0) { fCanvas->Close(); delete fCanvas; }
    if (fTrnHistTrue!=0) { delete fTrnHistTrue; }
    if (fTrnHistFalse!=0) { delete fTrnHistFalse; }
    if (fTstHistTrue!=0) { delete fTstHistTrue; }
    if (fTstHistFalse!=0) { delete fTstHistFalse; }
    if (fTrnCurve!=0) { delete fTrnCurve; }
    if (fTstCurve!=0) { delete fTstCurve; }
    delete [] fXtrn;
    delete [] fYtrn;
    delete [] fXtst;
    delete [] fYtst;
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

    TString trnLabel("Train");
    fTrnPlot = trnLabel + GetName(); 
    fTrnHistTrue = new TH1D(fTrnPlot,fTrnPlot,100,-1.1,1.1);
    fTrnHistTrue->SetFillColor(kGreen);
    fTrnHistFalse = new TH1D(fTrnPlot,fTrnPlot,100,-1.1,1.1);
    fTrnHistFalse->SetFillColor(kRed);
    TString tstLabel("Recall");
    fTstPlot = tstLabel + GetName(); 
    fTstHistTrue = new TH1D(fTstPlot,fTstPlot,100,-1.1,1.1);
    fTstHistTrue->SetFillColor(kGreen);
    fTstHistFalse = new TH1D(fTstPlot,fTstPlot,100,-1.1,1.1);
    fTstHistFalse->SetFillColor(kRed);

    fPlots = kTRUE;
}

void TSimpleNeuralNetPlotter::AddTrainSample(Double_t trn, Bool_t good)
{
    if (good)
	fTrnHistTrue->Fill(trn);
    else
	fTrnHistFalse->Fill(trn);
}

void TSimpleNeuralNetPlotter::AddTestSample(Double_t tst, Bool_t good)
{
    if (good)
	fTstHistTrue->Fill(tst);
    else
	fTstHistFalse->Fill(tst);
}

void TSimpleNeuralNetPlotter::AddTrainGraph(Double_t trn) 
{
    Int_t n = fNtrn%NPMAX;
    fXtrn[n] = (Double_t) fNtrn;
    fYtrn[n] = trn;
    fNtrn++;
    if (fTrnCurve==0) {
	fTrnCurve = new TGraph;
	fTrnCurve->SetFillColor(19);
	fTrnCurve->SetLineColor(3);
	fTrnCurve->SetLineWidth(3);
	fTrnCurve->SetMarkerColor(3);
	fTrnCurve->SetMarkerStyle(20);
	fTrnCurve->SetTitle(fTrnPlot+": Error on training data");
    }
}

void TSimpleNeuralNetPlotter::AddTestGraph(Double_t tst) 
{
    Int_t n = fNtst%NPMAX;
    fXtst[n] = (Double_t) fNtst;
    fYtst[n] = tst;
    fNtst++;
    if (fTstCurve==0) {
	fTstCurve = new TGraph;
	fTstCurve->SetFillColor(19);
	fTstCurve->SetLineColor(5);
	fTstCurve->SetLineWidth(3);
	fTstCurve->SetMarkerColor(5);
	fTstCurve->SetMarkerStyle(20);
	fTstCurve->SetTitle(fTstPlot+": Error on test data");
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
