#include <math.h>
#include <stdio.h>
#include "TMath.h"

gSystem->Load("libRooFit");
using namespace RooFit;

void ml_fit()
{

int aws;
cout<<"Use RooFit? (0-sim, 1-nao)"<<endl;
cin >> aws;

TFile* file = new TFile("HISTOGRAMS.root");

if(aws==0){ //unbinned fit
TTree* tree = file->Get("RD_el_WW_tree");
RooRealVar mWjj("mWjj","mWjj",2,306);
RooRealVar side_band("side_band","side_band",2,306);
side_band.setRange("low",50,66);
side_band.setRange("high",98,300);

RooArgSet ntupleVarSet(side_band);
RooDataSet data("data","data",tree,ntupleVarSet);
RooArgSet all(mWjj);
RooDataSet all_data("all_data","all_data",tree,all);

RooRealVar ngaus("ngaus","ngaus",0,150);
RooRealVar mgaus("mgaus","mgaus",90,130);
RooRealVar sigaus("sigaus","sigaus",60,90);
RooRealVar nlandau("nlandau","nlandau",0,150);
RooRealVar mlandau("mlandau","mlandau",60,150);
RooRealVar silandau("silandau","silandau",10,20);
RooRealVar a("a","a",-10,40);
RooRealVar b("b","b",0,150);

RooGenericPdf FVjets("FVjets","ngaus*TMath::Gaus(side_band,mgaus,sigaus)+nlandau*TMath::Landau(side_band,mlandau,silandau)+a*pow(side_band,2)+b",
		     RooArgSet(side_band,ngaus,mgaus,sigaus,nlandau,mlandau,silandau,a,b));
FVjets.fitTo(data,Range("low,high"),SumW2Error(kTRUE));

RooPlot* mWjjframe = side_band.frame();
data.plotOn(mWjjframe,Binning(38));
FVjets.plotOn(mWjjframe,LineColor(kGreen),Range(50,300));

FVjets.fitTo(all_data,SumW2Error(kTRUE));
RooPlot* mWjjframe2 = mWjj.frame();
RooHist* residuo = mWjjframe2->residHist();
all_data.plotOn(mWjjframe2,Binning(38),MarkerColor(kBlue));
FVjets.plotOn(mWjjframe,LineColor(kBlue),Range(50,300));

mWjjframe2->Draw();
mWjjframe->Draw("same");

cout<< "=========== Fit Result ==========="<<endl;
cout<< "Media Gauss: " << mgaus.getVal() <<endl;
cout<< "Sigma Gauss: " << sigaus.getVal() <<endl;
cout<< "Media Landau: " << mlandau.getVal() <<endl;
cout<< "Sigma Landau: " << silandau.getVal() <<endl;
cout<< "a: " << a.getVal() <<endl;
cout<< "b: " << b.getVal() <<endl;
cout<< "==================================" <<endl;

}

else if(aws==1){ //binned fit
TF1* FVjets = new TF1("FVjets","[6]*TMath::Gaus(x,[0],[1])+[7]*TMath::Landau(x,[2],[3])+[4]*x*x+[5]",0,300);
FVjets->SetParName(0,"media Gauss"); FVjets->SetParLimits(0,90,130);
FVjets->SetParName(1,"sigma Gauss"); FVjets->SetParLimits(1,60,90);
FVjets->SetParName(2,"media Landau"); FVjets->SetParLimits(2,60,150);
FVjets->SetParName(3,"sigma Landau"); FVjets->SetParLimits(3,10,20);
FVjets->SetParName(4,"cang Parabo"); FVjets->SetParLimits(4,-10,40);
FVjets->SetParName(5,"clin Parabo");
FVjets->SetParName(6,"norma Gauss");
FVjets->SetParName(7,"norma Landau");

TH1D *all = (TH1D*)file->Get("m_{W} - TL");
all->SetMarkerColor(kBlue); all->SetFillColor(0);
TH1D *hist = (TH1D*)file->Get("SB8");
hist->Fit("FVjets","","L",50,300);
all->Fit("FVjets","","L",50,300);
hist->SetMarkerStyle(4);
hist->Draw("E");
all->Draw("same");
}

else cout<<"Argumento invalido!!"<<endl;
}