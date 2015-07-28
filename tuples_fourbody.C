#include "../weights.h"
#include "TMath.h"
#include <vector>
#include <iostream>


void tuples_fourbody() {

//------------------------ Defining TFile ----------------------
TFile* file_wz     = TFile::Open("histograms/histogram_wz350.root");
TFile* file_ww     = TFile::Open("histograms/histogram_ww350.root");
TFile* file_tt     = TFile::Open("histograms/histogram_tt350.root");

TFile* file_stopstbar  = TFile::Open("histograms/histogram_stopstbar350.root");
TFile* file_stopst     = TFile::Open("histograms/histogram_stopst350.root");
TFile* file_stoptwtbar = TFile::Open("histograms/histogram_stoptwtbar350.root");
TFile* file_stoptwt    = TFile::Open("histograms/histogram_stoptwt350.root");
TFile* file_stopttbar  = TFile::Open("histograms/histogram_stopttbar350.root");
TFile* file_stoptt     = TFile::Open("histograms/histogram_stoptt350.root");

TFile* file_h170   = TFile::Open("histograms/histogram_350.root");
TFile* file_data   = TFile::Open("histograms/histogram_data350.root");

//--------------------- Defining histograms SB -----------------------------
TH1D *hist_wz = (TH1D*)file_wz->Get("MlvjjinSB"); hist_wz->Scale(weight_wz);
TH1D *hist_ww = (TH1D*)file_ww->Get("MlvjjinSB"); hist_ww->Scale(weight_ww);
TH1D *hist_tt = (TH1D*)file_tt->Get("MlvjjinSB"); hist_tt->Scale(weight_tt);

TH1D *hist_stopstbar  = (TH1D*)file_stopstbar->Get("MlvjjinSB");  hist_stopstbar->Scale(weight_stopstbar);
TH1D *hist_stopst     = (TH1D*)file_stopst->Get("MlvjjinSB");     hist_stopst->Scale(weight_stopst);
TH1D *hist_stoptwtbar = (TH1D*)file_stoptwtbar->Get("MlvjjinSB"); hist_stoptwtbar->Scale(weight_stoptwtbar);
TH1D *hist_stoptwt    = (TH1D*)file_stoptwt->Get("MlvjjinSB");    hist_stoptwt->Scale(weight_stoptwt);
TH1D *hist_stopttbar  = (TH1D*)file_stopttbar->Get("MlvjjinSB");  hist_stopttbar->Scale(weight_stopttbar);
TH1D *hist_stoptt     = (TH1D*)file_stoptt->Get("MlvjjinSB");     hist_stoptt->Scale(weight_stoptt);

TH1D *hist_data = (TH1D*)file_data->Get("MlvjjinSB");

//----------------------- Wjets = Data - MC_VV-TTbar -----------------------
int sc=0, xi_fit=200, xf_fit=700, bin=25, nbin=120, xf=3000;
TH1D *sideband= new TH1D("sideband","",nbin,0,xf);
for(int i=0; i<nbin; i++){
  float vdata = hist_data->GetBinContent(i+1);
  float vwz   = hist_wz->GetBinContent(i+1);
  float vww   = hist_ww->GetBinContent(i+1);
  float vtt   = hist_tt->GetBinContent(i+1);
  float vstopstbar = hist_stopstbar->GetBinContent(i+1);
  float vstopst = hist_stopst->GetBinContent(i+1);
  float vstoptwtbar = hist_stoptwtbar->GetBinContent(i+1);
  float vstopstwt = hist_stoptwt->GetBinContent(i+1);
  float vstopsttbar = hist_stopttbar->GetBinContent(i+1);
  float vstopstt = hist_stoptt->GetBinContent(i+1);
  
  float value = vdata-vwz-vww-vtt-vstopstbar-vstopst-vstoptwtbar-vstopstwt-vstopsttbar-vstopstt;
  sideband->SetBinContent(i+1,value);
}

double N4body, N2body=56807.8, scale;
cout<<"Original 4body Normalization:  "<<sideband->Integral(1,nbin)<<endl;
for(int sc=1; sc<200000; sc++){
    scale = sc/100000.;
    TH1D *temp = new TH1D("temp","",nbin,0,xf);
    temp->Add(sideband);
    temp->Scale(scale);
    N4body = temp->Integral(1,nbin);
    if(fabs(N4body-N2body)<0.5) break;
    temp->Delete();  
}

sideband->Scale(scale);
cout<<"4body Renormalization:         "<<N4body<<endl;
cout<<"Scale used:		      "<<scale<<endl;
    
//------------------ Soma dos MC com THStack (MlvjjinSB) ----------
TH1D *top_stack = new TH1D("top_stack","",nbin,0,xf); 
top_stack->SetFillColor(kGreen);
top_stack->SetLineColor(kGreen);
for(int i=0; i<nbin; i++){
  float vtt   = hist_tt->GetBinContent(i+1);
  float vstopstbar = hist_stopstbar->GetBinContent(i+1);
  float vstopst = hist_stopst->GetBinContent(i+1);
  float vstoptwtbar = hist_stoptwtbar->GetBinContent(i+1);
  float vstopstwt = hist_stoptwt->GetBinContent(i+1);
  float vstopsttbar = hist_stopttbar->GetBinContent(i+1);
  float vstopstt = hist_stoptt->GetBinContent(i+1);
  
  float value = vtt+vstopstbar+vstopst+vstoptwtbar+vstopstwt+vstopsttbar+vstopstt;
  top_stack->SetBinContent(i+1,value);
}

TH1D *w2jsb_stack = new TH1D("w2jsb_stack","",nbin,0,xf);
w2jsb_stack->SetFillColor(kRed);
w2jsb_stack->SetLineColor(kRed);
for(int i=0; i<nbin; i++){
  float vtt   = hist_tt->GetBinContent(i+1);
  float vstopstbar = hist_stopstbar->GetBinContent(i+1);
  float vstopst = hist_stopst->GetBinContent(i+1);
  float vstoptwtbar = hist_stoptwtbar->GetBinContent(i+1);
  float vstopstwt = hist_stoptwt->GetBinContent(i+1);
  float vstopsttbar = hist_stopttbar->GetBinContent(i+1);
  float vstopstt = hist_stoptt->GetBinContent(i+1);
  float vsideband = sideband->GetBinContent(i+1);
  
  float value = vtt+vstopstbar+vstopst+vstoptwtbar+vstopstwt+vstopsttbar+vstopstt+vsideband;
  w2jsb_stack->SetBinContent(i+1,value);
}

TH1D *vv_stack = new TH1D("vv_stack","",nbin,0,xf);
vv_stack->SetFillColor(kCyan);
vv_stack->SetLineColor(kCyan);
for(int i=0; i<nbin; i++){
  float vwz   = hist_wz->GetBinContent(i+1);
  float vww   = hist_ww->GetBinContent(i+1);
  float vtt   = hist_tt->GetBinContent(i+1);
  float vstopstbar = hist_stopstbar->GetBinContent(i+1);
  float vstopst = hist_stopst->GetBinContent(i+1);
  float vstoptwtbar = hist_stoptwtbar->GetBinContent(i+1);
  float vstopstwt = hist_stoptwt->GetBinContent(i+1);
  float vstopsttbar = hist_stopttbar->GetBinContent(i+1);
  float vstopstt = hist_stoptt->GetBinContent(i+1);
  float vsideband = sideband->GetBinContent(i+1);
  
  float value = vwz+vww+vtt+vstopstbar+vstopst+vstoptwtbar+vstopstwt+vstopsttbar+vstopstt+vsideband;
  vv_stack->SetBinContent(i+1,value);
}


TCanvas *a = new TCanvas("a","THStack");
vv_stack->Draw();
w2jsb_stack->Draw("same");
top_stack->Draw("same");
hist_data->Draw("E1,same"); hist_data->SetMarkerStyle(20); hist_data->SetMarkerColor(kBlack); hist_data-> SetLineColor(kBlack); hist_data->SetLineWidth(2);
TH1D *hist_h170 = (TH1D*)file_h170->Get("MlvjjinSB"); hist_h170->Scale(weight_h350); hist_h170->SetFillColor(0); hist_h170->SetLineColor(kBlue); hist_h170->SetLineWidth(2);
//for(int i=0; i<3; i++) hist_h170->Add(hist_h170); hist_h170->Draw("same");

TLegend *leg1 = new TLegend(0.70,0.58,0.85,0.81);
leg1->AddEntry(vv_stack,"WW+WZ","f");
leg1->AddEntry(w2jsb_stack,"W2Jets","f");
leg1->AddEntry(top_stack,"Top","f");
leg1->AddEntry(hist_data,"Data","P");
leg1->AddEntry(hist_h170,"H170x4","L");
leg1->SetBorderSize(0); leg1->SetFillColor(0);leg1->SetTextSize(0.03);
leg1->Draw();


//For the limit
int nbin2=(xf_fit-xi_fit)/bin, ri=xi_fit/bin;
TH1D *data = new TH1D("data","",nbin2,xi_fit,xf_fit);
TH1D *signal = new TH1D("signal","",nbin2,xi_fit,xf_fit);
TH1D *background = new TH1D("background","",nbin2,xi_fit,xf_fit);
for(int m=0; m<nbin2; m++){
  float vd = hist_data->GetBinContent(m+ri);
  data->SetBinContent(m+1,vd);
  float vs = hist_h170->GetBinContent(m+ri);
  signal->SetBinContent(m+1,vs);
  float vb = vv_stack->GetBinContent(m+ri);
  background->SetBinContent(m+1,vb); 
}
TFile *hists_to_limit170 = TFile::Open("hists_to_limit350.root","recreate");
data->Write();
signal->Write();
background->Write();
hists_to_limit170->Close();


//Making the input file to combine (limit extraction)
FILE *tolimit;
tolimit = fopen("limits350.txt","w");
fprintf(tolimit,"%s\t%i\t%s\n","imax",nbin2,"number of channels");
fprintf(tolimit,"%s\t%i\t%s\n","jmax",1,"number of backgrounds");
fprintf(tolimit,"%s\t%i\t%s\n","kmax",2,"number of nuisance parameters");
fprintf(tolimit,"%s\n","------------------------------------------------------------------------------");
fprintf(tolimit,"%s\t\t\t","bin"); for(int t=0; t<nbin2; t++) fprintf(tolimit,"%i\t",t+1);
fprintf(tolimit,"\n%s\t\t","observation"); for(int t=0; t<nbin2; t++) fprintf(tolimit,"%i\t",data->GetBinContent(t+1));
fprintf(tolimit,"\n%s","------------------------------------------------------------------------------");
fprintf(tolimit,"\n%s\t\t\t","bin"); for(int t=0; t<nbin2; t++) fprintf(tolimit,"%i\t%i\t",t+1,t+1);
fprintf(tolimit,"\n%s\t\t\t","process"); for(int t=0; t<nbin2; t++) fprintf(tolimit,"%s\t%s\t","ggH","others");
fprintf(tolimit,"\n%s\t\t\t","process"); for(int t=0; t<nbin2; t++) fprintf(tolimit,"%i\t%i\t",0,1);
fprintf(tolimit,"\n%s\t\t\t","rate"); for(int t=0; t<nbin2; t++) fprintf(tolimit,"%.2f\t%.2f\t",signal->GetBinContent(t+1),background->GetBinContent(t+1));
fprintf(tolimit,"\n%s\n","------------------------------------------------------------------------------");
fprintf(tolimit,"%s\t%s\t","CMS_eff_e","lnN"); for(int t=0; t<2*nbin2; t++) fprintf(tolimit,"%.2f\t",1.046);
fprintf(tolimit,"\n%s\t%s\t","lumi_8TeV","lnN"); for(int t=0; t<2*nbin2; t++) fprintf(tolimit,"%.2f\t",1.026);
fclose(tolimit);

}
