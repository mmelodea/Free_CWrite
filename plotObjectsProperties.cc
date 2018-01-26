#include <iostream>
#include <fstream>
#include <vector>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TArrow.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TSystem.h>


const unsigned int nhists = 51;


void loopSamples(std::vector<TString> tags, std::vector<std::vector<TString>> samples){
  ofstream outfile;
  //outfile.open("/home/micah/temp/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_jets.h");
  //outfile.open("/home/micah/temp/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8_jets.h");
  outfile.open("/home/micah/temp/ZZTo4L_13TeV_powheg_pythia8_jets.h");
  
  Float_t f_weight, f_int_weight, f_pu_weight, f_eff_weight, f_lept1_pt, f_lept1_eta, f_lept1_phi, f_lept1_charge, f_lept1_pfx, f_lept1_sip, f_lept1_mvaid, f_lept2_pt, f_lept2_eta, f_lept2_phi, f_lept2_charge, f_lept2_pfx, f_lept2_sip, f_lept2_mvaid, f_lept3_pt, f_lept3_eta, f_lept3_phi, f_lept3_charge, f_lept3_pfx, f_lept3_sip, f_lept3_mvaid, f_lept4_pt, f_lept4_eta, f_lept4_phi, f_lept4_charge, f_lept4_pfx, f_lept4_sip, f_lept4_mvaid, f_iso_max, f_sip_max, f_Z1mass, f_Z2mass, f_angle_costhetastar, f_angle_costheta1, f_angle_costheta2, f_angle_phi, f_angle_phistar1, f_eta4l, f_pt4l, f_mass4l, f_mass4lErr, f_njets_pass, f_deltajj, f_massjj, f_D_jet, f_jet1_pt, f_jet1_eta, f_jet1_phi, f_jet1_e, f_jet2_pt, f_jet2_eta, f_jet2_phi, f_jet2_e;
  Float_t f_jets_dnn_pt[6], f_jets_dnn_eta[6], f_jets_dnn_phi[6], f_jets_dnn_e[6];
  Float_t f_D_bkg_kin,f_D_bkg,f_D_gg,f_D_g4,f_Djet_VAJHU; 
  Float_t f_genmet, f_pfmet,f_mT,f_dphi;
  Int_t f_lept1_pdgid,f_lept2_pdgid,f_lept3_pdgid,f_lept4_pdgid;
  Int_t f_category,f_Ngood, f_Nbjets;
  Int_t f_run, f_lumi, f_event;
  
  std::vector<TString> histos_name;
  histos_name.push_back("lept1_pt");
  histos_name.push_back("lept1_eta");
  histos_name.push_back("lept1_phi");
  histos_name.push_back("lept1_pfx");
  histos_name.push_back("lept1_sip");
  histos_name.push_back("lept2_pt");
  histos_name.push_back("lept2_eta");
  histos_name.push_back("lept2_phi");
  histos_name.push_back("lept2_pfx");
  histos_name.push_back("lept2_sip");
  histos_name.push_back("lept3_pt");
  histos_name.push_back("lept3_eta");
  histos_name.push_back("lept3_phi");
  histos_name.push_back("lept3_pfx");
  histos_name.push_back("lept3_sip");
  histos_name.push_back("lept3_pt");
  histos_name.push_back("lept3_eta");
  histos_name.push_back("lept3_phi");
  histos_name.push_back("lept3_pfx");
  histos_name.push_back("lept3_sip");
  histos_name.push_back("Z1mass");
  histos_name.push_back("Z2mass");
  histos_name.push_back("costhetastar");
  histos_name.push_back("costheta1");
  histos_name.push_back("costheta2");
  histos_name.push_back("phi");
  histos_name.push_back("phistar1");
  histos_name.push_back("pt4l");
  histos_name.push_back("eta4l");
  histos_name.push_back("mass4l");
  histos_name.push_back("njets_pass");
  histos_name.push_back("deltajj");
  histos_name.push_back("massjj");
  histos_name.push_back("D_jet");
  histos_name.push_back("jet1_pt");
  histos_name.push_back("jet1_eta");
  histos_name.push_back("jet1_phi");
  histos_name.push_back("jet1_e");
  histos_name.push_back("jet2_pt");
  histos_name.push_back("jet2_eta");
  histos_name.push_back("jet2_phi");
  histos_name.push_back("jet2_e");
  histos_name.push_back("jet3_pt");
  histos_name.push_back("jet3_eta");
  histos_name.push_back("jet3_phi");
  histos_name.push_back("jet3_e");
  histos_name.push_back("D_bkg_kin");
  histos_name.push_back("D_bkg");
  histos_name.push_back("D_gg");
  histos_name.push_back("D_g4");
  histos_name.push_back("Djet_VAJHU");
  
    
  const unsigned int ntags = tags.size();
  TH1D *histos[ntags][nhists];
  for(unsigned int itag=0; itag<ntags; ++itag){
  /*  
  histos[itag][0] = new TH1D("f_lept1_pt","",100,0,250);
  histos[itag][1] = new TH1D("f_lept1_eta","",50,2.5,2.5);
  histos[itag][2] = new TH1D("f_lept1_phi","",64,-3.2,3.2);
  histos[itag][3] = new TH1D("f_lept1_pfx","",100,0,0.4);
  histos[itag][4] = new TH1D("f_lept1_sip","",100,-5,5);
  histos[itag][5] = new TH1D("f_lept2_pt","",100,0,250);
  histos[itag][6] = new TH1D("f_lept2_eta","",50,2.5,2.5);
  histos[itag][7] = new TH1D("f_lept2_phi","",64,-3.2,3.2);
  histos[itag][8] = new TH1D("f_lept2_pfx","",100,0,0.4);
  histos[itag][9] = new TH1D("f_lept2_sip","",100,-5,5);
  histos[itag][10] = new TH1D("f_lept3_pt","",100,0,250);
  histos[itag][11] = new TH1D("f_lept3_eta","",50,2.5,2.5);
  histos[itag][12] = new TH1D("f_lept3_phi","",64,-3.2,3.2);
  histos[itag][13] = new TH1D("f_lept3_pfx","",100,0,0.4);
  histos[itag][14] = new TH1D("f_lept3_sip","",100,-5,5);
  histos[itag][15] = new TH1D("f_lept4_pt","",100,0,250);
  histos[itag][16] = new TH1D("f_lept4_eta","",50,2.5,2.5);
  histos[itag][17] = new TH1D("f_lept4_phi","",64,-3.2,3.2);
  histos[itag][18] = new TH1D("f_lept4_pfx","",100,0,0.4);
  histos[itag][19] = new TH1D("f_lept4_sip","",100,-5,5);
  histos[itag][20] = new TH1D("f_Z1mass","",54,12,120);
  histos[itag][21] = new TH1D("f_Z2mass","",54,12,120);
  histos[itag][22] = new TH1D("f_angle_costhetastar","",50,-1,1);
  histos[itag][23] = new TH1D("f_angle_costheta1","",50,-1,1);
  histos[itag][24] = new TH1D("f_angle_costheta2","",50,-1,1);
  histos[itag][25] = new TH1D("f_angle_phi","",50,-1,1);
  histos[itag][26] = new TH1D("f_angle_phistar1","",50,-1,1);
  histos[itag][27] = new TH1D("f_pt4l","",250,0,500);
  histos[itag][28] = new TH1D("f_eta4l","",100,-5,5);
  histos[itag][29] = new TH1D("f_mass4l","",110,80,300);
  histos[itag][30] = new TH1D("f_njets_pass","",10,0,10);
  histos[itag][31] = new TH1D("f_deltajj","",100,0,8);
  histos[itag][32] = new TH1D("f_massjj","",150,0,1500);
  histos[itag][33] = new TH1D("f_D_jet","",100,0,3);
  histos[itag][34] = new TH1D("f_jet1_pt","",175,0,350);
  histos[itag][35] = new TH1D("f_jet1_eta","",50,-5,5);
  histos[itag][36] = new TH1D("f_jet1_phi","",64,-3.2,3.2);
  histos[itag][37] = new TH1D("f_jet1_e","",100,0,500);
  histos[itag][38] = new TH1D("f_jet2_pt","",175,0,350);
  histos[itag][39] = new TH1D("f_jet2_eta","",50,-5,5);
  histos[itag][40] = new TH1D("f_jet2_phi","",64,-3.2,3.2);
  histos[itag][41] = new TH1D("f_jet2_e","",100,0,500);
  histos[itag][42] = new TH1D("f_jet3_pt","",175,0,350);
  histos[itag][43] = new TH1D("f_jet3_eta","",50,-5,5);
  histos[itag][44] = new TH1D("f_jet3_phi","",64,-3.2,3.2);
  histos[itag][45] = new TH1D("f_jet3_e","",100,0,500);
  histos[itag][46] = new TH1D("f_D_bkg_kin","",60,0,1);
  histos[itag][47] = new TH1D("f_D_bkg","",60,0,1);
  histos[itag][48] = new TH1D("f_D_gg","",60,0,1);
  histos[itag][49] = new TH1D("f_D_g4","",60,0,1);
  histos[itag][50] = new TH1D("f_Djet_VAJHU","",60,0,1);
  */  
  unsigned int nsamples = samples[itag].size();
  for(unsigned int is=0; is<nsamples; ++is){
    TString ifile_name_root = "OrganizedSamples/"+samples[itag][is];
    if(gSystem->AccessPathName(ifile_name_root)){
      cout<<"File "<<ifile_name_root<<" doesn't exist!"<<endl;
      continue;
    }
    else cout<<"Loading file "<<ifile_name_root<<endl;
  
    TFile *ofile = TFile::Open(ifile_name_root);
    TTree *otree = (TTree*)ofile->Get("HZZ4LeptonsAnalysisReduced");
    otree->SetBranchAddress("f_run", &f_run);
    otree->SetBranchAddress("f_lumi", &f_lumi);    
    otree->SetBranchAddress("f_event", &f_event);    
    otree->SetBranchAddress("f_weight", &f_weight);
    otree->SetBranchAddress("f_int_weight", &f_int_weight);
    otree->SetBranchAddress("f_pu_weight", &f_pu_weight);
    otree->SetBranchAddress("f_eff_weight", &f_eff_weight);
    otree->SetBranchAddress("f_lept1_pt", &f_lept1_pt);
    otree->SetBranchAddress("f_lept1_eta", &f_lept1_eta);
    otree->SetBranchAddress("f_lept1_phi", &f_lept1_phi);
    otree->SetBranchAddress("f_lept1_charge", &f_lept1_charge);
    otree->SetBranchAddress("f_lept1_pfx", &f_lept1_pfx);
    otree->SetBranchAddress("f_lept1_sip", &f_lept1_sip);
    otree->SetBranchAddress("f_lept1_pdgid", &f_lept1_pdgid);
    otree->SetBranchAddress("f_lept2_pt", &f_lept2_pt);
    otree->SetBranchAddress("f_lept2_eta", &f_lept2_eta);
    otree->SetBranchAddress("f_lept2_phi", &f_lept2_phi);
    otree->SetBranchAddress("f_lept2_charge", &f_lept2_charge);
    otree->SetBranchAddress("f_lept2_pfx", &f_lept2_pfx);
    otree->SetBranchAddress("f_lept2_sip", &f_lept2_sip);
    otree->SetBranchAddress("f_lept2_pdgid", &f_lept2_pdgid);
    otree->SetBranchAddress("f_lept3_pt", &f_lept3_pt);
    otree->SetBranchAddress("f_lept3_eta", &f_lept3_eta);
    otree->SetBranchAddress("f_lept3_phi", &f_lept3_phi);
    otree->SetBranchAddress("f_lept3_charge", &f_lept3_charge);
    otree->SetBranchAddress("f_lept3_pfx", &f_lept3_pfx);
    otree->SetBranchAddress("f_lept3_sip", &f_lept3_sip);
    otree->SetBranchAddress("f_lept3_pdgid", &f_lept3_pdgid);
    otree->SetBranchAddress("f_lept4_pt", &f_lept4_pt);
    otree->SetBranchAddress("f_lept4_eta", &f_lept4_eta);
    otree->SetBranchAddress("f_lept4_phi", &f_lept4_phi);
    otree->SetBranchAddress("f_lept4_charge", &f_lept4_charge);
    otree->SetBranchAddress("f_lept4_pfx", &f_lept4_pfx);
    otree->SetBranchAddress("f_lept4_sip", &f_lept4_sip);
    otree->SetBranchAddress("f_lept4_pdgid", &f_lept4_pdgid);
    otree->SetBranchAddress("f_iso_max", &f_iso_max);
    otree->SetBranchAddress("f_sip_max", &f_sip_max);
    otree->SetBranchAddress("f_Z1mass", &f_Z1mass);
    otree->SetBranchAddress("f_Z2mass", &f_Z2mass);
    otree->SetBranchAddress("f_angle_costhetastar", &f_angle_costhetastar);
    otree->SetBranchAddress("f_angle_costheta1", &f_angle_costheta1);
    otree->SetBranchAddress("f_angle_costheta2", &f_angle_costheta2);
    otree->SetBranchAddress("f_angle_phi", &f_angle_phi);
    otree->SetBranchAddress("f_angle_phistar1", &f_angle_phistar1);
    otree->SetBranchAddress("f_pt4l", &f_pt4l);
    otree->SetBranchAddress("f_eta4l", &f_eta4l);
    otree->SetBranchAddress("f_mass4l", &f_mass4l);
    otree->SetBranchAddress("f_mass4lErr", &f_mass4lErr);
    otree->SetBranchAddress("f_njets_pass", &f_njets_pass);
    otree->SetBranchAddress("f_deltajj", &f_deltajj);
    otree->SetBranchAddress("f_massjj", &f_massjj);
    otree->SetBranchAddress("f_D_jet", &f_D_jet);
    otree->SetBranchAddress("f_jet1_pt", &f_jet1_pt);
    otree->SetBranchAddress("f_jet1_eta", &f_jet1_eta);
    otree->SetBranchAddress("f_jet1_phi", &f_jet1_phi);
    otree->SetBranchAddress("f_jet1_e", &f_jet1_e);
    otree->SetBranchAddress("f_jet2_pt", &f_jet2_pt);
    otree->SetBranchAddress("f_jet2_eta", &f_jet2_eta);
    otree->SetBranchAddress("f_jet2_phi", &f_jet2_phi);
    otree->SetBranchAddress("f_jet2_e", &f_jet2_e);
    otree->SetBranchAddress("f_jets_dnn_pt", &f_jets_dnn_pt);
    otree->SetBranchAddress("f_jets_dnn_eta", &f_jets_dnn_eta);
    otree->SetBranchAddress("f_jets_dnn_phi", &f_jets_dnn_phi);
    otree->SetBranchAddress("f_jets_dnn_e", &f_jets_dnn_e);
    otree->SetBranchAddress("f_D_bkg_kin", &f_D_bkg_kin);
    otree->SetBranchAddress("f_D_bkg", &f_D_bkg);
    otree->SetBranchAddress("f_D_gg", &f_D_gg);
    otree->SetBranchAddress("f_D_g4", &f_D_g4);
    otree->SetBranchAddress("f_Djet_VAJHU", &f_Djet_VAJHU);
    otree->SetBranchAddress("f_pfmet", &f_pfmet);
    otree->SetBranchAddress("f_genmet", &f_genmet);
    otree->SetBranchAddress("f_mT", &f_mT);
    otree->SetBranchAddress("f_dphi", &f_dphi);
    otree->SetBranchAddress("f_category", &f_category);
    otree->SetBranchAddress("f_Ngood", &f_Ngood);
    otree->SetBranchAddress("f_Nbjets", &f_Nbjets);
    Int_t Nentries = otree->GetEntries();

    for(int ievent=0; ievent<Nentries; ievent++){
      otree->GetEntry(ievent);
      
      //if(Nentries > 10 && ievent % (Nentries/10) == 0)
	//cout<<"iEvent -------------- "<<ievent<<endl;
	
      if( f_mass4l >= 118 && f_mass4l <= 130 && (((f_njets_pass == 2 || f_njets_pass == 3) && f_Nbjets <= 1) || (f_njets_pass > 3 && f_Nbjets == 0)) ){
	outfile << "{";
        outfile << f_run;
        outfile << ", ";
	outfile << f_lumi;
	outfile << ", ";
        outfile << f_event;
        outfile << ", ";
        outfile << f_jets_dnn_eta[0];
        outfile << ", ";
        outfile << f_jets_dnn_phi[0];
        outfile << ", ";
	outfile << f_jets_dnn_eta[1];
	outfile << ", ";
	outfile << f_jets_dnn_phi[1];
	outfile << ", ";
	outfile << f_jets_dnn_eta[2];
	outfile << ", ";
	outfile << f_jets_dnn_phi[2];
	outfile << "},\n";

	  /*
	histos[itag][0]->Fill(f_lept1_pt, f_weight);
	histos[itag][1]->Fill(f_lept1_eta, f_weight);
	histos[itag][2]->Fill(f_lept1_phi, f_weight);
	histos[itag][3]->Fill(f_lept1_pfx, f_weight);
	histos[itag][4]->Fill(f_lept1_sip, f_weight);
	histos[itag][5]->Fill(f_lept2_pt, f_weight);
	histos[itag][6]->Fill(f_lept2_eta, f_weight);
	histos[itag][7]->Fill(f_lept2_phi, f_weight);
	histos[itag][8]->Fill(f_lept2_pfx, f_weight);
	histos[itag][9]->Fill(f_lept2_sip, f_weight);
	histos[itag][10]->Fill(f_lept3_pt, f_weight);
	histos[itag][11]->Fill(f_lept3_eta, f_weight);
	histos[itag][12]->Fill(f_lept3_phi, f_weight);
	histos[itag][13]->Fill(f_lept3_pfx, f_weight);
	histos[itag][14]->Fill(f_lept3_sip, f_weight);
	histos[itag][15]->Fill(f_lept4_pt, f_weight);
	histos[itag][16]->Fill(f_lept4_eta, f_weight);
	histos[itag][17]->Fill(f_lept4_phi, f_weight);
	histos[itag][18]->Fill(f_lept4_pfx, f_weight);
	histos[itag][19]->Fill(f_lept4_sip, f_weight);
	histos[itag][20]->Fill(f_Z1mass, f_weight);
	histos[itag][21]->Fill(f_Z2mass, f_weight);
	histos[itag][22]->Fill(f_angle_costhetastar, f_weight);
	histos[itag][23]->Fill(f_angle_costheta1, f_weight);
	histos[itag][24]->Fill(f_angle_costheta2, f_weight);
	histos[itag][25]->Fill(f_angle_phi, f_weight);
	histos[itag][26]->Fill(f_angle_phistar1, f_weight);
	histos[itag][27]->Fill(f_pt4l, f_weight);
	histos[itag][28]->Fill(f_eta4l, f_weight);
	histos[itag][29]->Fill(f_mass4l, f_weight);
	histos[itag][30]->Fill(f_njets_pass, f_weight);
	histos[itag][31]->Fill(f_deltajj, f_weight);
	histos[itag][32]->Fill(f_massjj, f_weight);
	histos[itag][33]->Fill(f_D_jet, f_weight);
	histos[itag][34]->Fill(f_jets_dnn_pt[0], f_weight);
	histos[itag][35]->Fill(f_jets_dnn_eta[0], f_weight);
	histos[itag][36]->Fill(f_jets_dnn_phi[0], f_weight);
	histos[itag][37]->Fill(f_jets_dnn_e[0], f_weight);
	histos[itag][38]->Fill(f_jets_dnn_pt[1], f_weight);
	histos[itag][39]->Fill(f_jets_dnn_eta[1], f_weight);
	histos[itag][40]->Fill(f_jets_dnn_phi[1], f_weight);
	histos[itag][41]->Fill(f_jets_dnn_e[1], f_weight);
	histos[itag][42]->Fill(f_jets_dnn_pt[2], f_weight);
	histos[itag][43]->Fill(f_jets_dnn_eta[2], f_weight);
	histos[itag][44]->Fill(f_jets_dnn_phi[2], f_weight);
	histos[itag][45]->Fill(f_jets_dnn_e[2], f_weight);
	histos[itag][46]->Fill(f_D_bkg_kin, f_weight);
	histos[itag][47]->Fill(f_D_bkg, f_weight);
	histos[itag][48]->Fill(f_D_gg, f_weight);
	histos[itag][49]->Fill(f_D_g4, f_weight);
	histos[itag][50]->Fill(f_Djet_VAJHU, f_weight);
	*/
      }//End of "if" for VBF selection
    }            
    ofile->Close();  
  }//end samples loop
  }//end tags loop
  
  gROOT->SetBatch();
  
  ////Plot to do
  //TPaveText *cms_tag = new TPaveText(.34,.99,.9,1,"NDC");
  //cms_tag->AddText("CMS #bf{Preliminary  #sqrt{s} = 13 TeV, L = 35.6fb^{-1}}");
  //cms_tag->SetFillStyle(0);
  //cms_tag->SetBorderSize(0);
  //cms_tag->SetTextSize(0.14);
  /*
  TLegend *leg = new TLegend(0.72,0.63,0.88,0.97);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetEntrySeparation(0.1);

  TCanvas *cv = new TCanvas("cv","",10,10,700,700);
  TPad *pad1 = new TPad("pad1","",0.05,0.74,0.95,0.98);
  TPad *pad2 = new TPad("pad2","",0.05,0.02,0.95,0.75);
  pad1->SetFillColor(0);
  pad1->SetTopMargin(0.98);
  pad2->SetFillColor(0);
  pad2->SetTopMargin(0);
  pad1->Draw();
  pad2->Draw();

  THStack *shistos[nhists];
  TGraphAsymmErrors *data[nhists];
  TGraphErrors *ratio[nhists];
  TF1 *f1 = new TF1("f1","[0]",-10000,10000);
  unsigned int color[] = {2,3,4,5,6,7,8,9,28,30,41,46};
  for(unsigned int ihist=0; ihist<nhists; ++ihist){
      shistos[ihist] = new THStack();
      data[ihist] = new TGraphAsymmErrors();
      ratio[ihist] = new TGraphErrors();
      for(unsigned int itag=0; itag<ntags; ++itag){
	//stack MCs and Data
	if(tags[itag] == "Data"){
	  for(unsigned int ib=0; ib<histos[itag][ihist]->GetNbinsX(); ++ib){
	  float x = histos[itag][ihist]->GetBinCenter(ib+1);
	  float y = histos[itag][ihist]->GetBinContent(ib+1);
	  data[ihist]->SetPoint(ib,x,y);
	  
	  float emy = 0, epy = 0;
	  if(y != 0){
	    emy = -0.5 + sqrt(y+0.25);
	    epy = +0.5 + sqrt(y+0.25);
	  }
	  data[ihist]->SetPointError(ib,0,0,emy,epy);
	}
	if(ihist == 0)
	  leg->AddEntry(data[ihist],tags[itag],"p");
      }
      else{
	histos[itag][ihist]->SetLineColor(color[itag]);
	histos[itag][ihist]->SetFillColor(color[itag]);
	shistos[ihist]->Add( histos[itag][ihist] );
	if(ihist == 0)
	  leg->AddEntry(histos[itag][ihist],tags[itag],"f");
      }
    }
    
    //Computes data-mc ratio
    unsigned int ipoint = 0;
    unsigned int Nbins = ((TH1*)(shistos[ihist]->GetStack()->Last()))->GetNbinsX();
    for(unsigned int ib=0; ib<Nbins; ++ib){
      float stack_bin_content = ((TH1*)(shistos[ihist]->GetStack()->Last()))->GetBinContent(ib+1);
      if(stack_bin_content > 0){
	double dx, dy;
	data[ihist]->GetPoint(ib,dx,dy);
	if(dy == 0) continue;
	float rat = dy/stack_bin_content;
	float rat_err = sqrt((1./(stack_bin_content*stack_bin_content))*dy + (dy*dy/pow(stack_bin_content,4))*stack_bin_content);
	ratio[ihist]->SetPoint(ipoint,dx,rat);
	ratio[ihist]->SetPointError(ipoint,0,rat_err);
        ++ipoint;
	//cout<<"dx: "<<dx<<", rat: "<<rat<<endl;
      }
    }
    
    
    cout<<"Fiting "<<histos_name[ihist]<<endl;
    ratio[ihist]->Fit(f1,"","e",data[ihist]->GetXaxis()->GetXmin(),data[ihist]->GetXaxis()->GetXmax());
    ratio[ihist]->GetYaxis()->SetRangeUser(0,2);
    ratio[ihist]->GetXaxis()->SetLimits(data[ihist]->GetXaxis()->GetXmin(),data[ihist]->GetXaxis()->GetXmax());
    ratio[ihist]->SetMarkerColor(kBlack);
    ratio[ihist]->GetXaxis()->SetLabelSize(0);
    ratio[ihist]->GetYaxis()->SetLabelSize(0.13);
    ratio[ihist]->GetYaxis()->SetTitle("Data/MC");
    ratio[ihist]->GetYaxis()->SetTitleSize(0.13);
    ratio[ihist]->GetYaxis()->SetTitleOffset(0.3);
    
    cout<<"Drawing histograms..."<<endl;
    pad1->cd();
    ratio[ihist]->Draw("aep");
    ratio[ihist]->GetYaxis()->SetRangeUser(0,2);
    TPaveText *fit_res = new TPaveText(.55,.85,.9,.88,"NDC");
    fit_res->AddText(Form("%.3f #pm %.3f",f1->GetParameter(0),f1->GetParError(0)));
    fit_res->SetFillStyle(0);
    fit_res->SetBorderSize(0);
    fit_res->SetTextSize(0.14);
    fit_res->SetTextColor(kRed);
    fit_res->Draw();    
    
    pad2->cd();
    data[ihist]->SetMarkerColor(kBlack);
    data[ihist]->SetMarkerStyle(20);
    data[ihist]->GetXaxis()->SetTitle(histos_name[ihist]);
    TH1D *last = shistos[ihist]->GetStack()->Last();
    float binning = last->GetBinWidth(1);
    data[ihist]->GetYaxis()->SetTitle(Form("Events/%.3f",binning));
    data[ihist]->GetYaxis()->SetTitleOffset(1);
    data[ihist]->Draw("aep");
    gPad->SetLogy(true);
    //data[ihist]->GetYaxis()->SetRangeUser(0.1,1000);
    data[ihist]->GetYaxis()->SetRangeUser(0.001,100);
    shistos[ihist]->Draw("hist,same");
    data[ihist]->Draw("ep,same");
    leg->Draw();
    cv->Update();
    cv->Print(Form("Temp/histo_%i_log.png",ihist));
    
  }      
  */
      
  return;
}//End plot function


void plotObjectsProperties(void){
  
  std::vector<TString> Data;
  std::vector<TString> VBF;
  std::vector<TString> HJJ;
  std::vector<TString> ggH;
  std::vector<TString> qqZZ;
  std::vector<TString> ggZZ;
  std::vector<TString> ttH;
  std::vector<TString> VH;
  std::vector<TString> ttX;
  std::vector<TString> WX;
  std::vector<TString> ZZZ;
  
  
  //4mu
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016B-23Sep2016-v3_miniaod_1.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016B-23Sep2016-v3_miniaod_2.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016C-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016D-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016E-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016F-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016G-23Sep2016-v1_miniaod_1.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016G-23Sep2016-v1_miniaod_2.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016H-PromptReco-v2_miniaod.root");
  Data.push_back("histos4mu_25ns/output_DoubleMuon_Run2016H-PromptReco-v3_miniaod.root");
  VBF.push_back("histos4mu_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  HJJ.push_back("histos4mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  ggH.push_back("histos4mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  //qqZZ.push_back("histos4mu_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos4mu_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");  
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  ttH.push_back("histos4mu_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  VH.push_back("histos4mu_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos4mu_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos4mu_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  ttX.push_back("histos4mu_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  ttX.push_back("histos4mu_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  ttX.push_back("histos4mu_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos4mu_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  WX.push_back("histos4mu_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  WX.push_back("histos4mu_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos4mu_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  WX.push_back("histos4mu_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  ZZZ.push_back("histos4mu_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
 
  //4e
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016B-23Sep2016-v3_miniaod_1.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016B-23Sep2016-v3_miniaod_2.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016C-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016D-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016E-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016F-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016G-23Sep2016-v1_miniaod.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016H-PromptReco-v2_miniaod_1.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016H-PromptReco-v2_miniaod_2.root");
  Data.push_back("histos4e_25ns/output_DoubleEG_Run2016H-PromptReco-v3_miniaod.root");
  VBF.push_back("histos4e_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  HJJ.push_back("histos4e_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  ggH.push_back("histos4e_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  //qqZZ.push_back("histos4e_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos4e_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  ttH.push_back("histos4e_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");  
  VH.push_back("histos4e_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos4e_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos4e_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");  
  ttX.push_back("histos4e_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  ttX.push_back("histos4e_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  ttX.push_back("histos4e_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");  
  WX.push_back("histos4e_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  WX.push_back("histos4e_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  WX.push_back("histos4e_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos4e_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  WX.push_back("histos4e_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  ZZZ.push_back("histos4e_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");

  
  //2e2mu
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016B-23Sep2016-v3_miniaod_1.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016B-23Sep2016-v3_miniaod_2.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016C-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016D-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016E-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016F-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016G-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016H-PromptReco-v2_miniaod_1.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016H-PromptReco-v2_miniaod_2.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleEG_Run2016H-PromptReco-v3_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016B-23Sep2016-v3_miniaod_1.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016B-23Sep2016-v3_miniaod_2.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016C-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016D-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016E-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016F-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016G-23Sep2016-v1_miniaod_1.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016G-23Sep2016-v1_miniaod_2.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016H-PromptReco-v2_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_DoubleMuon_Run2016H-PromptReco-v3_miniaod.root");      
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016B-23Sep2016-v3_miniaod_1.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016B-23Sep2016-v3_miniaod_2.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016C-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016D-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016E-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016F-23Sep2016-v1_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016G-23Sep2016-v1_miniaod_1.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016G-23Sep2016-v1_miniaod_2.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016H-PromptReco-v2_miniaod.root");
  Data.push_back("histos2e2mu_25ns/output_MuonEG_Run2016H-PromptReco-v3_miniaod.root");  
  VBF.push_back("histos2e2mu_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");   
  HJJ.push_back("histos2e2mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");  
  ggH.push_back("histos2e2mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");  
  //qqZZ.push_back("histos2e2mu_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos2e2mu_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  ttH.push_back("histos2e2mu_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  VH.push_back("histos2e2mu_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos2e2mu_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");  
  VH.push_back("histos2e2mu_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  ttX.push_back("histos2e2mu_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  ttX.push_back("histos2e2mu_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  ttX.push_back("histos2e2mu_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  WX.push_back("histos2e2mu_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  ZZZ.push_back("histos2e2mu_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  
  
  std::vector<std::vector<TString>> samples;
  std::vector<TString> tags;

  //samples.push_back(Data);	tags.push_back("Data");
  //samples.push_back(VBF);	tags.push_back("VBF");
  //samples.push_back(ggH);	tags.push_back("ggH");
  //samples.push_back(HJJ);	tags.push_back("HJJ");
  samples.push_back(qqZZ);	tags.push_back("qqZZ");
  //samples.push_back(ggZZ);	tags.push_back("ggZZ");
  //samples.push_back(ttH);	tags.push_back("ttH");
  //samples.push_back(VH);	tags.push_back("VH");
  //samples.push_back(ttX);	tags.push_back("ttX");
  //samples.push_back(WX);	tags.push_back("WX");
  //samples.push_back(ZZZ);	tags.push_back("ZZZ");

  loopSamples(tags, samples);
     
  //Ends program
  return;
}
