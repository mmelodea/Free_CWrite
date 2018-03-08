#include <iostream>
#include <fstream>
#include <vector>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TF1.h>
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

//#include "/home/micah/cernbox/MonoHiggsHZZ4L/model963_7.cc"

const unsigned int nhists = 60;
TString samples_path = "OriginalSamples/";


void loopSamples(std::vector<TString> tags, std::vector< std::vector<TString> > samples){
  TH1::SetDefaultSumw2();
  
  //ofstream outfile;
  //outfile.open("/home/micah/temp/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_jets.h");
  //outfile.open("/home/micah/temp/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8_jets.h");
  //outfile.open("/home/micah/temp/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_jets.h");
  //outfile.open("/home/micah/temp/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_jets.h");
  //outfile.open("/home/micah/temp/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8_jets.h");
  //outfile.open("/home/micah/temp/ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8_jets.h");
  //outfile.open("/home/micah/temp/GluGluToContinToZZTo4L_13TeV_MCFM701_pythia8_jets.h");
  //outfile.open("/home/micah/temp/ZZTo4L_13TeV_powheg_pythia8_jets.h");
  
  
  Float_t f_weight, f_int_weight, f_pu_weight, f_eff_weight, f_lept1_pt, f_lept1_eta, f_lept1_phi, f_lept1_charge, f_lept1_pfx, f_lept1_sip, f_lept1_mvaid, f_lept2_pt, f_lept2_eta, f_lept2_phi, f_lept2_charge, f_lept2_pfx, f_lept2_sip, f_lept2_mvaid, f_lept3_pt, f_lept3_eta, f_lept3_phi, f_lept3_charge, f_lept3_pfx, f_lept3_sip, f_lept3_mvaid, f_lept4_pt, f_lept4_eta, f_lept4_phi, f_lept4_charge, f_lept4_pfx, f_lept4_sip, f_lept4_mvaid, f_iso_max, f_sip_max, f_Z1mass, f_Z2mass, f_angle_costhetastar, f_angle_costheta1, f_angle_costheta2, f_angle_phi, f_angle_phistar1, f_eta4l, f_pt4l, f_mass4l, f_mass4lErr, f_njets_pass, f_deltajj, f_massjj, f_D_jet;
  Float_t f_jet1_highpt_pt, f_jet1_highpt_eta, f_jet1_highpt_phi, f_jet1_highpt_e, f_jet2_highpt_pt, f_jet2_highpt_eta, f_jet2_highpt_phi, f_jet2_highpt_e, f_jet3_highpt_pt, f_jet3_highpt_eta, f_jet3_highpt_phi, f_jet3_highpt_e;
  Float_t f_jet4_highpt_pt, f_jet4_highpt_eta, f_jet4_highpt_phi, f_jet4_highpt_e, f_jet5_highpt_pt, f_jet5_highpt_eta, f_jet5_highpt_phi, f_jet5_highpt_e, f_jet6_highpt_pt, f_jet6_highpt_eta, f_jet6_highpt_phi, f_jet6_highpt_e;
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
  histos_name.push_back("lept4_pt");
  histos_name.push_back("lept4_eta");
  histos_name.push_back("lept4_phi");
  histos_name.push_back("lept4_pfx");
  histos_name.push_back("lept4_sip");
  histos_name.push_back("Z1mass");
  histos_name.push_back("Z2mass");
  histos_name.push_back("costhetastar");
  histos_name.push_back("costheta1");
  histos_name.push_back("costheta2");
  histos_name.push_back("thetastar");
  histos_name.push_back("theta1");
  histos_name.push_back("theta2");
  histos_name.push_back("phi");
  histos_name.push_back("phistar1");
  histos_name.push_back("pt4l");
  histos_name.push_back("eta4l");
  histos_name.push_back("mass4l");
  histos_name.push_back("mass4l_zoom");
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
  histos_name.push_back("pf_met");
  histos_name.push_back("nbjets");
  histos_name.push_back("dphi");
  histos_name.push_back("mT");  
  histos_name.push_back("Network_Discriminant");
  
    
  float fsig_err = 0, fbkg_err = 0;
  float fnsig = 0, fnbkg = 0, fndata = 0;
  const unsigned int ntags = tags.size();
  TH1D *histos[ntags][nhists];
  for(unsigned int itag=0; itag<ntags; ++itag){
    float events = 0, yields = 0;
    float fjet3 = 0, fjet4 = 0, fjet5 = 0, fjet6 = 0;
    
    histos[itag][0] = new TH1D("","",100,0,250);
    histos[itag][1] = new TH1D("","",50,-2.5,2.5);
    histos[itag][2] = new TH1D("","",64,-3.2,3.2);
    histos[itag][3] = new TH1D("","",100,0,0.4);
    histos[itag][4] = new TH1D("","",100,-5,5);
    histos[itag][5] = new TH1D("","",100,0,250);
    histos[itag][6] = new TH1D("","",50,-2.5,2.5);
    histos[itag][7] = new TH1D("","",64,-3.2,3.2);
    histos[itag][8] = new TH1D("","",100,0,0.4);
    histos[itag][9] = new TH1D("","",100,-5,5);
    histos[itag][10] = new TH1D("","",100,0,250);
    histos[itag][11] = new TH1D("","",50,-2.5,2.5);
    histos[itag][12] = new TH1D("","",64,-3.2,3.2);
    histos[itag][13] = new TH1D("","",100,0,0.4);
    histos[itag][14] = new TH1D("","",100,-5,5);
    histos[itag][15] = new TH1D("","",100,0,250);
    histos[itag][16] = new TH1D("","",50,-2.5,2.5);
    histos[itag][17] = new TH1D("","",64,-3.2,3.2);
    histos[itag][18] = new TH1D("","",100,0,0.4);
    histos[itag][19] = new TH1D("","",100,-5,5);
    histos[itag][20] = new TH1D("","",18,12,120);//
    histos[itag][21] = new TH1D("","",18,12,120);//
    histos[itag][22] = new TH1D("","",50,-1.05,1.05);
    histos[itag][23] = new TH1D("","",50,-1.05,1.05);
    histos[itag][24] = new TH1D("","",50,-1.05,1.05);
    histos[itag][25] = new TH1D("","",50,-0.05,3.16);
    histos[itag][26] = new TH1D("","",50,-0.05,3.16);    
    histos[itag][27] = new TH1D("","",50,-0.05,3.16);
    histos[itag][28] = new TH1D("","",50,-3.16,3.16);
    histos[itag][29] = new TH1D("","",50,-3.16,3.16);
    histos[itag][30] = new TH1D("","",250,0,500);
    histos[itag][31] = new TH1D("","",100,-5,5);
    histos[itag][32] = new TH1D("","",110,80,300);
    histos[itag][33] = new TH1D("","",7,111,139);
    histos[itag][34] = new TH1D("","",10,0,10);
    histos[itag][35] = new TH1D("","",8,0,12);//
    histos[itag][36] = new TH1D("","",5,0,4000);
    histos[itag][37] = new TH1D("","",5,0,2);
    histos[itag][38] = new TH1D("","",175,0,350);
    histos[itag][39] = new TH1D("","",50,-5,5);
    histos[itag][40] = new TH1D("","",64,-3.2,3.2);
    histos[itag][41] = new TH1D("","",100,0,500);
    histos[itag][42] = new TH1D("","",175,0,350);
    histos[itag][43] = new TH1D("","",50,-5,5);
    histos[itag][44] = new TH1D("","",64,-3.2,3.2);
    histos[itag][45] = new TH1D("","",100,0,500);
    histos[itag][46] = new TH1D("","",175,0,350);
    histos[itag][47] = new TH1D("","",50,-5,5);
    histos[itag][48] = new TH1D("","",64,-3.2,3.2);
    histos[itag][49] = new TH1D("","",100,0,500);
    histos[itag][50] = new TH1D("","",5,-0.05,1.05);
    histos[itag][51] = new TH1D("","",5,-0.05,1.05);
    histos[itag][52] = new TH1D("","",5,-0.05,1.05);
    histos[itag][53] = new TH1D("","",5,-0.05,1.05);
    histos[itag][54] = new TH1D("","",5,-0.05,1.05);
    histos[itag][55] = new TH1D("","",30,0,300);//
    histos[itag][56] = new TH1D("","",10,0,10);
    histos[itag][57] = new TH1D("","",100,-4,4);
    histos[itag][58] = new TH1D("","",50,0,500);
    histos[itag][59] = new TH1D("","",5,-0.05,1.05);
      
    unsigned int nsamples = samples[itag].size();
    for(unsigned int is=0; is<nsamples; ++is){
      TString ifile_name_root = samples_path+samples[itag][is];
      if(gSystem->AccessPathName(ifile_name_root)){
	cout<<"File "<<ifile_name_root<<" doesn't exist!"<<endl;
	continue;
      }
      //else cout<<"Loading file "<<ifile_name_root<<endl;
  
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
    otree->SetBranchAddress("f_jet1_highpt_pt", &f_jet1_highpt_pt);
    otree->SetBranchAddress("f_jet1_highpt_eta", &f_jet1_highpt_eta);
    otree->SetBranchAddress("f_jet1_highpt_phi", &f_jet1_highpt_phi);
    otree->SetBranchAddress("f_jet1_highpt_e", &f_jet1_highpt_e);
    otree->SetBranchAddress("f_jet2_highpt_pt", &f_jet2_highpt_pt);
    otree->SetBranchAddress("f_jet2_highpt_eta", &f_jet2_highpt_eta);
    otree->SetBranchAddress("f_jet2_highpt_phi", &f_jet2_highpt_phi);
    otree->SetBranchAddress("f_jet2_highpt_e", &f_jet2_highpt_e);
    otree->SetBranchAddress("f_jet3_highpt_pt", &f_jet3_highpt_pt);
    otree->SetBranchAddress("f_jet3_highpt_eta", &f_jet3_highpt_eta);
    otree->SetBranchAddress("f_jet3_highpt_phi", &f_jet3_highpt_phi);
    otree->SetBranchAddress("f_jet3_highpt_e", &f_jet3_highpt_e);
    otree->SetBranchAddress("f_jet4_highpt_pt", &f_jet4_highpt_pt);
    otree->SetBranchAddress("f_jet4_highpt_eta", &f_jet4_highpt_eta);
    otree->SetBranchAddress("f_jet4_highpt_phi", &f_jet4_highpt_phi);
    otree->SetBranchAddress("f_jet4_highpt_e", &f_jet4_highpt_e);
    otree->SetBranchAddress("f_jet5_highpt_pt", &f_jet5_highpt_pt);
    otree->SetBranchAddress("f_jet5_highpt_eta", &f_jet5_highpt_eta);
    otree->SetBranchAddress("f_jet5_highpt_phi", &f_jet5_highpt_phi);
    otree->SetBranchAddress("f_jet5_highpt_e", &f_jet5_highpt_e);
    otree->SetBranchAddress("f_jet6_highpt_pt", &f_jet6_highpt_pt);
    otree->SetBranchAddress("f_jet6_highpt_eta", &f_jet6_highpt_eta);
    otree->SetBranchAddress("f_jet6_highpt_phi", &f_jet6_highpt_phi);
    otree->SetBranchAddress("f_jet6_highpt_e", &f_jet6_highpt_e);
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

    int pass_events = 0;
    for(int ievent=0; ievent<Nentries; ievent++){
      otree->GetEntry(ievent);
      
      //if(Nentries > 10 && ievent % (Nentries/10) == 0)
	//cout<<"iEvent -------------- "<<ievent<<endl;
	
      if( (((f_njets_pass == 2 || f_njets_pass == 3) && f_Nbjets <= 1) || (f_njets_pass > 3 && f_Nbjets == 0)) && f_mass4l >= 118 && f_mass4l <= 130 ){
	//cout<<"Pass events: "<<pass_events<<endl;
	/*
	outfile << "{";
        outfile << f_run;
        outfile << ", ";
	outfile << f_lumi;
	outfile << ", ";
        outfile << f_event;
        outfile << ", ";
	outfile << f_jet1_highpt_eta;
	outfile << ", ";
        outfile << f_jet1_highpt_phi;
        outfile << ", ";
        outfile << f_jet2_highpt_eta;
        outfile << ", ";
	outfile << f_jet2_highpt_phi;
        outfile << ", ";
	outfile << f_jet3_highpt_eta;
	outfile << ", ";
	outfile << f_jet3_highpt_phi;
	outfile << "},\n";
        */
	
	std::vector<double> inputs;
	
	inputs.push_back( f_lept1_pt );
	inputs.push_back( f_lept1_eta );
	inputs.push_back( f_lept1_phi );
	inputs.push_back( f_lept2_pt );
	inputs.push_back( f_lept2_eta );
	inputs.push_back( f_lept2_phi );
	inputs.push_back( f_lept3_pt );
	inputs.push_back( f_lept3_eta );
	inputs.push_back( f_lept3_phi );
	inputs.push_back( f_lept4_pt );
	inputs.push_back( f_lept4_eta );
	inputs.push_back( f_lept4_phi );
	
	inputs.push_back( f_jet1_highpt_pt );
	inputs.push_back( f_jet1_highpt_eta );
	inputs.push_back( f_jet1_highpt_phi );
	//inputs.push_back( f_jet1_highpt_e*TMath::CosH(f_jet1_highpt_eta) );
	inputs.push_back( f_jet2_highpt_pt );
	inputs.push_back( f_jet2_highpt_eta );
	inputs.push_back( f_jet2_highpt_phi );
	//inputs.push_back( f_jet2_highpt_e*TMath::CosH(f_jet2_highpt_eta) );	
	
	if( f_jet3_highpt_pt != -999 ){
	  inputs.push_back( f_jet3_highpt_pt );
	  inputs.push_back( f_jet3_highpt_eta );
	  inputs.push_back( f_jet3_highpt_phi );
	  //inputs.push_back( f_jet3_highpt_e*TMath::CosH(f_jet3_highpt_eta) );	
	}else{
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	  //inputs.push_back( 0 );
	}
        /*
	if( f_jet4_highpt_pt != -999 ){
	  inputs.push_back( f_jet4_highpt_pt );
	  inputs.push_back( f_jet4_highpt_eta );
	  inputs.push_back( f_jet4_highpt_phi );
	}else{
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	}
	
	if( f_jet5_highpt_pt != -999 ){
	  inputs.push_back( f_jet5_highpt_pt );
	  inputs.push_back( f_jet5_highpt_eta );
	  inputs.push_back( f_jet5_highpt_phi );
	}else{
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	}
	
	if( f_jet6_highpt_pt != -999 ){
	  inputs.push_back( f_jet6_highpt_pt );
	  inputs.push_back( f_jet6_highpt_eta );
	  inputs.push_back( f_jet6_highpt_phi );
	}else{
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	  inputs.push_back( 0 );
	}
	
	inputs.push_back( f_pfmet );
	inputs.push_back( f_njets_pass );
	inputs.push_back( f_Nbjets );
	*/
	float dnn_discriminant = 0.5;//when not using nets
	//dnn_discriminant = dnn_303(inputs);
	//if(dnn_discriminant <= 0.660) continue;	
	//if(f_Djet_VAJHU <= 0.480) continue;
	//if(f_D_jet <= 0.600) continue;	

        ++events;
        yields += f_weight;
        ++pass_events;
	if(tags[itag] == "Data"){
	  fndata += f_weight;
	}else if(tags[itag] == "VBF"){
          fnsig += f_weight;
	  fsig_err += f_weight*f_weight;
        }else{
          fnbkg += f_weight;
 	  fbkg_err += f_weight*f_weight;
	}
	
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
	histos[itag][25]->Fill(TMath::ACos(f_angle_costhetastar), f_weight);
	histos[itag][26]->Fill(TMath::ACos(f_angle_costheta1), f_weight);
	histos[itag][27]->Fill(TMath::ACos(f_angle_costheta2), f_weight);
	histos[itag][28]->Fill(f_angle_phi, f_weight);
	histos[itag][29]->Fill(f_angle_phistar1, f_weight);
	histos[itag][30]->Fill(f_pt4l, f_weight);
	histos[itag][31]->Fill(f_eta4l, f_weight);
	histos[itag][32]->Fill(f_mass4l, f_weight);
	if(f_mass4l >= 118 && f_mass4l <=130){
	  histos[itag][33]->Fill(f_mass4l, f_weight);
	}
	histos[itag][34]->Fill(f_njets_pass, f_weight);
	histos[itag][35]->Fill(f_deltajj, f_weight);
	histos[itag][36]->Fill(f_massjj, f_weight);
	histos[itag][37]->Fill(f_D_jet, f_weight);
	histos[itag][38]->Fill(f_jet1_highpt_pt, f_weight);
	histos[itag][39]->Fill(f_jet1_highpt_eta, f_weight);
	histos[itag][40]->Fill(f_jet1_highpt_phi, f_weight);
	histos[itag][41]->Fill(f_jet1_highpt_e, f_weight);
	histos[itag][42]->Fill(f_jet2_highpt_pt, f_weight);
	histos[itag][43]->Fill(f_jet2_highpt_eta, f_weight);
	histos[itag][44]->Fill(f_jet2_highpt_phi, f_weight);
	histos[itag][45]->Fill(f_jet2_highpt_e, f_weight);
	histos[itag][46]->Fill(f_jet3_highpt_pt, f_weight);
	histos[itag][47]->Fill(f_jet3_highpt_eta, f_weight);
	histos[itag][48]->Fill(f_jet3_highpt_phi, f_weight);
	histos[itag][49]->Fill(f_jet3_highpt_e, f_weight);
	histos[itag][50]->Fill(f_D_bkg_kin, f_weight);
	histos[itag][51]->Fill(f_D_bkg, f_weight);
	histos[itag][52]->Fill(f_D_gg, f_weight);
	histos[itag][53]->Fill(f_D_g4, f_weight);
	histos[itag][54]->Fill(f_Djet_VAJHU, f_weight);
	histos[itag][55]->Fill(f_pfmet, f_weight);
	histos[itag][56]->Fill(f_Nbjets, f_weight);
	histos[itag][57]->Fill(f_dphi, f_weight);
	histos[itag][58]->Fill(f_mT, f_weight);
	histos[itag][59]->Fill(dnn_discriminant, f_weight);
	
	
	//See jets contributions
	if(f_jet3_highpt_pt != -999) ++fjet3;
	if(f_jet4_highpt_pt != -999) ++fjet4;
	if(f_jet5_highpt_pt != -999) ++fjet5;
	if(f_jet6_highpt_pt != -999) ++fjet6;
      }//End of "if" for VBF selection
    }            
    ofile->Close();
    //cout<<ifile_name_root<<" = "<<pass_events<<endl;
  }//end samples loop
  
    std::cout<<tags[itag]<<" --- events = "<<events<<" --- yields = "<<yields<<std::endl;
    //std::cout<<"fjet3 = "<<100*(fjet3/events)<<"%"<<std::endl;
    //std::cout<<"fjet4 = "<<100*(fjet4/events)<<"%"<<std::endl;
    //std::cout<<"fjet5 = "<<100*(fjet5/events)<<"%"<<std::endl;
    //std::cout<<"fjet6 = "<<100*(fjet6/events)<<"%"<<std::endl;
  }//end tags loop
  
  
  gROOT->SetBatch();
  ////Plot to do
  TPaveText *cms_tag = new TPaveText(.18,.97,.9,.99,"NDC");
  cms_tag->AddText("CMS #bf{Preliminary  #sqrt{s} = 13 TeV, L = 35.9fb^{-1}}");
  cms_tag->SetFillStyle(0);
  cms_tag->SetBorderSize(0);
  cms_tag->SetTextSize(0.12);
  
  TLegend *leg = new TLegend(0.70,0.63,0.86,0.97);
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
  TGraph *hdata[nhists];
  TGraphAsymmErrors *fdata[nhists];
  TGraphErrors *ratio[nhists];
  TF1 *f1 = new TF1("f1","[0]",-10000,10000);
  unsigned int color[] = {1,2,3,4,5,6,7,12,28,30,36,41,46};
  for(unsigned int ihist=0; ihist<(unsigned int)nhists; ++ihist){
    //if(histos_name[ihist] != "mass4l_zoom") continue;
    
      shistos[ihist] = new THStack();
      hdata[ihist] = new TGraph();
      fdata[ihist] = new TGraphAsymmErrors();
      ratio[ihist] = new TGraphErrors();
      for(unsigned int itag=0; itag<ntags; ++itag){
	if(histos[itag][ihist]->GetEntries() == 0) continue;
	//stack MCs and Data
	if(tags[itag] == "Data"){
	  int ip = 0;
	  for(unsigned int ib=0; ib<(unsigned int)histos[itag][ihist]->GetNbinsX(); ++ib){
	    float x = histos[itag][ihist]->GetBinCenter(ib+1);
	    float y = histos[itag][ihist]->GetBinContent(ib+1);
	  
	    hdata[ihist]->SetPoint(ib,x,y);
	    
	    if(y > 0){
	      float emy = 0, epy = 0;
	      emy = -0.5 + sqrt(y+0.25);
	      epy = +0.5 + sqrt(y+0.25);
	      fdata[ihist]->SetPoint(ip,x,y);
	      fdata[ihist]->SetPointError(ip,0,0,emy,epy);
	      ++ip;
	    }
	  }
	  if(ihist == 0)
	    leg->AddEntry(fdata[ihist],tags[itag],"p");
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
    TH1D *mc_sum = ((TH1D*)(shistos[ihist]->GetStack()->Last()));
    unsigned int Nbins = mc_sum->GetNbinsX();
    for(unsigned int ib=0; ib<Nbins; ++ib){
      float stack_bin_content = mc_sum->GetBinContent(ib+1);
      if(stack_bin_content > 0){
	double dx, dy;
	hdata[ihist]->GetPoint(ib,dx,dy);
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
    ratio[ihist]->Fit(f1,"","",hdata[ihist]->GetXaxis()->GetXmin(),hdata[ihist]->GetXaxis()->GetXmax());
    //ratio[ihist]->GetYaxis()->SetRangeUser(0,2);
    ratio[ihist]->GetXaxis()->SetLimits(hdata[ihist]->GetXaxis()->GetXmin(),hdata[ihist]->GetXaxis()->GetXmax());
    ratio[ihist]->SetMarkerColor(kBlack);
    ratio[ihist]->GetXaxis()->SetLabelSize(0);
    ratio[ihist]->GetYaxis()->SetLabelSize(0.13);
    ratio[ihist]->GetYaxis()->SetTitle("Data/MC");
    ratio[ihist]->GetYaxis()->SetTitleSize(0.13);
    ratio[ihist]->GetYaxis()->SetTitleOffset(0.3);
    
    cout<<"Drawing histograms..."<<endl;
    pad1->cd();
    ratio[ihist]->Draw("aep");
    //ratio[ihist]->GetYaxis()->SetRangeUser(0,2);
    TPaveText *fit_res = new TPaveText(.55,.85,.9,.88,"NDC");
    fit_res->AddText(Form("%.3f #pm %.3f",f1->GetParameter(0),f1->GetParError(0)));
    fit_res->SetFillStyle(0);
    fit_res->SetBorderSize(0);
    fit_res->SetTextSize(0.14);
    fit_res->SetTextColor(kRed);
    fit_res->Draw();    
    cms_tag->Draw();
    
    pad2->cd();
    fdata[ihist]->SetMarkerColor(kBlack);
    fdata[ihist]->SetMarkerStyle(20);
    fdata[ihist]->GetXaxis()->SetTitle(histos_name[ihist]);
    float binning = mc_sum->GetBinWidth(1);
    fdata[ihist]->GetYaxis()->SetTitle(Form("Events/%.3f",binning));
    fdata[ihist]->GetYaxis()->SetTitleOffset(1);
    fdata[ihist]->Draw("aep");
    fdata[ihist]->GetXaxis()->SetLimits(hdata[ihist]->GetXaxis()->GetXmin(),hdata[ihist]->GetXaxis()->GetXmax());
    
    //If wants to plot in log scale: comment first line and uncomment the others
    fdata[ihist]->SetMinimum(0);
    //if(histos_name[ihist] == "mass4l") fdata[ihist]->SetMaximum(6);
    //if(histos_name[ihist] == "Z1mass" || histos_name[ihist] == "Z2mass") fdata[ihist]->SetMaximum(4.5);
    
    //gPad->SetLogy(true);
    //fdata[ihist]->GetYaxis()->SetRangeUser(0.1,1000);
    //fdata[ihist]->GetYaxis()->SetRangeUser(0.001,100);
    
    shistos[ihist]->Draw("hist,same");
    TH1D fmc_sum = *mc_sum;
    fmc_sum.SetFillColor(kBlack);
    fmc_sum.SetLineColor(kBlack);
    fmc_sum.SetFillStyle(3018);
    fmc_sum.SetMarkerStyle(0);
    fmc_sum.Draw("e2,same");
    fdata[ihist]->Draw("ep,same");
    leg->Draw();
    
    TPaveText *count_tag = new TPaveText(.25,.85,.4,.98,"NDC");
    count_tag->AddText(Form("Data: %.0f",fndata));
    count_tag->AddText(Form("S: %.2f#pm%.2f",fnsig,sqrt(fsig_err)));
    count_tag->AddText(Form("B: %.2f#pm%.2f",fnbkg,sqrt(fbkg_err)));
    count_tag->SetFillStyle(0);
    count_tag->SetBorderSize(0);
    count_tag->SetTextSize(0.05);
    count_tag->Draw();
    
    cv->Update();
    //cv->Print("PlotsVBFnoDjetCutNetwork/histo_"+histos_name[ihist]+".png");
    cv->Print("mass4l_allNets/histo_"+histos_name[ihist]+".png");
    
  }      
  std::cout<<Form("S: %.2f#pm%.2f, B: %.2f#pm%.2f",fnsig,sqrt(fsig_err),fnbkg,sqrt(fbkg_err))<<std::endl;
      
  return;
}//End plot function


void plotObjectsProperties(void){
  gSystem->Exec("mkdir -p PlotsVBFnoDjet");
  
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
  std::vector<TString> DYJets;
  std::vector<TString> QCD;
  std::vector<TString> ZZJJ;
    

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
  DYJets.push_back("histos4mu_25ns/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  ggH.push_back("histos4mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  HJJ.push_back("histos4mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4mu_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4mu_25ns/output_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8.root");
  ttH.push_back("histos4mu_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  ttX.push_back("histos4mu_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  ttX.push_back("histos4mu_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  ttX.push_back("histos4mu_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  VBF.push_back("histos4mu_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  VH.push_back("histos4mu_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos4mu_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  WX.push_back("histos4mu_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  WX.push_back("histos4mu_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos4mu_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  WX.push_back("histos4mu_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos4mu_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  VH.push_back("histos4mu_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  qqZZ.push_back("histos4mu_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos4mu_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  ZZZ.push_back("histos4mu_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  ZZJJ.push_back("histos4mu_25ns/output_ZZJJTo4L_EWK_13TeV-madgraph-pythia8.root");
  
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
  DYJets.push_back("histos4e_25ns/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");           
  ggH.push_back("histos4e_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");  
  HJJ.push_back("histos4e_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");               
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");              
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");                 
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");             
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos4e_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos4e_25ns/output_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8.root");
  ttH.push_back("histos4e_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  ttX.push_back("histos4e_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  ttX.push_back("histos4e_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  ttX.push_back("histos4e_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  VBF.push_back("histos4e_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  VH.push_back("histos4e_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos4e_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  WX.push_back("histos4e_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  WX.push_back("histos4e_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos4e_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  WX.push_back("histos4e_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos4e_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  VH.push_back("histos4e_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  qqZZ.push_back("histos4e_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos4e_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  ZZZ.push_back("histos4e_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  ZZJJ.push_back("histos4e_25ns/output_ZZJJTo4L_EWK_13TeV-madgraph-pythia8.root");
  
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
  DYJets.push_back("histos2e2mu_25ns/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  ggH.push_back("histos2e2mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  HJJ.push_back("histos2e2mu_25ns/output_GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root");
  ggZZ.push_back("histos2e2mu_25ns/output_GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8.root");
  QCD.push_back("histos2e2mu_25ns/output_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8.root");
  ttH.push_back("histos2e2mu_25ns/output_ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  ttX.push_back("histos2e2mu_25ns/output_TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  ttX.push_back("histos2e2mu_25ns/output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root");
  ttX.push_back("histos2e2mu_25ns/output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  VBF.push_back("histos2e2mu_25ns/output_VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
  VH.push_back("histos2e2mu_25ns/output_WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  VH.push_back("histos2e2mu_25ns/output_WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WWTo2L2Nu_13TeV-powheg.root");
  WX.push_back("histos2e2mu_25ns/output_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  WX.push_back("histos2e2mu_25ns/output_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  VH.push_back("histos2e2mu_25ns/output_ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root");
  qqZZ.push_back("histos2e2mu_25ns/output_ZZTo2L2Nu_13TeV_powheg_pythia8.root");
  qqZZ.push_back("histos2e2mu_25ns/output_ZZTo4L_13TeV_powheg_pythia8.root");
  ZZZ.push_back("histos2e2mu_25ns/output_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root");
  ZZJJ.push_back("histos2e2mu_25ns/output_ZZJJTo4L_EWK_13TeV-madgraph-pythia8.root");
  
  std::vector< std::vector<TString> > samples;
  std::vector<TString> tags;

  samples.push_back(Data);	tags.push_back("Data");
  samples.push_back(VBF);	tags.push_back("VBF");
  //samples.push_back(ggH);	tags.push_back("ggH");
  samples.push_back(HJJ);	tags.push_back("HJJ");
  samples.push_back(qqZZ);	tags.push_back("qqZZ");
  samples.push_back(ggZZ);	tags.push_back("ggZZ");
  samples.push_back(ttH);	tags.push_back("ttH");
  samples.push_back(VH);	tags.push_back("VH");
  samples.push_back(ttX);	tags.push_back("ttX");
  samples.push_back(WX);	tags.push_back("WX");
  samples.push_back(ZZZ);	tags.push_back("ZZZ");
  samples.push_back(QCD);	tags.push_back("QCD");
  samples.push_back(DYJets);	tags.push_back("DYJets");
  samples.push_back(ZZJJ);	tags.push_back("ZZ+2jets");

  loopSamples(tags, samples);
     
  //Ends program
  return;
}
