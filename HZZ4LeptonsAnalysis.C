  //By Miqueias M. A.

  #include <time.h>
//  #include <op.h>
  #define HZZ4LeptonsAnalysis_cxx
  #include "HZZ4LeptonsAnalysis.h"
  #include "TStyle.h"
  #include "TGraph.h"
  #include "TLine.h"
  #include "TH1D.h"
  #include "TH2D.h"
  #include "TCanvas.h"
  #include "THStack.h"
  #include "TLorentzVector.h"
  #include "TMath.h"
  #include "TVector3.h"
  #include "TTree.h"
  #include <stdio.h>
  #include <stdlib.h>
  #include <iostream>
  #include <string>
  #include <sstream>
  #include <vector>
  #define prad 3.141592654
  #define angle 180/3.141592654
  #define mass_Z 91.1876
  #define Lumi7TeV 5.051   //1/fb
  #define Lumi8TeV 19.712  //1/fb
//  #define cs_ggZZ_4l_7TeV 
//  #define cs_ggZZ_4l_8TeV
//  #define cs_ggZZ_2l2l_7TeV 3.48
//  #define cs_ggZZ_2l2l_8TeV 12.03
//  #define cs_qqZZ_7TeV
//  #define cs_qqZZ_8TeV
//  #define cs_ZX 

  using std::cout;
  using std::cin;
  using std::endl;


  //Z 4vector function
  TLorentzVector Z(Float_t lep1_pt, Float_t lep1_eta, Float_t lep1_phi, Float_t lep1_mass, Float_t lep2_pt, Float_t lep2_eta, Float_t lep2_phi, Float_t lep2_mass)
  {
    TLorentzVector el1, el2, Z_Candidate;
    el1.SetPtEtaPhiE(lep1_pt, lep1_eta, lep1_phi, lep1_mass);
    el2.SetPtEtaPhiE(lep2_pt, lep2_eta, lep2_phi, lep2_mass);
    Z_Candidate = el1 + el2;
    return Z_Candidate;
  }

  //Leptons 4vector function
  TLorentzVector q_vetor(Float_t lep_pt, Float_t lep_eta, Float_t lep_phi, Float_t lep_mass)
  {
    TLorentzVector qp;
    qp.SetPtEtaPhiE(lep_pt, lep_eta, lep_phi, lep_mass);
    return qp;
  }


void HZZ4LeptonsAnalysis::Loop()
{
    //--------- Setup for grahpics appearance ---------------
    TStyle* myStyle = new TStyle("myStyle","My ROOT Style");
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetPadBorderMode(0);
    myStyle->SetPadColor(0);
    myStyle->SetCanvasColor(0);
    myStyle->SetTitleColor(0);
    myStyle->SetStatColor(0);
    myStyle->SetTitleSize(0.05,"x");
    myStyle->SetTitleColor(1,"x");
    myStyle->SetTitleOffset(1,"x");
    myStyle->SetLabelSize(0.05,"x");
    myStyle->SetTitleSize(0.05,"y");
    myStyle->SetTitleColor(1,"y");
    myStyle->SetTitleOffset(1.2,"y");
    myStyle->SetLabelSize(0.05,"y");
    //myStyle->SetHistLineColor(kBlue);
    //myStyle->SetPadLeftMargin(0.15);	

    gROOT->SetStyle("myStyle");
    //-------------------------------------------------------
   
    if (fChain == 0) return;
//#######################################################  Set-Up  ###################################################################

    //Some variables
    Int_t ent=-1; 
    Int_t H_4e=0, H_4u=0, H_hbd=0, H_ind=0;
    
    //Range and #bins of histograms
    Int_t z1_xi=40, z1_xf=120, z2_xi=12, z2_xf=120, hg_xi=100, hg_xf=850, pT_i=0, pT_f=200;
    Float_t z1_bin=40, z2_bin=54, hg_bin=75, pT_bin=100;
   
   
    //----------------------- Histograms for distribution pT, Eta, Phi --------------------------------------------------------
    TH1D* hPtEle = new TH1D("hPtEle","hPtEle",150,0,150);   hPtEle->SetFillColor(4);   hPtEle->GetXaxis()->SetTitle("p_{T}(GeV)");
    TH1D* hEtaEle = new TH1D("hEtaEle","hEtaEle",100,-4,4);   hEtaEle->SetFillColor(4);   hEtaEle->GetXaxis()->SetTitle("#eta");
    TH1D* hPhiEle = new TH1D("hPhiEle","hPhiEle",100,-4,4);   hPhiEle->SetFillColor(4);   hPhiEle->GetXaxis()->SetTitle("#phi");
    TH1D* hPtMu = new TH1D("hPtMu","hPtMu",150,0,150);   hPtMu->SetFillColor(4);   hPtMu->GetXaxis()->SetTitle("p_{T}(GeV)");
    TH1D* hEtaMu = new TH1D("hEtaMu","hEtaMu",100,-4,4);   hEtaMu->SetFillColor(4);   hEtaMu->GetXaxis()->SetTitle("#eta");
    TH1D* hPhiMu = new TH1D("hPhiMu","hPhiMu",100,-4,4);   hPhiMu->SetFillColor(4);   hPhiMu->GetXaxis()->SetTitle("#phi");
    //-----------------------------------------------------------------------------------------------------------------------------
   
    //-------------------- Histograms for Topologic Variables -------------------------------------------------
    TH1D* hcos0star = new TH1D("hcos0star","cos #theta^{*}",20,-1,1);
    TH1D* hcos01 = new TH1D("hcos01","cos #theta_{1}",20,-1,1);
    TH1D* hcos02 = new TH1D("hcos02","cos #theta_{2}",20,-1,1);
    TH1D* hphi = new TH1D("hphi","#Phi",32,-3.14,3.14);
    TH1D* hphi1 = new TH1D("hphi1","#Phi_{1}",32,-3.14,3.14);
   
    //-------------------- Histograms for particles reconstruction --------------------------------------------
    TString type = "Eventos/", unit = "GeV", scale;
    std::stringstream sac;
    
    //---------------------- Z's ---------------------
    sac << (z1_xf - z1_xi)/z1_bin;   scale = sac.str();
    TH1D* h_Z1ee = new TH1D("h_Z1ee","Z1 #rightarrow 2e",z1_bin,z1_xi,z1_xf);   h_Z1ee->SetFillColor(0);   
    h_Z1ee->GetXaxis()->SetTitle("m_{2e} (GeV)");   h_Z1ee->GetYaxis()->SetTitle(type+scale+unit);
    TH1D* h_Z1uu = new TH1D("h_Z1uu","Z1 #rightarrow 2#mu",z1_bin,z1_xi,z1_xf);   h_Z1uu->SetFillColor(0);
    h_Z1uu->GetXaxis()->SetTitle("m_{2#mu} (GeV)");   h_Z1uu->GetYaxis()->SetTitle(type+scale+unit);
    TH1D* h_Z1 = new TH1D("h_Z1","Z1 #rightarrow ll",z1_bin,z1_xi,z1_xf);   h_Z1->SetFillColor(0);
    h_Z1->GetXaxis()->SetTitle("m_{ll} (GeV)");   h_Z1->GetYaxis()->SetTitle(type+scale+unit);
    
    sac.str("");   sac << (z2_xf - z2_xi)/z2_bin;   scale = sac.str();
    TH1D* h_Z2ee = new TH1D("h_Z2ee","Z2 #rightarrow 2e",z2_bin,z2_xi,z2_xf);   h_Z2ee->SetFillColor(0);   
    h_Z2ee->GetXaxis()->SetTitle("m_{2e} (GeV)");   h_Z2ee->GetYaxis()->SetTitle(type+scale+unit);
    TH1D* h_Z2uu = new TH1D("h_Z2uu","Z2 #rightarrow 2#mu",z2_bin,z2_xi,z2_xf);   h_Z2uu->SetFillColor(0);   
    h_Z2uu->GetXaxis()->SetTitle("m_{2#mu} (GeV)");   h_Z2uu->GetYaxis()->SetTitle(type+scale+unit);
    TH1D* h_Z2 = new TH1D("h_Z2","Z2 #rightarrow ll",z2_bin,z2_xi,z2_xf);   h_Z2->SetFillColor(0);
    h_Z2->GetXaxis()->SetTitle("m_{ll} (GeV)");   h_Z2->GetYaxis()->SetTitle(type+scale+unit);

    //-------------------- Higgs's mass -------------------
    sac.str("");   sac << (hg_xf - hg_xi)/hg_bin;   scale = sac.str();
    TH1D* h_Higgs_4u = new TH1D("h_Higgs_4u","H #rightarrow 4#mu",hg_bin, hg_xi, hg_xf);
    h_Higgs_4u->SetFillColor(kBlue);
    h_Higgs_4u->GetXaxis()->SetTitle("m_{4#mu} (GeV)");
    h_Higgs_4u->GetYaxis()->SetTitle(type+scale+unit);
   
    TH1D* h_Higgs_4e = new TH1D("h_Higgs_4e","H #rightarrow 4e",hg_bin, hg_xi, hg_xf);
    h_Higgs_4e->SetFillColor(kBlack);
    h_Higgs_4e->GetXaxis()->SetTitle("m_{4e} (GeV)");
    h_Higgs_4e->GetYaxis()->SetTitle(type+scale+unit);
   
    TH1D* h_Higgs_hbd = new TH1D("h_Higgs_hbd","H #rightarrow 2e2#mu/2#mu2e",hg_bin, hg_xi, hg_xf);
    h_Higgs_hbd->SetFillColor(kGreen);
    h_Higgs_hbd->GetXaxis()->SetTitle("m_{2e2#mu/2#mu2e} (GeV)");
    h_Higgs_hbd->GetYaxis()->SetTitle(type+scale+unit);
    
    TH1D* h_Higgs_ind = new TH1D("h_Higgs_correct","H #rightarrow 4l",hg_bin, hg_xi, hg_xf);
    //h_Higgs_ind->SetFillColor();
    h_Higgs_ind->GetXaxis()->SetTitle("m_{4l} (GeV)");
    h_Higgs_ind->GetYaxis()->SetTitle(type+scale+unit);
    
    TH2D* h_Higgs_pT = new TH2D("h_Higgs_pT","H #rightarrow 4l (p_{T})",hg_bin, hg_xi, hg_xf, pT_bin, pT_i, pT_f);
    //h_Higgs_ind->SetFillColor();
    h_Higgs_pT->GetXaxis()->SetTitle("m_{4l} (GeV)");
    h_Higgs_pT->GetYaxis()->SetTitle("p_{T}^{4l} (2GeV^{-1})");
    
    TH1D* h_Higgs_mela = new TH1D("h_Higgs_mela","H #rightarrow 4l (MELA)",hg_bin, hg_xi, hg_xf);
    h_Higgs_mela->GetXaxis()->SetTitle("m_{4l} (GeV)");
    h_Higgs_mela->GetYaxis()->SetTitle(type+scale+unit);
    
    TH2D* h_Higgs_mela1 = new TH2D("h_Higgs_mela1","H #rightarrow 4l (MELA) cos #theta *",hg_bin, hg_xi, hg_xf, 200, -1, 1);
    h_Higgs_mela1->GetXaxis()->SetTitle("m_{4l} (GeV)");
    h_Higgs_mela1->GetYaxis()->SetTitle("cos #theta *");
    
    TH2D* h_Higgs_mela2 = new TH2D("h_Higgs_mela2","H #rightarrow 4l (MELA) cos #theta",hg_bin, hg_xi, hg_xf, 200, -1, 1);
    h_Higgs_mela2->GetXaxis()->SetTitle("m_{4l} (GeV)");
    h_Higgs_mela2->GetYaxis()->SetTitle("cos #theta");
    
    TH2D* h_Higgs_mela3 = new TH2D("h_Higgs_mela3","H #rightarrow 4l (MELA) #Phi",hg_bin, hg_xi, hg_xf, 620, -3.14, 3.14);
    h_Higgs_mela3->GetXaxis()->SetTitle("m_{4l} (GeV)");
    h_Higgs_mela3->GetYaxis()->SetTitle("#Phi");
    
    TH2D* h_Higgs_mela4 = new TH2D("h_Higgs_mela4","H #rightarrow 4l (MELA) #Phi_{1}",hg_bin, hg_xi, hg_xf, 620, -3.14, 3.14);
    h_Higgs_mela4->GetXaxis()->SetTitle("m_{4l} (GeV)");
    h_Higgs_mela4->GetYaxis()->SetTitle("#Phi_{1}");
         
    //------------------------------ Z mass scatter -------------------------------------------------------
    TH2D* Zmass_scatter = new TH2D("Zmass_scatter","Z Mass Scatter", 80, 0, 200, 80, 0, 200);
    Zmass_scatter->GetXaxis()->SetTitle("Z on shell mass");
    Zmass_scatter->GetYaxis()->SetTitle("Z off shell mass");
    Zmass_scatter->SetMarkerStyle(6);
    //-----------------------------------------------------------------------------------------------------
   
    
    //For MELA
    std::vector<float> vcos_star, vcos_theta1, vcos_theta2, vphi, vphi1, vhiggs_mass;
   
    TTree *tree = new TTree("VBF","");
    Double_t Zon_mass, Zoff_mass, Higgs_mass, cos_star, cos_theta1, cos_theta2, Phi, Phi1;
    Float_t RECO_PARTICLE[4][4]; Int_t lep[4];
    TLorentzVector q11, q12, q21, q22;
    std::vector<Double_t> Zs_mass;
    //tree->Branch("Zs_mass",&Zs_mass);
    tree->Branch("Zon_mass",&Zon_mass);  //nome da brach, variavel usada, tipo da variavel (D:double, F:float, I:inteiro)
    tree->Branch("Zoff_mass",&Zoff_mass);
    tree->Branch("Higgs_mass",&Higgs_mass);
    //tree->Branch("cos0star",&cos_star);
    //tree->Branch("costheta1",&cos_theta1);
    //tree->Branch("costheta2",&cos_theta2);
    //tree->Branch("phi",&Phi);
    //tree->Branch("phi1",&Phi1);
    tree->Branch("lep",&lep,"lep[4]/I");
    tree->Branch("RECO_PARTICLE",&RECO_PARTICLE,"RECO_PARTICLE[4][4]/F");
    tree->SetDirectory(0);
    
   
    cout<< "\n============================================================================" << endl;
    
    Long64_t nentries = fChain->GetEntries();
    cout << "Analysing "<< nentries<<" events..." << endl;
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries ;jentry++) 
     {
       //cout<< "jentry= " << jentry << endl;
       ent = ent + 1;
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break; 
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       // if (Cut(ientry) < 0) continue;
           
//##########################################################  Loop analyzer of user #######################################################  
       
      ////////// pT, Eta, Phi /////////////
      for ( Int_t iel=0; iel < RECO_NELE ; iel++ )
       {if(RECOELE_CHARGE[iel] != -999){
	 hPtEle->Fill(RECOELE_PT[iel]);
	 hEtaEle->Fill(RECOELE_ETA[iel]);
	 hPhiEle->Fill(RECOELE_PHI[iel]);}
       }
     
      ////////// pT, Eta, Phi ///////////////   
      for ( Int_t imu=0 ; imu < RECO_NMU ; imu++ )
       {if(RECOMU_CHARGE[imu] != -999){
	 hPtMu->Fill(RECOMU_PT[imu]);
	 hEtaMu->Fill(RECOMU_ETA[imu]);
	 hPhiMu->Fill(RECOMU_PHI[imu]);}
       }

       
      ////////////////////// Loop on leptons pairs //////////////////////
//      cout<<"********** Run "<<jentry<<" *************"<<endl;
//================================================================== On shell Z ===============================================================================
      Int_t minimo_e=(z1_xf-z1_xi), minimo_u=(z1_xf-z1_xi), e1=-1, e2=-1, e3=-1, e4=-1, u1=-1, u2=-1, u3=-1, u4=-1;
      Int_t lep1=-1, lep2=-1;
      Double_t mass_pele=0, mass_pmu=0;
      TLorentzVector Z_uu, Z_ee, Zuu_off, Zee_off;
      TLorentzVector Z_on, Z_off;
      q11.SetPtEtaPhiE(-99,0,0,0), q12.SetPtEtaPhiE(-99,0,0,0), q21.SetPtEtaPhiE(-99,0,0,0), q22.SetPtEtaPhiE(-99,0,0,0);
      
      Z_on.SetPtEtaPhiE(0.,0.,0.,0.);   Z_off.SetPtEtaPhiE(0.,0.,0.,0.);
      for( Int_t i = 0; i < RECO_NMU; i++ ) 
      for( Int_t j=i+1; j < RECO_NMU; j++ )
       	{	  	  
	 if( RECOMU_CHARGE[i]*RECOMU_CHARGE[j] < 0 && RECOMU_CHARGE[i] != -999 && RECOMU_CHARGE[j] != -999 )
	  {
	   if( RECOMU_PT[i] <= 5 || RECOMU_PT[j] <= 5 || fabs(RECOMU_ETA[i]) >= 2.4 || fabs(RECOMU_ETA[j]) >= 2.4 ) continue;
	  	    	
	   mass_pmu = Z(RECOMU_PT[i], RECOMU_ETA[i], RECOMU_PHI[i], RECOMU_E[i], RECOMU_PT[j], RECOMU_ETA[j], RECOMU_PHI[j], RECOMU_E[j]).M();
	   if( fabs(mass_pmu - mass_Z) < minimo_u )
	    {
	     //Verification of closeness Z candidate mass 
	     minimo_u = fabs(mass_pmu - mass_Z);         u1 = i; u2 = j;
	     Z_uu = Z(RECOMU_PT[i], RECOMU_ETA[i], RECOMU_PHI[i], RECOMU_E[i], RECOMU_PT[j], RECOMU_ETA[j], RECOMU_PHI[j], RECOMU_E[j]);
	    }
	  }
	}
	 
      for( Int_t i = 0; i < RECO_NELE; i++ ) 
      for( Int_t j=i+1; j < RECO_NELE; j++ )
       	{
	 if( RECOELE_CHARGE[i]*RECOELE_CHARGE[j] < 0 && RECOELE_CHARGE[i] !=-999 && RECOELE_CHARGE[j] !=-999 )
	  {
	   if( RECOELE_PT[i] <= 7 || RECOELE_PT[j] <= 7 || fabs(RECOELE_ETA[i]) >= 2.5 || fabs(RECOELE_ETA[j]) >= 2.5 ) continue ;	
	   mass_pele = Z(RECOELE_PT[i], RECOELE_ETA[i], RECOELE_PHI[i], RECOELE_E[i], RECOELE_PT[j], RECOELE_ETA[j], RECOELE_PHI[j], RECOELE_E[j]).M();
	   if( fabs(mass_pele - mass_Z) < minimo_e ) 
	    { 
	     //Verification of closeness Z candidate mass 
	     minimo_e = fabs(mass_pele - mass_Z);        e1 = i; e2 = j;
	     Z_ee = Z(RECOELE_PT[i], RECOELE_ETA[i], RECOELE_PHI[i], RECOELE_E[i], RECOELE_PT[j], RECOELE_ETA[j], RECOELE_PHI[j], RECOELE_E[j]);
	    }  
	  }
	}	

//============================================ Define Z1 and searches Z2 ======================================================	
     
     //************* Case 1: On shell Z, built with muons, more close to nominal Z mass *********************
     if( fabs(Z_uu.M()-mass_Z)<fabs(Z_ee.M()-mass_Z) )
      {
	Z_on = Z_uu;          h_Z1->Fill(Z_on.M());		  //Finely save candidate Z1
	lep[0] = (RECOMU_CHARGE[u1]>0)?13:-13;		lep[1] = (RECOMU_CHARGE[u2]>0)?13:-13;
	RECO_PARTICLE[0][0] = RECOMU_PT[u1];		RECO_PARTICLE[1][0] = RECOMU_PT[u2];
	RECO_PARTICLE[0][1] = RECOMU_ETA[u1];		RECO_PARTICLE[1][1] = RECOMU_ETA[u2];
	RECO_PARTICLE[0][2] = RECOMU_PHI[u1];		RECO_PARTICLE[1][2] = RECOMU_PHI[u2];
	RECO_PARTICLE[0][3] = RECOMU_E[u1];		RECO_PARTICLE[1][3] = RECOMU_E[u2];
	if( RECOMU_CHARGE[u1] < 0 && RECOMU_CHARGE[u2] > 0 ){ lep1 = u1; lep2 = u2; }
	else if( RECOMU_CHARGE[u1] > 0 && RECOMU_CHARGE[u2] < 0 ){ lep1 = u2; lep2 = u1; }
	q11.SetPtEtaPhiE(RECOMU_PT[lep1], RECOMU_ETA[lep1], RECOMU_PHI[lep1], RECOMU_E[lep1]);
	q12.SetPtEtaPhiE(RECOMU_PT[lep2], RECOMU_ETA[lep2], RECOMU_PHI[lep2], RECOMU_E[lep2]);
     //--------------------------- Selection of 1st lepton of off shell Z -----------------------------------------------------
      Float_t max_ptu=0; 
      for( Int_t i = 0; i < RECO_NMU; i++ ) 
     	{
	 if( i != u1 && i != u2 && RECOMU_CHARGE[i] !=-999 )
	  { 
           if( RECOMU_PT[i] <= 5 || fabs(RECOMU_ETA[i]) >= 2.4 ) continue;	  
           if( RECOMU_PT[i] > max_ptu )
	    {  
	     //Searches 1st maximum pT muon 
	     max_ptu = RECOMU_PT[i];  u3 = i;
	    }
	  }
	}
      Float_t max_pte=0;
      for(Int_t i=0; i< RECO_NELE; i++)
        {
         if( (RECOELE_PT[i] <= 7 || fabs(RECOELE_ETA[i]) >= 2.5) && RECOELE_CHARGE[i] ==-999 ) continue;
         if( RECOELE_PT[i] > max_pte ) 
	  { 
	   // Searches 1st maximum pT electron 
	   max_pte = RECOELE_PT[i];  e3 = i; 
	  }
	}
	
      //-------------------- Selection of 2nd lepton of off shell Z --------------------------------------------------------------------------
      max_ptu=0;
      for( Int_t i = 0; i < RECO_NMU; i++ ) 
        {
	 if( u3 != -1 && RECOMU_CHARGE[i]*RECOMU_CHARGE[u3] < 0 && i!=u1 && i!=u2 && i!=u3 )
	  {
           if( (RECOMU_PT[i] <= 5 || fabs(RECOMU_ETA[i]) >= 2.4) && RECOMU_CHARGE[i] ==-999 ) continue;
	   if( RECOMU_PT[i] > max_ptu )
	    {  
	     // Searches 2nd maximum pT muon 
	     max_ptu = RECOMU_PT[i];  u4 = i;
	    }
	  }
        }
      max_pte=0;
      for(Int_t i=0; i<RECO_NELE; i++)
        {
	 if( e3 != -1 && RECOELE_CHARGE[i]*RECOELE_CHARGE[e3] < 0 && i!=e3 )
	  {
	   if( (RECOELE_PT[i] <= 7 || fabs(RECOELE_ETA[i]) >= 2.5) && RECOELE_CHARGE[i] ==-999 ) continue;
	   if( RECOELE_PT[i] > max_pte ) 
	    { 
	     // Searches 2nd maximum pT electron 
	     max_pte = RECOELE_PT[i];  e4 = i;
	    }
	  }
        }
        
        //------- Choice of Z2 Case 1 -------------------
        if( (u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && RECOMU_PT[u3]+RECOMU_PT[u4]>RECOELE_PT[e3]+RECOELE_PT[e4]) || e3==-1 || e4==-1 )
	{ 
	  Z_off = Z(RECOMU_PT[u3], RECOMU_ETA[u3], RECOMU_PHI[u3], RECOMU_E[u3], RECOMU_PT[u4], RECOMU_ETA[u4], RECOMU_PHI[u4], RECOMU_E[u4]);
	  h_Z2->Fill(Z_off.M()); 
	  //-------------------------------------------------------------------------------------------------------------------
	  if( RECOMU_CHARGE[u3] < 0 && RECOMU_CHARGE[u4] > 0 ){ lep1 = u3; lep2 = u4; }
	  else if( RECOMU_CHARGE[u3] > 0 && RECOMU_CHARGE[u4] < 0 ){ lep1 = u4; lep2 = u3; }
	  q21.SetPtEtaPhiE(RECOMU_PT[lep1], RECOMU_ETA[lep1], RECOMU_PHI[lep1], RECOMU_E[lep1]);
	  q22.SetPtEtaPhiE(RECOMU_PT[lep2], RECOMU_ETA[lep2], RECOMU_PHI[lep2], RECOMU_E[lep2]);
	  
	  lep[2] = (RECOMU_CHARGE[u3]>0)?13:-13;	lep[3] = (RECOMU_CHARGE[u4]>0)?13:-13;
	  RECO_PARTICLE[2][0] = RECOMU_PT[u3];		RECO_PARTICLE[3][0] = RECOMU_PT[u4];
	  RECO_PARTICLE[2][1] = RECOMU_ETA[u3];		RECO_PARTICLE[3][1] = RECOMU_ETA[u4];
	  RECO_PARTICLE[2][2] = RECOMU_PHI[u3];		RECO_PARTICLE[3][2] = RECOMU_PHI[u4];
	  RECO_PARTICLE[2][3] = RECOMU_E[u3];		RECO_PARTICLE[3][3] = RECOMU_E[u4];
	}
	if( (u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && RECOMU_PT[u3]+RECOMU_PT[u4]<RECOELE_PT[e3]+RECOELE_PT[e4]) || u3==-1 || u4==-1 )
	{
	  Z_off = Z(RECOELE_PT[e3], RECOELE_ETA[e3], RECOELE_PHI[e3], RECOELE_E[e3], RECOELE_PT[e4], RECOELE_ETA[e4], RECOELE_PHI[e4], RECOELE_E[e4]);
	  h_Z2->Fill(Z_off.M());
	  //-------------------------------------------------------------------------------------------------------------------------
	  if( RECOELE_CHARGE[e3] < 0 && RECOELE_CHARGE[e4] > 0 ){ lep1 = e3; lep2 = e4; }
	  else if( RECOELE_CHARGE[e3] > 0 && RECOELE_CHARGE[e4] < 0 ){ lep1 = e4; lep2 = e3; }
	  q21.SetPtEtaPhiE(RECOELE_PT[lep1], RECOELE_ETA[lep1], RECOELE_PHI[lep1], RECOELE_E[lep1]);
	  q22.SetPtEtaPhiE(RECOELE_PT[lep2], RECOELE_ETA[lep2], RECOELE_PHI[lep2], RECOELE_E[lep2]);
	  
	  lep[2] = (RECOELE_CHARGE[e3]>0)?11:-11;	lep[2] = (RECOELE_CHARGE[e4]>0)?11:-11;
	  RECO_PARTICLE[2][0] = RECOELE_PT[e3];		RECO_PARTICLE[3][0] = RECOELE_PT[e4];
	  RECO_PARTICLE[2][1] = RECOELE_ETA[e3];	RECO_PARTICLE[3][1] = RECOELE_ETA[e4];
	  RECO_PARTICLE[2][2] = RECOELE_PHI[e3];	RECO_PARTICLE[3][2] = RECOELE_PHI[e4];
	  RECO_PARTICLE[2][3] = RECOELE_E[e3];		RECO_PARTICLE[3][3] = RECOELE_E[e4];
	}
	
      }//End choice process for Case 1 (muonicZ +~ nominalZ)
      
      
     //************* Case 2: On shell Z, built with electrons, more close to nominal Z mass **********************
     if( fabs(Z_uu.M()-mass_Z)>fabs(Z_ee.M()-mass_Z) )
     {
       Z_on = Z_ee;           h_Z1->Fill(Z_ee.M());
       //------------------------------------------------------------------------------------------------------------------------
       if( RECOELE_CHARGE[e1] < 0 && RECOELE_CHARGE[e2] > 0 ){ lep1 = e1; lep2 = e2; }
       else if( RECOELE_CHARGE[e1] > 0 && RECOELE_CHARGE[e2] < 0 ){ lep1 = e2; lep2 = e1; }
       q11.SetPtEtaPhiE(RECOELE_PT[lep1], RECOELE_ETA[lep1], RECOELE_PHI[lep1], RECOELE_E[lep1]);
       q12.SetPtEtaPhiE(RECOELE_PT[lep2], RECOELE_ETA[lep2], RECOELE_PHI[lep2], RECOELE_E[lep2]);
       
       lep[0] = (RECOELE_CHARGE[e1]>0)?11:-11;		lep[1] = (RECOELE_CHARGE[e2]>0)?11:-11;
       RECO_PARTICLE[0][0] = RECOELE_PT[e1];		RECO_PARTICLE[1][0] = RECOELE_PT[e2];
       RECO_PARTICLE[0][1] = RECOELE_ETA[e1];		RECO_PARTICLE[1][1] = RECOELE_ETA[e2];
       RECO_PARTICLE[0][2] = RECOELE_PHI[e1];		RECO_PARTICLE[1][2] = RECOELE_PHI[e2];
       RECO_PARTICLE[0][3] = RECOELE_E[e1];		RECO_PARTICLE[1][3] = RECOELE_E[e2];
      //------------------------------------ Selection of 1st lepton of off shell Z ---------------------------------------------------
      Float_t max_ptu=0; 
      for( Int_t i = 0; i < RECO_NMU; i++ ) 
     	{
         if( RECOMU_PT[i] <= 5 || fabs(RECOMU_ETA[i]) >= 2.4 ) continue;	  
         if( RECOMU_PT[i] > max_ptu && RECOMU_CHARGE[i] !=-999 )
	  {  
	   // Searchs 1st maximum pT muon 
	   max_ptu = RECOMU_PT[i];  u3 = i;
	  }
	}	
      Float_t max_pte=0;
      for(Int_t i=0; i< RECO_NELE; i++)
        {
         if( i != e1 && i != e2 )
	  {
           if( RECOELE_PT[i] <= 7 || fabs(RECOELE_ETA[i]) >= 2.5 ) continue;
           if( RECOELE_PT[i] > max_pte && RECOELE_CHARGE[i] !=-999 ) 
	    { 
	     // Searches 1st maximum pT electron 
	     max_pte = RECOELE_PT[i];  e3 = i; 
	    }
	  }
	}
	
      //------------------------- Selection of 2nd lepton of off shell Z --------------------------------------------------------------------------
      max_ptu=0;
      for( Int_t i = 0; i < RECO_NMU; i++ ) 
        {
	 if( u3 != -1 && RECOMU_CHARGE[i]*RECOMU_CHARGE[u3] < 0 && i!=u3 )
	  {
           if( RECOMU_PT[i] <= 5 || fabs(RECOMU_ETA[i]) >= 2.4 ) continue;
	   if( RECOMU_PT[i] > max_ptu && RECOMU_CHARGE[i] !=-999 )
	    {  
	     // Searches 2nd maximum pT muon 
	     max_ptu = RECOMU_PT[i];  u4 = i;
	    }
	  }
        }
      max_pte=0;
      for(Int_t i=0; i<RECO_NELE; i++)
        {
	 if( e3 != -1 && RECOELE_CHARGE[i]*RECOELE_CHARGE[e3] < 0 && i!=e1 && i!=e2 && i!=e3 )
	  {
	   if( RECOELE_PT[i] <= 7 || fabs(RECOELE_ETA[i]) >= 2.5 ) continue;
	   if( RECOELE_PT[i] > max_pte && RECOELE_CHARGE[i] !=-999 ) 
	    { 
	     // Searches 2nd maximum pT electron 
	     max_pte = RECOELE_PT[i];  e4 = i;
	    }
	  }
        }
        
        //------- Choice of off shell Z Case 2 -------------------
        if( (u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && RECOMU_PT[u3]+RECOMU_PT[u4]>RECOELE_PT[e3]+RECOELE_PT[e4]) || e3==-1 || e4==-1 )
	{ 
	  Z_off = Z(RECOMU_PT[u3], RECOMU_ETA[u3], RECOMU_PHI[u3], RECOMU_E[u3], RECOMU_PT[u4], RECOMU_ETA[u4], RECOMU_PHI[u4], RECOMU_E[u4]);
	  h_Z2->Fill(Z_off.M()); 
	  //-------------------------------------------------------------------------------------------------------------------
	  if( RECOMU_CHARGE[u3] < 0 && RECOMU_CHARGE[u4] > 0 ){ lep1 = u3; lep2 = u4; }
	  else if( RECOMU_CHARGE[u3] > 0 && RECOMU_CHARGE[u4] < 0 ){ lep1 = u4; lep2 = u3; }
	  q21.SetPtEtaPhiE(RECOMU_PT[lep1], RECOMU_ETA[lep1], RECOMU_PHI[lep1], RECOMU_E[lep1]);
	  q22.SetPtEtaPhiE(RECOMU_PT[lep1], RECOMU_ETA[lep1], RECOMU_PHI[lep1], RECOMU_E[lep1]);
	  
	  lep[2] = (RECOMU_CHARGE[u3]>0)?13:-13;	lep[3] = (RECOMU_CHARGE[u4]>0)?13:-13;
	  RECO_PARTICLE[2][0] = RECOMU_PT[u3];		RECO_PARTICLE[3][0] = RECOMU_PT[u4];
	  RECO_PARTICLE[2][1] = RECOMU_ETA[u3];		RECO_PARTICLE[3][1] = RECOMU_ETA[u4];
	  RECO_PARTICLE[2][2] = RECOMU_PHI[u3];		RECO_PARTICLE[3][2] = RECOMU_PHI[u4];
	  RECO_PARTICLE[2][3] = RECOMU_E[u3];		RECO_PARTICLE[3][3] = RECOMU_E[u4];
	}
	if( (u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && RECOMU_PT[u3]+RECOMU_PT[u4]<RECOELE_PT[e3]+RECOELE_PT[e4]) || u3==-1 || u4==-1 )
	{
	  Z_off = Z(RECOELE_PT[e3], RECOELE_ETA[e3], RECOELE_PHI[e3], RECOELE_E[e3], RECOELE_PT[e4], RECOELE_ETA[e4], RECOELE_PHI[e4], RECOELE_E[e4]);
	  h_Z2->Fill(Z_off.M());
	  //-------------------------------------------------------------------------------------------------------------------------
	  if( RECOELE_CHARGE[e3] < 0 && RECOELE_CHARGE[e4] > 0 ){ lep1 = e3; lep2 = e4; }
	  else if( RECOELE_CHARGE[e3] > 0 && RECOELE_CHARGE[e4] < 0 ){ lep1 = e4; lep2 = e3; }
	  q21.SetPtEtaPhiE(RECOELE_PT[lep1], RECOELE_ETA[lep1], RECOELE_PHI[lep1], RECOELE_E[lep1]);
	  q22.SetPtEtaPhiE(RECOELE_PT[lep2], RECOELE_ETA[lep2], RECOELE_PHI[lep2], RECOELE_E[lep2]);
	  
	  lep[2] = (RECOELE_CHARGE[e3]>0)?11:-11;	lep[3] = (RECOELE_CHARGE[e4]>0)?11:-11;
	  RECO_PARTICLE[2][0] = RECOELE_PT[e3];		RECO_PARTICLE[3][0] = RECOELE_PT[e4];
	  RECO_PARTICLE[2][1] = RECOELE_ETA[e3];	RECO_PARTICLE[3][1] = RECOELE_ETA[e4];
	  RECO_PARTICLE[2][2] = RECOELE_PHI[e3];	RECO_PARTICLE[3][2] = RECOELE_PHI[e4];
	  RECO_PARTICLE[2][3] = RECOELE_E[e3];		RECO_PARTICLE[3][3] = RECOELE_E[e4];
	}
	
      }//End choice process for Case 2 (electronicZ +~ nominalZ)
//==================================================================================================================================================================	
//End selection process of Z's//       
//=============================================== Higgs candidates selection ========================================================	
	if( u1 != -1 && u2 != -1 && u3 != -1 && u4 != -1 && (Z_on+Z_off).M()>=hg_xi && (Z_on+Z_off).M()<=hg_xf )
	{ h_Higgs_4u->Fill((Z_on + Z_off).M()); H_4u = H_4u + 1; Zmass_scatter->Fill(Z_on.M(),Z_off.M()); }//4 muons
	
	if( e1 != -1 && e2 != -1 && e3 != -1 && e4 != -1 && (Z_on + Z_off).M() >= hg_xi && (Z_on + Z_off).M() <= hg_xf )
	{ h_Higgs_4e->Fill((Z_on + Z_off).M()); H_4e = H_4e + 1; Zmass_scatter->Fill(Z_on.M(), Z_off.M()); }//4 eletrons
	
	if( ((u1 != -1 && u2 != -1 && e3 != -1 && e4 != -1) || (u3 != -1 && u4 != -1 && e1 != -1 && e2 != -1)) && (Z_on + Z_off).M() >= hg_xi && (Z_on + Z_off).M() <= hg_xf)
	{ h_Higgs_hbd->Fill((Z_on + Z_off).M()); H_hbd = H_hbd +1; Zmass_scatter->Fill(Z_on.M(), Z_off.M()); }//2 muons e 2 eletrons

        //  All <<-- This is the most correct now!!
	if( (Z_on+Z_off).M()>=hg_xi && (Z_on+Z_off).M()<=hg_xf )
	 {
	  H_ind = H_ind + 1;
          h_Higgs_ind->Fill( (Z_on+Z_off).M() );
	  
	  h_Higgs_pT->Fill( (Z_on+Z_off).M(),(Z_on+Z_off).Pt() );
	  //For Tree						
	  Higgs_mass = (Z_on + Z_off).M();			
/*	  cos_star = cos0star(qHiggs, Z_on, Z_off);		h_Higgs_mela1->Fill(Higgs_mass,cos_star);		 
	  cos_theta1 = theta1(Z_on, q11, q12, q21, q22);	h_Higgs_mela2->Fill(Higgs_mass,cos_theta1);
	  cos_theta2 = theta2(Z_on, q11, q12, q21, q22);	
	  Phi = phi(qHiggs, q11, q12, q21, q22);		h_Higgs_mela3->Fill(Higgs_mass,Phi);
	  Phi1 = phi1(q11, q12);				h_Higgs_mela4->Fill(Higgs_mass,Phi1);
*/	 }
	 
	 /*Save a Tree in the Root file created for store the other histograms;
	  For save vector, the way is the same for handle to any vector (I mean, can be used a c++ vector)*/
	 Zs_mass.clear();
	 Zon_mass = Z_on.M(); Zoff_mass = Z_off.M();
	 Zs_mass.push_back(Z_on.M());  Zs_mass.push_back(Z_off.M());
	 tree->Fill();
	 
	 
     }//End general loop
      //fclose(Z_masses);
      int Sum_cos_star=0, Sum_Phi1=0;
      for(int r=1; r<=20; r++){
	Sum_cos_star = Sum_cos_star + hcos0star->GetBinContent(r);
	Sum_Phi1 = Sum_Phi1 + hphi1->GetBinContent(r);
      }
//#######################################################  Wrap-Up  #################################################################### 
  
    cout<< "Entries analised: " << ent << endl; //Observar: deve mudar quando alterados os arquivos adicionados ao TChain no arquivo HZZ4LeptonsAnalysis.h!
    cout<< "----------------------------------------------------------------------------" << endl;
    //cout<< "Z1_uu.size= " << Z1_uu.size() << "   |Z2_uu.size= " << Z2_uu.size() << endl;
    //cout<< "Z1_ee.size= " << Z1_ee.size() << "   |Z2_ee.size= " << Z2_ee.size() << endl;
    cout<< "H_4u= " << H_4u << "  |H_4e= " << H_4e << "  |H_hbd= " << H_hbd << endl;
    cout<< "cos0star_medio: "<< Sum_cos_star/20. << "  Phi1_medio: " << Sum_Phi1/20. << endl;
  
  
    /////////////////// Higgs histograms for all channels sum /////////////////////////
    THStack Higgs_stack("Higgs","Soma dos Canais (Fit desabilitado)");
    Higgs_stack.Add(h_Higgs_4e);
    Higgs_stack.Add(h_Higgs_4u);
    Higgs_stack.Add(h_Higgs_hbd);
    
   
    TH1D* hsum_fit = new TH1D("hsum_fit","Soma dos Canais (Fit habilitado)",hg_bin, hg_xi, hg_xf);
    hsum_fit->GetXaxis()->SetTitle("m_{4l} (GeV)");
    hsum_fit->GetYaxis()->SetTitle(type+scale+unit);
    hsum_fit->Add(h_Higgs_4e);
    hsum_fit->Add(h_Higgs_4u);
    hsum_fit->Add(h_Higgs_hbd);
    ///////////////////////////////////////  
    //h_Higgs_all_pT->Add(h_Higgs_4e_pT);
    //h_Higgs_all_pT->Add(h_Higgs_4mu_pT);
    //h_Higgs_all_pT->Add(h_Higgs_hbd_pT);
      
      //Salve informations in rootple
      TFile f("My_MCRECO.root","RECREATE");
      hPtEle->Write();
      hEtaEle->Write();
      hPhiEle->Write();
      hPtMu->Write();
      hEtaMu->Write();
      hPhiMu->Write();
      h_Z1->Write();
      h_Z2->Write();
      Zmass_scatter->Write();
      h_Higgs_4u->Write();
      h_Higgs_4e->Write();
      h_Higgs_hbd->Write();
      hsum_fit->Write();
      Higgs_stack.Write();
      h_Higgs_ind->Write();
      h_Higgs_pT->Write();
      hcos0star->Write();
      hcos01->Write();
      hcos02->Write();
      hphi->Write();
      hphi1->Write();
      h_Higgs_mela1->Write();
      h_Higgs_mela2->Write();
      h_Higgs_mela3->Write();
      h_Higgs_mela4->Write();
      tree->Write();
      f.Close();
      cout<<"File root created..."<<endl;
              
      cout << "____________________________________________________________________________" << endl;
      cout << "! Processo finalizado !" << endl;
      cout << "============================================================================\n" << endl;
      
}