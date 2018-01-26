
// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer2
// Class:      DemoAnalyzer2
// 
/**\class DemoAnalyzer2 DemoAnalyzer2.cc Demo/DemoAnalyzer2/plugins/DemoAnalyzer2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Miqueias Melo De Almeida
//         Created:  Tue, 30 May 2017 19:01:58 GMT
//
//


#include <math.h>
#include <vector>
#include <TFile.h>
#include <TTree.h>

//vectors with the events of interest
//#include "VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_jets.h"
#include "GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8_jets.h"
//#include "ZZTo4L_13TeV_powheg_pythia8_jets.h"



// system include files
#include <memory>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


bool decresing (int i, int j) { return (i>j); }

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DemoAnalyzer2 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer2(const edm::ParameterSet&);
      ~DemoAnalyzer2();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<pat::Jet>> PFJets_;
  TFile *ofile;
  TTree *tree;
  Int_t event;
  Int_t lumi;
  Int_t run;
  Float_t jet_eta[3], jet_phi[3], jet_subjetness[3], jet_ptD[3], jet_closure[3];
  Float_t jet_photonEnergy[3], jet_electronEnergy[3], jet_muonEnergy[3], jet_chargedEmEnergy[3], jet_neutralEmEnergy[3], jet_chargedHadronEnergy[3], jet_neutralHadronEnergy[3];
  Float_t subjet_pt[3][100], subjet_eta[3][100], subjet_phi[3][100];  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer2::DemoAnalyzer2(const edm::ParameterSet& iConfig):
  PFJets_ ( consumes<std::vector<pat::Jet>> (iConfig.getParameter<edm::InputTag> ("PFJets")) )
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   tree = new TTree("HZZ4LeptonsAnalysisReduced","Properties of jets in selected events HZZ4L analysis");
   tree->Branch("run",&run,"run/I");
   tree->Branch("lumi",&lumi,"lumi/I");
   tree->Branch("event",&event,"event/I");
   tree->Branch("jet_eta",&jet_eta,"jet_eta[3]/F");
   tree->Branch("jet_phi",&jet_phi,"jet_phi[3]/F");
   tree->Branch("jet_closure",&jet_closure,"jet_closure[3]/F");
   tree->Branch("jet_subjetness",&jet_subjetness,"jet_subjetness[3]/F");
   tree->Branch("jet_ptD",&jet_ptD,"jet_ptD[3]/F");
   tree->Branch("jet_photonEnergy",&jet_photonEnergy,"jet_photonEnergy[3]/F");
   tree->Branch("jet_electronEnergy",&jet_electronEnergy,"jet_electronEnergy[3]/F");
   tree->Branch("jet_muonEnergy",&jet_muonEnergy,"jet_muonEnergy[3]/F");
   tree->Branch("jet_chargedEmEnergy",&jet_chargedEmEnergy,"jet_chargedEmEnergy[3]/F");
   tree->Branch("jet_neutralEmEnergy",&jet_neutralEmEnergy,"jet_neutralEmEnergy[3]/F");
   tree->Branch("jet_chargedHadronEnergy",&jet_chargedHadronEnergy,"jet_chargedHadronEnergy[3]/F");
   tree->Branch("jet_neutralHadronEnergy",&jet_neutralHadronEnergy,"jet_neutralHadronEnergy[3]/F");
   tree->Branch("jet_component_pt",&subjet_pt,"jet_component_pt[3][100]/F");
   tree->Branch("jet_component_eta",&subjet_eta,"jet_component_eta[3][100]/F");
   tree->Branch("jet_component_phi",&subjet_phi,"jet_component_phi[3][100]/F");
}


DemoAnalyzer2::~DemoAnalyzer2()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
    //ofile = new TFile("VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_JetProperties.root","recreate");
    ofile = new TFile("GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8_JetProperties.root","recreate");
    //ofile = new TFile("ZZTo4L_13TeV_powheg_pythia8_JetProperties.root","recreate");

    tree->Write();
    ofile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif


   //------------------ My analyzer itself ------------------

   //Check if event is the one looked for                                                                                                                   
   // Event info                                                                                                                                            
   run = iEvent.id().run();
   lumi = iEvent.id().luminosityBlock();
   event = iEvent.id().event();
   
   int Registry = -1;
   //sel_events is a vector in header files containing the selected events
   for(unsigned int iSelEv=0; iSelEv<(unsigned int)sel_events.size(); ++iSelEv)
     if(sel_events.at(iSelEv).at(0) == run && sel_events.at(iSelEv).at(1) == lumi && sel_events.at(iSelEv).at(2) == event)
       Registry = iSelEv;

   
   if(Registry != -1){
     std::cout<< Form("Found registry --> Run: %i, Lumi: %i, Event: %i",run,lumi,event) <<std::endl;
     //The sel_events vector conatins the jets eta,phi
     float selJet1Eta = sel_events.at(Registry).at(3);
     float selJet1Phi = sel_events.at(Registry).at(4);
     float selJet2Eta = sel_events.at(Registry).at(5);
     float selJet2Phi = sel_events.at(Registry).at(6);
     float selJet3Eta = sel_events.at(Registry).at(7);
     float selJet3Phi = sel_events.at(Registry).at(8);
     
     for(unsigned int iel=0; iel<3; ++iel){
      jet_eta[iel] = -999;
      jet_phi[iel] = -999;
      jet_subjetness[iel] = -999;
      jet_ptD[iel] = -999;
      jet_closure[iel] = -999;
      jet_photonEnergy[iel] = -999;
      jet_electronEnergy[iel] = -999;
      jet_muonEnergy[iel] = -999;
      jet_chargedEmEnergy[iel] = -999;
      jet_neutralEmEnergy[iel] = -999;
      jet_chargedHadronEnergy[iel] = -999;
      jet_neutralHadronEnergy[iel] = -999;
      
      for(unsigned int isj=0; isj<100; ++isj){
	subjet_pt[iel][isj] = -999;
	subjet_eta[iel][isj] = -999;
	subjet_phi[iel][isj] = -999;
      }
     }

     edm::Handle<std::vector<pat::Jet>> pfjetH;
     iEvent.getByToken(PFJets_, pfjetH);

     //Accessing the PFJets to find the closest jets
     float minDR1 = 99, minDR2 = 99, minDR3 = 99;
     std::vector<float> spt1, seta1, sphi1;
     std::vector<float> spt2, seta2, sphi2;
     std::vector<float> spt3, seta3, sphi3;

     for ( std::vector<pat::Jet>::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ){
       float pfJetEta = jet->eta();
       float pfJetPhi = jet->phi();
       unsigned int Nsubjets = jet->numberOfDaughters();
       
       float DR1 = std::sqrt( std::pow(selJet1Eta-pfJetEta,2) + std::pow(selJet1Phi-pfJetPhi,2) );
       if(DR1 < minDR1 && DR1 >= 0){
	 minDR1 = DR1;
	 jet_closure[0] = DR1;
	 jet_eta[0] = pfJetEta;
	 jet_phi[0] = pfJetPhi;
         jet_subjetness[0] = Nsubjets;
         jet_photonEnergy[0] = jet->photonEnergy();
         jet_electronEnergy[0] = jet->electronEnergy();
         jet_muonEnergy[0] = jet->muonEnergy();
         jet_chargedEmEnergy[0] = jet->chargedEmEnergy();
         jet_neutralEmEnergy[0] = jet->neutralEmEnergy();
         jet_chargedHadronEnergy[0] = jet->chargedHadronEnergy();
         jet_neutralHadronEnergy[0] = jet->neutralHadronEnergy();
	 
	 spt1.clear();
	 seta1.clear();
	 sphi1.clear();
	 float sumpt = 0, sumpt2 = 0;
 	 for(unsigned int idau=0; idau<Nsubjets; ++idau){
	   const reco::Candidate *pfCand = jet->daughter(idau);
	   sumpt += pfCand->pt();
	   sumpt2 += pfCand->pt()*pfCand->pt();
	   spt1.push_back( pfCand->pt() );
	   seta1.push_back( pfCand->eta() );
	   sphi1.push_back( pfCand->phi() );
	 }//End loop over PF subjets
	 jet_ptD[0] = sqrt(sumpt2)/sumpt;
       }
	 
       float DR2 = std::sqrt( std::pow(selJet2Eta-pfJetEta,2) + std::pow(selJet2Phi-pfJetPhi,2) );
       if(DR2 < minDR2 && DR2 >= 0){
	 minDR2 = DR2;
	 jet_closure[1] = DR2;
	 jet_eta[1] = pfJetEta;
	 jet_phi[1] = pfJetPhi;
         jet_subjetness[1] = Nsubjets;
         jet_photonEnergy[1] = jet->photonEnergy();
         jet_electronEnergy[1] = jet->electronEnergy();
         jet_muonEnergy[1] = jet->muonEnergy();
         jet_chargedEmEnergy[1] = jet->chargedEmEnergy();
         jet_neutralEmEnergy[1] = jet->neutralEmEnergy();
         jet_chargedHadronEnergy[1] = jet->chargedHadronEnergy();
         jet_neutralHadronEnergy[1] = jet->neutralHadronEnergy();
	 
	 spt2.clear();
	 seta2.clear();
	 sphi2.clear();
	 float sumpt = 0, sumpt2 = 0;
 	 for(unsigned int idau=0; idau<Nsubjets; ++idau){
           const reco::Candidate *pfCand = jet->daughter(idau);
	   sumpt += pfCand->pt();
	   sumpt2 += pfCand->pt()*pfCand->pt();
	   spt2.push_back( pfCand->pt() );
	   seta2.push_back( pfCand->eta() );
	   sphi2.push_back( pfCand->phi() );
	 }//End loop over PF subjets
	 jet_ptD[1] = sqrt(sumpt2)/sumpt;
       }

       if(selJet3Eta != -999){
	 float DR3 = std::sqrt( std::pow(selJet3Eta-pfJetEta,2) + std::pow(selJet3Phi-pfJetPhi,2) );
         if(DR3 < minDR3 && DR1 >= 0){
	   minDR3 = DR3;
	   jet_closure[2] = DR3;
	   jet_eta[2] = pfJetEta;
	   jet_phi[2] = pfJetPhi;
           jet_subjetness[2] = Nsubjets;
           jet_photonEnergy[2] = jet->photonEnergy();
           jet_electronEnergy[2] = jet->electronEnergy();
           jet_muonEnergy[2] = jet->muonEnergy();
           jet_chargedEmEnergy[2] = jet->chargedEmEnergy();
           jet_neutralEmEnergy[2] = jet->neutralEmEnergy();
           jet_chargedHadronEnergy[2] = jet->chargedHadronEnergy();
           jet_neutralHadronEnergy[2] = jet->neutralHadronEnergy();
	 
	   spt3.clear();
	   seta3.clear();
	   sphi3.clear();
	   float sumpt = 0, sumpt2 = 0;
 	   for(unsigned int idau=0; idau<Nsubjets; ++idau){
             const reco::Candidate *pfCand = jet->daughter(idau);
	     sumpt += pfCand->pt();
	     sumpt2 += pfCand->pt()*pfCand->pt();
	     spt3.push_back( pfCand->pt() );
	     seta3.push_back( pfCand->eta() );
	     sphi3.push_back( pfCand->phi() );
	   }//End loop over PF subjets
	   jet_ptD[2] = sqrt(sumpt2)/sumpt;
         }
       }

     }//End PFJet loop
     
     
     //Sort jet components by pT
     unsigned int ncomponents1 = spt1.size();
     if(ncomponents1 > 0){
       int ijet1 = 0;
       for(unsigned int i1=0; i1 < ncomponents1; ++i1){
         float maxpt = -1, ci2 = -1;
         for(unsigned int i2=0; i2 < ncomponents1; ++i2){
	   if(spt1[i2] > maxpt){
	     maxpt = spt1[i2];
	     ci2 = i2;
	   }
	 }
	 
         subjet_pt[0][ijet1] = spt1[ci2];
         subjet_eta[0][ijet1] = seta1[ci2];
         subjet_phi[0][ijet1] = sphi1[ci2]; 

         spt1[ci2] = -2;
         seta1[ci2] = -2;
         sphi1[ci2] = -2;
         ++ijet1;
       }
     }
       
     unsigned int ncomponents2 = spt2.size();
     if(ncomponents2 > 0){
       int ijet2 = 0;
       for(unsigned int i1=0; i1 < ncomponents2; ++i1){
         float maxpt = -1, ci2 = -1;
         for(unsigned int i2=0; i2 < ncomponents2; ++i2){
	   if(spt2[i2] > maxpt){
	     maxpt = spt2[i2];
	     ci2 = i2;
	   }
	 }
      
         subjet_pt[1][ijet2] = spt2[ci2];
         subjet_eta[1][ijet2] = seta2[ci2];
         subjet_phi[1][ijet2] = sphi2[ci2]; 
       
         spt2[ci2] = -2;
         seta2[ci2] = -2;
         sphi2[ci2] = -2;
         ++ijet2;
       }
     }

     unsigned int ncomponents3 = spt3.size();
     if(ncomponents3 > 0){
       int ijet3 = 0;
       for(unsigned int i1=0; i1 < ncomponents3; ++i1){
	 float maxpt = -1, ci2 = -1;
	 for(unsigned int i2=0; i2 < ncomponents3; ++i2){
	   if(spt3[i2] > maxpt){
	     maxpt = spt3[i2];
	     ci2 = i2;
	   }
	 }
	 
	 subjet_pt[2][ijet3] = spt3[ci2];
	 subjet_eta[2][ijet3] = seta3[ci2];
	 subjet_phi[2][ijet3] = sphi3[ci2]; 

         spt3[ci2] = -2;
         seta3[ci2] = -2;
         sphi3[ci2] = -2;
	 ++ijet3;
       }
     }
     
     //Fills the tree
     tree->Fill();
   }//Event filter

}//End of the analyzer member function


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer2::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer2::endJob() 
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer2);
