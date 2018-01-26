// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

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
//#include "VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.h"
//#include "GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8.h"
#include "ZZTo4L_13TeV_powheg_pythia8.h"



// system include files
#include <memory>

// user include files
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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<reco::PFJetCollection> PFJets_;
  TFile *ofile;
  TTree *tree;
  int event;
  int lumi;
  int run;
  std::vector<float> jet_eta, jet_phi, subjetness, ptD, jet_closure;
  std::vector<std::vector<float>> subjet_pt, subjet_eta, subjet_phi;  
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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):
  PFJets_ ( consumes<reco::PFJetCollection> (iConfig.getParameter<edm::InputTag>("PFJets")) )
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   tree = new TTree("HZZ4LeptonsAnalysisReduced","Properties of jets in selected events HZZ4L analysis");
   tree->Branch("run",&run,"run/I");
   tree->Branch("lumi",&lumi,"lumi/I");
   tree->Branch("event",&event,"event/I");
   tree->Branch("jet_eta","std::vector<float>",&jet_eta);
   tree->Branch("jet_phi","std::vector<float>",&jet_phi);
   tree->Branch("jet_closure","std::vector<float>",&jet_closure);
   tree->Branch("subjetness","std::vector<float>",&subjetness);
   tree->Branch("ptD","std::vector<float>",&ptD);
   tree->Branch("subjet_pt","std::vector<std::vector<float>>",&subjet_pt);
   tree->Branch("subjet_eta","std::vector<std::vector<float>>",&subjet_eta);
   tree->Branch("subjet_phi","std::vector<std::vector<float>>",&subjet_phi);
   
}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
    ofile = new TFile("JetProperties.root","recreate");
    tree->Write();
    ofile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   lumi = iEvent.luminosityBlock();
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
     
     
     jet_eta.clear();
     jet_phi.clear();
     subjetness.clear();
     ptD.clear();
     jet_closure.clear();
     subjet_pt.clear();
     subjet_eta.clear();
     subjet_phi.clear();  

     if(selJet3Eta == -999){
       jet_eta = {-999, -999};
       jet_phi = {-999, -999};
       subjetness = {-999, -999};
       ptD = {-999, -999};
       jet_closure = {-999, -999};
     }
     else{
       jet_eta = {-999, -999, -999};
       jet_phi = {-999, -999, -999};
       subjetness = {-999, -999, -999};
       ptD = {-999, -999, -999};
       jet_closure = {-999, -999, -999};
     }
     
     edm::Handle<reco::PFJetCollection> pfjetH;
     iEvent.getByToken(PFJets_, pfjetH);
     //unsigned int Njets = pfjetH->size();
     //std::cout<< "---------------- Found Registry --------------------" <<std::endl;
     //std::cout<< Form("Run: %i, Lumi: %i, Event: %i, Njets: %i",Run,Lumi,Event,Njets) <<std::endl;

     //Accessing the PFJets to find the closest jets
     float minDR1 = 99, minDR2 = 99, minDR3 = 99;
     std::vector<float> spt1, seta1, sphi1;
     std::vector<float> spt2, seta2, sphi2;
     std::vector<float> spt3, seta3, sphi3;
     for ( reco::PFJetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ){
       float pfJetEta = jet->eta();
       float pfJetPhi = jet->phi();
       int Nsubjets = jet->numberOfDaughters();
       
       float DR1 = std::sqrt( std::pow(selJet1Eta-pfJetEta,2) + std::pow(selJet1Phi-pfJetPhi,2) );
       if(DR1 < minDR1 && DR1 >= 0){
	 minDR1 = DR1;
	 jet_closure[0] = DR1;
	 jet_eta[0] = pfJetEta;
	 jet_phi[0] = pfJetPhi;
	 
	 spt1.clear();
	 seta1.clear();
	 sphi1.clear();
	 float sumpt = 0, sumpt2 = 0;
 	 for(int idau=0; idau<Nsubjets; ++idau){
	   reco::Candidate const * pfCand = jet->daughter(idau);
	   sumpt += pfCand->pt();
	   sumpt2 += pfCand->pt()*pfCand->pt();
	   spt1.push_back( pfCand->pt() );
	   seta1.push_back( pfCand->eta() );
	   sphi1.push_back( pfCand->phi() );
	 }//End loop over PF subjets
	 ptD[0] = sqrt(sumpt2)/sumpt; 
	 subjetness[0] = Nsubjets;
       }
	 
       float DR2 = std::sqrt( std::pow(selJet2Eta-pfJetEta,2) + std::pow(selJet2Phi-pfJetPhi,2) );
       if(DR2 < minDR2 && DR2 >= 0){
	 minDR2 = DR2;
	 jet_closure[1] = DR2;
	 jet_eta[1] = pfJetEta;
	 jet_phi[1] = pfJetPhi;
	 
	 spt2.clear();
	 seta2.clear();
	 sphi2.clear();
	 float sumpt = 0, sumpt2 = 0;
 	 for(int idau=0; idau<Nsubjets; ++idau){
	   reco::Candidate const * pfCand = jet->daughter(idau);
	   sumpt += pfCand->pt();
	   sumpt2 += pfCand->pt()*pfCand->pt();
	   spt2.push_back( pfCand->pt() );
	   seta2.push_back( pfCand->eta() );
	   sphi2.push_back( pfCand->phi() );
	 }//End loop over PF subjets
	 ptD[1] = sqrt(sumpt2)/sumpt; 
	 subjetness[1] = Nsubjets;
       }

       if(selJet3Eta != -999){
	 float DR3 = std::sqrt( std::pow(selJet3Eta-pfJetEta,2) + std::pow(selJet3Phi-pfJetPhi,2) );
         if(DR3 < minDR3 && DR1 >= 0){
	   minDR3 = DR3;
	   jet_closure[2] = DR3;
	   jet_eta[2] = pfJetEta;
	   jet_phi[2] = pfJetPhi;
	 
	   spt3.clear();
	   seta3.clear();
	   sphi3.clear();
	   float sumpt = 0, sumpt2 = 0;
 	   for(int idau=0; idau<Nsubjets; ++idau){
	     reco::Candidate const * pfCand = jet->daughter(idau);
	     sumpt += pfCand->pt();
	     sumpt2 += pfCand->pt()*pfCand->pt();
	     spt3.push_back( pfCand->pt() );
	     seta3.push_back( pfCand->eta() );
	     sphi3.push_back( pfCand->phi() );
	   }//End loop over PF subjets
	   ptD[2] = sqrt(sumpt2)/sumpt; 
	   subjetness[2] = Nsubjets;
         }
       }
     }
     
     subjet_pt.push_back(spt1);
     subjet_eta.push_back(seta1);
     subjet_phi.push_back(sphi1); 

     subjet_pt.push_back(spt2);
     subjet_eta.push_back(seta2);
     subjet_phi.push_back(sphi2);

     if(selJet3Eta != -999){
       subjet_pt.push_back(spt3);
       subjet_eta.push_back(seta3);
       subjet_phi.push_back(sphi3);
     }
       
       
     //Fills the tree
     tree->Fill();
   }//Event filter

}//End of the analyzer member function


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
