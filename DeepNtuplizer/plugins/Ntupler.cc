// -*- C++ -*-
//
// Package:    NtuPler/Ntupler
// Class:      Ntupler
// 
/**\class Ntupler Ntupler.cc Ntupler/Ntupler/plugins/Ntupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emil Sorensen Bols
//         Created:  Tue, 16 Jan 2018 10:15:06 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TRandom3.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//ROOT includes
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include <TRandom.h>
#include <TH2F.h>
#include <TH1F.h>


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetfwd.h"


#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include "../interface/ntuple_bTagVars1.h"
#include "../interface/ntuple_bTagVars2.h"


#include <algorithm>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Ntupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Ntupler(const edm::ParameterSet&);
      ~Ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::map<std::string,TH2F*> histContainer_;
      std::map<std::string,TH1F*> histoContainer_;
      edm::EDGetTokenT<reco::PFJetCollection>      jetToken1_;
      edm::EDGetTokenT<reco::PFJetCollection>      jetToken2_;
      edm::EDGetTokenT<reco::CaloJetCollection>      jetToken3_;
      edm::EDGetTokenT< std::vector<reco::ShallowTagInfo> > OnShallowSrc_;
      edm::EDGetTokenT< std::vector<reco::ShallowTagInfo> > OffShallowSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OffCSVtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnCSVtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnCSVCalotagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OffBtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnBtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OffCtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnCtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OffUDSGtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnUDSGtagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnBCalotagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnCCalotagSrc_;
      edm::EDGetTokenT<reco::JetTagCollection> OnUDSGCalotagSrc_;
      edm::InputTag onJetsrc_;
      edm::InputTag offJetsrc_;
      edm::Service<TFileService> fs;
      TTree *tree_;
      ntuple_bTagVars1 go;
      ntuple_bTagVars2 go1;
      // ----------member data ---------------------------
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
Ntupler::Ntupler(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   jetToken1_=(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("offJetsrc")));
   jetToken2_=(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("onJetsrc")));
   jetToken3_=(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("onJetCalosrc")));
   OnShallowSrc_=(consumes< std::vector<reco::ShallowTagInfo> >(iConfig.getParameter<edm::InputTag>("onShallowsrc")));
   OffShallowSrc_=(consumes< std::vector<reco::ShallowTagInfo> >(iConfig.getParameter<edm::InputTag>("offShallowsrc")));
   OffBtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("offBtagsrc")));
   OnBtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onBtagsrc")));
   OffCtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("offCtagsrc")));
   OnCtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onCtagsrc")));
   OffUDSGtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("offUDSGtagsrc")));
   OnUDSGtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onUDSGtagsrc")));
   OnCSVtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onCSVtagsrc")));
   OnCSVCalotagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onCSVCalotagsrc")));
   OffCSVtagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("offCSVtagsrc")));
   OnUDSGCalotagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onUDSGCalotagsrc")));
   OnBCalotagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onBCalotagsrc")));
   OnCCalotagSrc_=(consumes< reco::JetTagCollection >(iConfig.getParameter<edm::InputTag>("onCCalotagsrc")));
}


Ntupler::~Ntupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle< std::vector<reco::ShallowTagInfo> > OntagInfos;
   edm::Handle< std::vector<reco::ShallowTagInfo> > OfftagInfos;
   edm::Handle<reco::PFJetCollection> pfjetOnline;
   edm::Handle<reco::CaloJetCollection> calojetOnline;
   edm::Handle<reco::PFJetCollection> pfjetOffline;
   edm::Handle<reco::JetTagCollection> onbtagdisc;
   edm::Handle<reco::JetTagCollection> onbcalotagdisc;
   edm::Handle<reco::JetTagCollection> offbtagdisc;
   edm::Handle<reco::JetTagCollection> oncsvtagdisc;
   edm::Handle<reco::JetTagCollection> oncsvcalotagdisc;
   edm::Handle<reco::JetTagCollection> offcsvtagdisc;
   edm::Handle<reco::JetTagCollection> onctagdisc;
   edm::Handle<reco::JetTagCollection> onccalotagdisc;
   edm::Handle<reco::JetTagCollection> offctagdisc;
   edm::Handle<reco::JetTagCollection> onudsgtagdisc;
   edm::Handle<reco::JetTagCollection> onudsgcalotagdisc;
   edm::Handle<reco::JetTagCollection> offudsgtagdisc;

   iEvent.getByToken(jetToken1_, pfjetOffline);
   iEvent.getByToken(jetToken2_, pfjetOnline);
   iEvent.getByToken(jetToken3_, calojetOnline);
   iEvent.getByToken(OnBtagSrc_, onbtagdisc);
   iEvent.getByToken(OffBtagSrc_, offbtagdisc);
   iEvent.getByToken(OnCSVtagSrc_, oncsvtagdisc);
   iEvent.getByToken(OnCSVCalotagSrc_, oncsvcalotagdisc);
   iEvent.getByToken(OffCSVtagSrc_, offcsvtagdisc);
   iEvent.getByToken(OnCtagSrc_, onctagdisc);
   iEvent.getByToken(OffCtagSrc_, offctagdisc);
   iEvent.getByToken(OnUDSGtagSrc_, onudsgtagdisc);
   iEvent.getByToken(OffUDSGtagSrc_, offudsgtagdisc);
   iEvent.getByToken(OnShallowSrc_, OntagInfos);
   iEvent.getByToken(OffShallowSrc_, OfftagInfos);
   iEvent.getByToken(OnBCalotagSrc_, onbcalotagdisc);
   iEvent.getByToken(OnCCalotagSrc_, onccalotagdisc);
   iEvent.getByToken(OnUDSGCalotagSrc_, onudsgcalotagdisc);
   std::vector<int> matches;
   std::vector<int> matches1;
   
   if(OntagInfos.isValid() & OfftagInfos.isValid() & offbtagdisc.isValid() & onbtagdisc.isValid() ){
     for(size_t n = 0; n < onbtagdisc->size() ; n++){
       double dRMin = 99.0;
       int index = -1;
       for( size_t z = 0; z < offbtagdisc->size() ; z++){
	 double deta = (*offbtagdisc)[z].first->eta() - (*onbtagdisc)[n].first->eta();
	 double dphi = abs( (*offbtagdisc)[z].first->phi() - (*onbtagdisc)[n].first->phi());
	 if(dphi > 3.1415){
	   dphi = 2*3.1415-dphi;
	 }
	 double dR = sqrt(deta*deta + dphi*dphi);
	 if(dR < dRMin){
	   dRMin = dR;
	   index = z;
	 }
	 if(dRMin > 0.3){
	     index = -1;
	 }
       }
       matches.push_back(index);	 
     }
   }
 
   if(OntagInfos.isValid() & OfftagInfos.isValid() & offbtagdisc.isValid() & oncsvtagdisc.isValid() & offcsvtagdisc.isValid() & oncsvcalotagdisc.isValid() & onbtagdisc.isValid() & offctagdisc.isValid() & onctagdisc.isValid() & offudsgtagdisc.isValid() & onudsgtagdisc.isValid() & onudsgcalotagdisc.isValid() & onbcalotagdisc.isValid() & onccalotagdisc.isValid()){
     for(size_t p = 0; p < OntagInfos->size(); p++){
       if(matches.at(p) == -1){continue;}
       const reco::ShallowTagInfo OntagInfo = OntagInfos->at(p);
       const reco::ShallowTagInfo OfftagInfo = OfftagInfos->at(matches.at(p));
       std::vector<float> OnlineDisc;
       std::vector<float> OfflineDisc;
       OnlineDisc.push_back((*onbtagdisc)[p].second);
       OfflineDisc.push_back((*offbtagdisc)[matches.at(p)].second);
       OnlineDisc.push_back((*onctagdisc)[p].second);
       OfflineDisc.push_back((*offctagdisc)[matches.at(p)].second);
       OnlineDisc.push_back((*onudsgtagdisc)[p].second);
       OfflineDisc.push_back((*offudsgtagdisc)[matches.at(p)].second);
       OnlineDisc.push_back((*oncsvtagdisc)[p].second);
       OfflineDisc.push_back((*offcsvtagdisc)[matches.at(p)].second);
       if( onbcalotagdisc->size() == onudsgcalotagdisc->size() == onccalotagdisc->size() == onbtagdisc->size() ){
	 OnlineDisc.push_back((*onbcalotagdisc)[p].second);
	 OnlineDisc.push_back((*onccalotagdisc)[p].second);
	 OnlineDisc.push_back((*onudsgcalotagdisc)[p].second);
	 OnlineDisc.push_back((*oncsvcalotagdisc)[p].second);
       }
       else{
	 OnlineDisc.push_back(-1.0);
	 OnlineDisc.push_back(-1.0);
	 OnlineDisc.push_back(-1.0);
	 OnlineDisc.push_back(-1.0);       
       }
       bool writeit = true;
       if(!go.fillBranches(OntagInfo,iEvent,OnlineDisc)){writeit = false;}
       if(!go1.fillBranches(OfftagInfo,iEvent,OfflineDisc)){writeit = false;}
       if(writeit){tree_->Fill();}
       histContainer_["btag"] ->Fill((*offbtagdisc)[matches.at(p)].second,(*onbtagdisc)[p].second);
       histContainer_["pt"] ->Fill((*offbtagdisc)[matches.at(p)].first->pt(),(*onbtagdisc)[p].first->pt());
     }
   }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntupler::beginJob()
{
  histContainer_["btag"]=fs->make<TH2F>("btag", "online vs offline btag",    90, 0, 1.0, 90, 0, 1.0);
  histContainer_["pt"]=fs->make<TH2F>("pt", "online vs offline pt",    100, 0, 1000.0, 100, 0, 1000.0);
  histoContainer_["AllJetpt"]=fs->make<TH1F>("AllJetpt", "PFJetCollection pt",    100, 0, 500.0);
  histoContainer_["ShallowTagpt"]=fs->make<TH1F>("ShallowTagpt", "ShallowTagpt",    100, 0, 500.0);
  tree_=(fs->make<TTree>("tree" ,"tree" ));
  go.initBranches(tree_);
  go1.initBranches(tree_);
  //tree_->Branch("trackJetPt", &trackJetPt_ , "trackJetPt_/F");

  //addBranch(tree_,"trackJetPt"             , &trackJetPt_             , "trackJetPt_/F"             );
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntupler::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntupler);
