// -*- C++ -*-
//
// Package:    /DeepNTuplizer
// Class:      DeepNTuplizer
// 
/**\class ttsemilepFilter ttsemilepFilter.cc DeepNTuplizer/plugins/ttsemilepFilter.cc

 Description: This filter only takes the events which fulfill the following criteria:
    - exactly one muon in src_muons
    - no electron in src_electrons
    - exactly four jets in src_jets
    - transverse mass of muon and met > minMT_muonMETpair

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Walter
//         Created:  Thu, 12 Apr 2018 17:00:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"


//
// class declaration
//

class ttsemilepFilter : public edm::stream::EDFilter<> {
   public:
      explicit ttsemilepFilter(const edm::ParameterSet&);
      ~ttsemilepFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      edm::EDGetTokenT<edm::View<pat::Muon> >    src_muons;
      edm::EDGetTokenT<edm::View<pat::Electron> >    src_electrons;
      edm::EDGetTokenT<edm::View<pat::Jet> >    src_jets;
      edm::EDGetTokenT<edm::View<pat::MET> >    src_mets;


      edm::Handle<edm::View<pat::Muon> >    muons;
      edm::Handle<edm::View<pat::Electron> >    electrons;
      edm::Handle<edm::View<pat::Jet> >    jets;
      edm::Handle<edm::View<pat::MET> >    mets;

      double minMT_muonMETpair = 0;
      unsigned int nElectrons = 0;
      unsigned int nMuons = 0;

      /*
      double min_chi2 = 0;
      double m_w = 80.4;
      double m_t = 173.1;

      int first_q = 0;
      int second_q = 0;
      int first_b = 0;
      int second_b = 0;

      const std::string btag = "pfCSV_v2JetTags:probb";
      */

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
ttsemilepFilter::ttsemilepFilter(const edm::ParameterSet& iConfig):
    src_muons(consumes<edm::View<pat::Muon> >    (iConfig.getParameter<edm::InputTag>("src_muons"))),
    src_electrons(consumes<edm::View<pat::Electron> >    (iConfig.getParameter<edm::InputTag>("src_electrons"))),
    src_jets(consumes<edm::View<pat::Jet> >    (iConfig.getParameter<edm::InputTag>("src_jets"))),
    src_mets(consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("src_mets"))){
   //now do what ever initialization is needed

    //minPt_muon=(iConfig.getParameter<double>("minPt_muon"));
    //minPt_jet=(iConfig.getParameter<double>("minPt_jet"));
    minMT_muonMETpair=(iConfig.getParameter<double>("cut_minMT_leptonMETpair"));
    nElectrons=(iConfig.getParameter<unsigned int>("nElectrons"));
    nMuons=(iConfig.getParameter<unsigned int>("nMuons"));

}


ttsemilepFilter::~ttsemilepFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ttsemilepFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    iEvent.getByToken(src_muons,  muons);
    iEvent.getByToken(src_electrons,  electrons);
    iEvent.getByToken(src_jets, jets);
    iEvent.getByToken(src_mets, mets);


    if(muons->size()!=nMuons)
        return false;

    if(electrons->size()!=nElectrons)
        return false;

    if(jets->size()!=4)
        return false;

    if(mets->size()==0){ //this should normally not happen, is the file corrupted??
        std::cout<<"there was an event with no met object, strange ..."<<std::endl;
        return false;

    }

    const auto met = mets->ptrAt(0);
    double mWt = -1.;

    if(nMuons==1 && nElectrons==0){
        const auto muon = muons->ptrAt(0);
        mWt = std::sqrt(2*(muon->pt()*met->pt() - muon->px()*met->px() - muon->py()*met->py()));
    }
    else if (nElectrons==1 && nMuons==0){
        const auto electron = electrons->ptrAt(0);
        mWt = std::sqrt(2*(electron->pt()*met->pt() - electron->px()*met->px() - electron->py()*met->py()));
    }




    if(mWt == -1.){                     //no W reconstruction
    }
    else if(mWt < minMT_muonMETpair){   //leptonic decaying W reconstruction
        return false;
    }



    //hadronically decaying t reconstruction
    /*
    for(size_t i = 0; i < jets->size(); ++i){
        for(size_t j = i+1; j < jets->size(); ++j){
            if(i!=j){
                double m_2jet = (jets->ptrAt(i)->p4() + jets->ptrAt(j)->p4()).M();

                for(size_t k = 0; k < jets->size(); ++k){
                    if(i!= j && i!=k && j!=k){
                        double m_3jet = (jets->ptrAt(i)->p4() + jets->ptrAt(j)->p4() + jets->ptrAt(k)->p4()).M();

                        double chi2 = - std::log(std::pow(m_w - m_2jet,2)/m_2jet + std::pow(m_t - m_3jet,2)/m_3jet);
                        std::cout<<"combination ("<<i<<j<<k<<") has m_2jet = "<< m_2jet<<" , m_3jet = "<<m_3jet<<" chi2 = "<<chi2<<std::endl;

                        if(chi2 > min_chi2){
                            min_chi2 = chi2;
                            first_q = i;
                            second_q = j;
                            first_b = k;
                            second_b = 6 - i - j - k;
                        }
                    }
                }
            }
        }
    }

    if(jets->ptrAt(first_b)->bDiscriminator()
    */

    return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ttsemilepFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ttsemilepFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ttsemilepFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ttsemilepFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ttsemilepFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ttsemilepFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ttsemilepFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ttsemilepFilter);
