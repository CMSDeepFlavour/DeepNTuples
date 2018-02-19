// -*- C++ -*-
//
// Package:    DeepNTuples/DeepNtuplizer
// Class:      puAnalyzer
//
/**\class lheAnalyzer lheAnalyzer.cc DeepNTuples/DeepNtuplizer/plugins/lheAnalyzer.cc
 Description: This analyzer is used in LHEcounter.py to produce a histogram with the lhe weights of the input files.
 This is needed to calculate the effective event number of mc events which is the total number of events minus the once with negative lhe weights
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Walter
//         Created:  Wed, 22 Nov 2017 09:40:19 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class lheAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      lheAnalyzer(const edm::ParameterSet&);

      ~lheAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<LHEEventProduct> t_LHEInfo;
      TH1F * hist_lheWeight_0;
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
lheAnalyzer::lheAnalyzer(const edm::ParameterSet& iConfig):
   t_LHEInfo(consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>("lheInfo"))){

   edm::Service<TFileService> fs;
   hist_lheWeight_0 = fs->make<TH1F>( "lheweight"  , "lhe weight [0]", 10,  -2., 2. );

}


lheAnalyzer::~lheAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
lheAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   edm::Handle<LHEEventProduct>  h_LHEWeights;

   iEvent.getByToken(t_LHEInfo, h_LHEWeights);

   double theWeight;

   theWeight = h_LHEWeights->weights()[0].wgt/std::abs(h_LHEWeights->weights()[0].wgt);



   hist_lheWeight_0->Fill(theWeight);

}


// ------------ method called once each job just before starting event loop  ------------
void
lheAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
lheAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
lheAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(lheAnalyzer);