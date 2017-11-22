// -*- C++ -*-
//
// Package:    DeepNTuples/DeepNtuplizer
// Class:      puAnalyzer
// 
/**\class puAnalyzer puAnalyzer.cc DeepNTuples/DeepNtuplizer/plugins/puAnalyzer.cc

 Description: This analyzer is used in pupMCHistProd.py to calculate the pileup from the number of interactions per bunch crossing

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Walter
//         Created:  Mon, 20 Nov 2017 18:49:19 GMT
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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class puAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      puAnalyzer(const edm::ParameterSet&);

      ~puAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT< std::vector<PileupSummaryInfo> > EDMPUInfoToken;
      TH1F * h_pileup;
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
puAnalyzer::puAnalyzer(const edm::ParameterSet& iConfig):
   EDMPUInfoToken(consumes<std::vector<PileupSummaryInfo> > (iConfig.getParameter<edm::InputTag>("pileupInfo"))){

   edm::Service<TFileService> fs;
   h_pileup = fs->make<TH1F>( "pileup"  , "Number of interactions per bunch crossing", 60,  0., 60 );

}


puAnalyzer::~puAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
puAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //using namespace edm;

   edm::Handle< std::vector<PileupSummaryInfo> >  PupInfo;

   iEvent.getByToken(EDMPUInfoToken, PupInfo);


   std::vector<PileupSummaryInfo>::const_iterator PVI;
   float Tnpv = -1;
   for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

       int BX = PVI->getBunchCrossing();

       if(BX == 0) {
         Tnpv = PVI->getTrueNumInteractions();
         continue;
       }
   }
   h_pileup->Fill(Tnpv);

}


// ------------ method called once each job just before starting event loop  ------------
void 
puAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
puAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
puAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(puAnalyzer);
