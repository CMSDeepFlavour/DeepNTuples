// -*- C++ -*-
//
// Package:    DeepNTuples/DeepNtuplizer
// Class:      globalInfo
//
/**\class globalInfo globalInfo.cc DeepNTuples/DeepNtuplizer/plugins/globalInfo.cc
 Description: This analyzer is used to count the (effective) number of initial events.
 In the LHE case, the events with negative lhe weights have to be substracted from these with positive lhe weights
*/
//
// Original Author:  David Walter
//         Created:  Tue, 20 Feb 2018 17:18:02 GMT
//
//

// system include files
#include <iostream>
#include <memory>
#include <vector>
#include <TList.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


class globalInfo : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        globalInfo(const edm::ParameterSet&);

        ~globalInfo();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::Service<TFileService> fs;
        edm::EDGetTokenT<LHEEventProduct> t_LHEInfo;
        edm::Handle<LHEEventProduct>  h_LHEWeights;

        long nEvents_ = 0;
        long nNegLHEEvents_ = 0;
        bool useLHEWeights_ = true;

        TList * infolist = new TList;


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
globalInfo::globalInfo(const edm::ParameterSet& iConfig):
    t_LHEInfo(consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>("lheInfo")))
    {
    useLHEWeights_ = (iConfig.getParameter<bool>("useLHEWeights"));

}


globalInfo::~globalInfo()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
globalInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    nEvents_++;

    if(useLHEWeights_){
        iEvent.getByToken(t_LHEInfo, h_LHEWeights);
        double lheWeight = h_LHEWeights->weights()[0].wgt/std::abs(h_LHEWeights->weights()[0].wgt);
        if(lheWeight < 0.){
            nNegLHEEvents_++;
        }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
globalInfo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
globalInfo::endJob()
{
    std::cout<<"total number of initial events is: "<<nEvents_<<std::endl;
    std::cout<<"number of initial events with negative lhe weight is: "<<nNegLHEEvents_<<std::endl;
    std::cout<<"effective event number is: "<<nEvents_ - nNegLHEEvents_<<std::endl;
    if( !fs ){
        throw edm::Exception( edm::errors::Configuration,
                "TFile Service is not registered in cfg file" );
    }

    infolist = (fs->make<TList>());

    infolist->Add(new TNamed("nInitialEvents",std::to_string(nEvents_ - nNegLHEEvents_).c_str()));

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
globalInfo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(globalInfo);