#include <vector>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"



class ElectronIdAdder : public edm::EDProducer {
    public:
    explicit ElectronIdAdder(const edm::ParameterSet& iConfig):
    src_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src"))),
    idMap_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("idMap"))),
    vSrc_(consumes<edm::View<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vSrc"))){
        produces<std::vector<pat::Electron> >();
    }
    ~ElectronIdAdder(){}
    void produce(edm::Event & iEvent, const edm::EventSetup& iSetup) override;

    private:
        edm::EDGetTokenT<edm::View<pat::Electron> > src_;
        edm::EDGetTokenT<edm::ValueMap<bool> > idMap_;
        edm::EDGetTokenT<edm::View<reco::Vertex> > vSrc_;
};

void ElectronIdAdder::produce(edm::Event & iEvent, const edm::EventSetup & es) {
    edm::Handle<edm::View<pat::Electron> > electrons;
    edm::Handle<edm::ValueMap<bool> > map;
    edm::Handle<edm::View<reco::Vertex> > vertex;

    iEvent.getByToken(src_, electrons);
    iEvent.getByToken(idMap_, map);
    iEvent.getByToken(vSrc_, vertex);

    std::vector<pat::Electron> * out = new std::vector<pat::Electron>; out->reserve(electrons->size());

    for (size_t i = 0; i < electrons->size(); ++i){
        const auto iElectron = electrons->ptrAt(i);
        bool idresult = (*map)[iElectron];
        bool notInCrack = false;
        double absd0 = 1.;
        double absdz = 1.;
        bool inAbsD0 = false;
        bool inAbsDz = false;
        double global_default_value = 1.5;
        double SCeta = (iElectron->superCluster().isAvailable()) ? iElectron->superCluster()->position().eta() : global_default_value;
        double absSCeta = fabs(SCeta);


        if( iElectron->superCluster().isAvailable() ){
            notInCrack = ( absSCeta<1.4442 || absSCeta>1.5660 );
        }
        if( iElectron->gsfTrack().isAvailable() ){
            absd0 = fabs(iElectron->gsfTrack()->dxy(vertex->ptrAt(0)->position()));
            absdz = fabs(iElectron->gsfTrack()->dz(vertex->ptrAt(0)->position()));
        }
        if( iElectron->superCluster().isAvailable() ){
            if( absSCeta < 1.4442){
                inAbsD0 = absd0 < 0.05;
                inAbsDz = absdz < 0.1;
            }
            if( absSCeta > 1.5660){
                inAbsD0 = absd0 < 0.1;
                inAbsDz = absdz < 0.2;
            }
        }

        pat::Electron newele(*iElectron);
        newele.addUserFloat("tightcutbased",(float)idresult);
        newele.addUserFloat("notInEtaVetoRegion", (float)notInCrack);
        newele.addUserFloat("inAbsD0",(float)inAbsD0);
        newele.addUserFloat("inAbsDz",(float)inAbsDz);
        out->push_back(newele);
    }
  std::unique_ptr<std::vector<pat::Electron> > ptr(out);
  iEvent.put(std::move(ptr));

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronIdAdder);