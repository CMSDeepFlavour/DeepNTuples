#include <vector>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class MuonIdAdder : public edm::EDProducer {
    public:
    explicit MuonIdAdder(const edm::ParameterSet& iConfig):
    src_(consumes<edm::View<pat::Muon> >    (iConfig.getParameter<edm::InputTag>("src"))),
    vSrc_(consumes<edm::View<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vSrc"))){
        produces<std::vector<pat::Muon> >();
    }
    ~MuonIdAdder(){}
    void produce(edm::Event & iEvent, const edm::EventSetup& iSetup) override;

    private:
        edm::EDGetTokenT<edm::View<pat::Muon> >    src_;
        edm::EDGetTokenT<edm::View<reco::Vertex> > vSrc_;
};

void MuonIdAdder::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    edm::Handle<edm::View<pat::Muon> >    muons;
    edm::Handle<edm::View<reco::Vertex> > vertex;

    iEvent.getByToken(src_,  muons);
    iEvent.getByToken(vSrc_, vertex);

    std::vector<pat::Muon> * out = new std::vector<pat::Muon>; out->reserve(muons->size());

    for (size_t i = 0; i < muons->size(); ++i){
        const auto muon = muons->ptrAt(i);
        bool idresult = muon->isTightMuon(*(vertex->ptrAt(0)));

        pat::Muon newmuon(*muon);
        newmuon.addUserFloat("tightcutbased",(float)idresult);
        out->push_back(newmuon);
    }
  std::unique_ptr<std::vector<pat::Muon> > ptr(out);
  iEvent.put(std::move(ptr));

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIdAdder);