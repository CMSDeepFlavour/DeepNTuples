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
    explicit MuonIdAdder(const edm::ParameterSet&);
    ~MuonIdAdder(){}
    void produce(edm::Event & iEvent, const edm::EventSetup& iSetup) override;

    private:
        edm::EDGetTokenT<edm::View<pat::Muon> >    src_;
        edm::EDGetTokenT<edm::View<reco::Vertex> > vSrc_;

        edm::Handle<edm::View<pat::Muon> >    muons;
        edm::Handle<edm::View<reco::Vertex> > vertex;

        double minPt_;
        double maxAbsEta_;
        double maxRMI_;
};

MuonIdAdder::MuonIdAdder(const edm::ParameterSet& iConfig):
    src_(consumes<edm::View<pat::Muon> >    (iConfig.getParameter<edm::InputTag>("src"))),
    vSrc_(consumes<edm::View<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vSrc"))){
        produces<std::vector<pat::Muon> >();

        minPt_=(iConfig.getParameter<double>("minPt"));
        maxAbsEta_=(iConfig.getParameter<double>("maxAbsEta"));
        maxRMI_ = (iConfig.getParameter<double>("maxRMI"));
    }

void MuonIdAdder::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {


    iEvent.getByToken(src_,  muons);
    iEvent.getByToken(vSrc_, vertex);

    std::vector<pat::Muon> * out = new std::vector<pat::Muon>;
    out->reserve(muons->size());

    for (size_t i = 0; i < muons->size(); ++i){
        const auto muon = muons->ptrAt(i);
        bool idresult = muon->isTightMuon(*(vertex->ptrAt(0)));
        double rmi = (muon->pfIsolationR04().sumChargedHadronPt +std::max(0.0, (muon->pfIsolationR04().sumNeutralHadronEt +muon->pfIsolationR04().sumPhotonEt -0.5*muon->pfIsolationR04().sumPUPt)))/muon->pt();

        pat::Muon newmuon(*muon);

        if(idresult && rmi < maxRMI_ && muon->pt() > minPt_ && std::abs(muon->eta())<maxAbsEta_){
            out->push_back(newmuon);
        }

    }
    std::unique_ptr<std::vector<pat::Muon> > ptr(out);
    iEvent.put(std::move(ptr));
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIdAdder);