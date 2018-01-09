// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "../interface/ntuple_content.h"
#include "../interface/ntuple_SV.h"
#include "../interface/ntuple_JetInfo.h"
#include "../interface/ntuple_pfCands.h"
#include "../interface/ntuple_bTagVars.h"
#include "../interface/ntuple_FatJetInfo.h"
#include "../interface/ntuple_DeepVertex.h"

//ROOT includes
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include <TRandom.h>

//CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include <algorithm>

//trash?

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"



#if defined( __GXX_EXPERIMENTAL_CXX0X__)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#endif

struct MagneticField;


class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
    explicit DeepNtuplizer(const edm::ParameterSet&);
    ~DeepNtuplizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const;
    Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const ;
    float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const ;


    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
    edm::EDGetTokenT<edm::View<pat::Jet> >      jetToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
    edm::EDGetTokenT<double> rhoToken_;

    std::string t_qgtagger;

    edm::Service<TFileService> fs;
    TTree *tree_;

    size_t njetstotal_;
    size_t njetswithgenjet_;
    size_t njetsselected_;

	ntuple_content * addModule(ntuple_content *m, std::string name = ""){
        modules_.push_back(m);
				module_names_.push_back(name);
        return m;
    }
    std::vector<ntuple_content* > modules_;
	std::vector<std::string> module_names_;

    bool applySelection_;
};

DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
                            vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
                            svToken_(consumes<reco::VertexCompositePtrCandidateCollection>(
                                    iConfig.getParameter<edm::InputTag>("secVertices"))),
                                    jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
                                    puToken_(consumes<std::vector<PileupSummaryInfo >>(iConfig.getParameter<edm::InputTag>("pupInfo"))),
                                    rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInfo"))),
                                    t_qgtagger(iConfig.getParameter<std::string>("qgtagger"))
{
    /*
     *  Initialise the modules here
     *  Everything else does not need to be changed if
     *  modules don't interact.
     */

    // read configuration parameters
    const double jetR = iConfig.getParameter<double>("jetR");
    const bool  runFatJets_ = iConfig.getParameter<bool>("runFatJet");

    //not implemented yet
    const bool useHerwigCompatibleMatching=iConfig.getParameter<bool>("useHerwigCompatible");
    const bool isHerwig=iConfig.getParameter<bool>("isHerwig");

    ntuple_content::useoffsets = iConfig.getParameter<bool>("useOffsets");

    applySelection_=iConfig.getParameter<bool>("applySelection");

    ntuple_SV* svmodule=new ntuple_SV("", jetR);
    addModule(svmodule, "SVNtuple");

    //Loose IVF vertices
    //ntuple_SV* svmodule_LooseIVF=new ntuple_SV("LooseIVF_", jetR);
    //svmodule_LooseIVF->setSVToken(
    //        consumes<reco::VertexCompositePtrCandidateCollection>(
    //                iConfig.getParameter<edm::InputTag>("LooseSVs")));
    //removed LooseIVF module
    //addModule(svmodule_LooseIVF);

    // DeepVertex info
    // ntuple_DeepVertex* deepvertexmodule=new ntuple_DeepVertex(jetR);
    // deepvertexmodule->setCandidatesToken(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("candidates")));
    // addModule(deepvertexmodule);

    ntuple_JetInfo* jetinfo=new ntuple_JetInfo();
    jetinfo->setQglToken(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "qgLikelihood")));
    jetinfo->setPtDToken(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "ptD")));
    jetinfo->setAxis2Token(consumes<edm::ValueMap<float>>(edm::InputTag(t_qgtagger, "axis2")));
    jetinfo->setMultToken(consumes<edm::ValueMap<int>>(edm::InputTag(t_qgtagger, "mult")));

    jetinfo->setUseHerwigCompatibleMatching(useHerwigCompatibleMatching);
    jetinfo->setIsHerwig(isHerwig);

    jetinfo->setGenJetMatchReclusterToken(
            consumes<edm::Association<reco::GenJetCollection> >(
                    iConfig.getParameter<edm::InputTag>( "genJetMatchRecluster" )));
    jetinfo->setGenJetMatchWithNuToken(
            consumes<edm::Association<reco::GenJetCollection> >(
                    iConfig.getParameter<edm::InputTag>( "genJetMatchWithNu" )));

    jetinfo->setGenParticlesToken(
            consumes<reco::GenParticleCollection>(
                    iConfig.getParameter<edm::InputTag>("pruned")));

    jetinfo->setMuonsToken(
            consumes<pat::MuonCollection>(
                    iConfig.getParameter<edm::InputTag>("muons")));

    jetinfo->setElectronsToken(
            consumes<pat::ElectronCollection>(
                    iConfig.getParameter<edm::InputTag>("electrons")));

    addModule(jetinfo, "jetinfo");

    ntuple_pfCands * pfcands = new ntuple_pfCands();
    pfcands->setJetRadius(jetR);

    addModule(pfcands, "pfcands");

    addModule(new ntuple_bTagVars(), "bTagVars");

    if(runFatJets_){
        auto *fatjetinfo = new ntuple_FatJetInfo(jetR);
        fatjetinfo->setGenParticleToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned")));
        fatjetinfo->setFatJetToken(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets")));
        addModule(fatjetinfo, "fatJets");
    }
    /*
     *
     * Modules initialized
     *
     * parse the input parameters (if any)
     */
    for(auto& m: modules_)
        m->getInput(iConfig);

}


DeepNtuplizer::~DeepNtuplizer()
{
    return;
    for(auto& m:modules_)
        delete m;
}


// ------------ method called for each event  ------------
void
DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //global info

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found

    edm::Handle<std::vector<reco::VertexCompositePtrCandidate> > secvertices;
    iEvent.getByToken(svToken_, secvertices);

    edm::Handle<std::vector<PileupSummaryInfo> > pupInfo;
    iEvent.getByToken(puToken_, pupInfo);

    edm::Handle<double> rhoInfo;
    iEvent.getByToken(rhoToken_,rhoInfo);

    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByToken(jetToken_, jets);

    for(auto& m:modules_){
        m->setPrimaryVertices(vertices.product());
        m->setSecVertices(secvertices.product());
        m->setPuInfo(pupInfo.product());
        m->setRhoInfo(rhoInfo.product());
        m->readSetup(iSetup);
        m->readEvent(iEvent);
    }

    std::vector<size_t> indices(jets->size());
    for(size_t i=0;i<jets->size();i++)
        indices.at(i)=i;

    if(applySelection_)
        std::random_shuffle (indices.begin(),indices.end());

    edm::View<pat::Jet>::const_iterator jetIter;
    // loop over the jets
    //for (edm::View<pat::Jet>::const_iterator jetIter = jets->begin(); jetIter != jets->end(); ++jetIter) {
    for(size_t j=0;j<indices.size();j++){
        njetstotal_++;
        size_t jetidx=indices.at(j);
        jetIter = jets->begin()+jetidx;
        const pat::Jet& jet = *jetIter;

        if(jet.genJet())
            njetswithgenjet_++;

        bool writejet=true;
				size_t idx = 0;
        for(auto& m:modules_){
					//std::cout << module_names_[idx] << std::endl;
            if(! m->fillBranches(jet, jetidx, jets.product())){
                writejet=false;
                if(applySelection_) break;
            }
						idx++;
        }
				//std::cout << "Jet done" << std::endl;
        if( (writejet&&applySelection_) || !applySelection_ ){
            tree_->Fill();
            njetsselected_++;
        }
    } // end of looping over the jets
}


// ------------ method called once each job just before starting event loop  ------------
void
DeepNtuplizer::beginJob()
{
    if( !fs ){
        throw edm::Exception( edm::errors::Configuration,
                "TFile Service is not registered in cfg file" );
    }
    tree_=(fs->make<TTree>("tree" ,"tree" ));

    for(auto& m:modules_)
        m->initBranches(tree_);

    njetstotal_=0;
    njetswithgenjet_=0;
    njetsselected_=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
DeepNtuplizer::endJob()
{

    std::cout << "total number of processed jets: " << njetstotal_<<std::endl;
    std::cout << "total number of jets with gen:  " << njetswithgenjet_<<std::endl;
    std::cout << "total number of selected jets:  "<< njetsselected_<<std::endl;
    std::cout << "fraction of selected jets:      "<< (float)njetsselected_/(float)njetstotal_<<std::endl;
    std::cout << "fraction of selected jets with gen: "<< (float)njetsselected_/(float)njetswithgenjet_<<std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
