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


// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

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
	edm::EDGetTokenT<edm::View<pat::Jet> >      jetToken_;


	edm::Service<TFileService> fs;
	TTree *tree_;


	ntuple_content * addModule(ntuple_content *m){
		modules_.push_back(m);
		return m;
	}
	std::vector<ntuple_content* > modules_;

};


DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
										  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
										  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets")))


{
	/*
	 *  Initialise the modules here
	 *  Everything else does not need to be changed if
	 *  modules don't interact.
	 */

	ntuple_SV* svmodule=new ntuple_SV();
	svmodule->setSVToken(
			consumes<reco::VertexCompositePtrCandidateCollection>(
					iConfig.getParameter<edm::InputTag>("secVertices")));
	addModule(svmodule);

	ntuple_JetInfo* jetinfo=new ntuple_JetInfo();
	jetinfo->setQglToken(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood")));
	jetinfo->setPtDToken(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood")));
	jetinfo->setAxis2Token(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2")));
	jetinfo->setMultToken(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult")));

	jetinfo->setGenJetMatchReclusterToken(
			consumes<edm::Association<reco::GenJetCollection> >(
					iConfig.getParameter<edm::InputTag>( "genJetMatchRecluster" )));
	jetinfo->setGenJetMatchWithNuToken(
			consumes<edm::Association<reco::GenJetCollection> >(
					iConfig.getParameter<edm::InputTag>( "genJetMatchWithNu" )));
	addModule(jetinfo);


	ntuple_pfCands * pfcands = new ntuple_pfCands();
	addModule(pfcands);




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


	for(auto& m:modules_){
		m->setPrimaryVertices(vertices.product());
		m->readEvent(iEvent);
	}
	edm::Handle<edm::View<pat::Jet> > jets;
	iEvent.getByToken(jetToken_, jets);



	// loop over the jets
	for (edm::View<pat::Jet>::const_iterator jetIter = jets->begin(); jetIter != jets->end(); ++jetIter) {
		const pat::Jet& jet = *jetIter;
		size_t jetidx=jetIter-jets->begin();


		bool writejet=true;
		for(auto& m:modules_){
			if(! m->fillBranches(jet, jetidx, jets.product()))
				writejet=false;
		}
		if(writejet)
			tree_->Fill();
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
}

// ------------ method called once each job just after ending the event loop  ------------
void
DeepNtuplizer::endJob()
{

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
