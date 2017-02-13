/*
 * ntuple_JetInfo.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */




#include "../interface/ntuple_JetInfo.h"



void ntuple_JetInfo::getInput(const edm::ParameterSet& iConfig){

	gluonReduction_=(iConfig.getParameter<double>("gluonReduction"));
	jetPtMin_=(iConfig.getParameter<double>("jetPtMin"));
	jetPtMax_=(iConfig.getParameter<double>("jetPtMax"));
	jetAbsEtaMin_=(iConfig.getParameter<double>("jetAbsEtaMin"));
	jetAbsEtaMax_=(iConfig.getParameter<double>("jetAbsEtaMax"));
}
void ntuple_JetInfo::initBranches(TTree* tree){

	// truthe labels
	tree->Branch("gen_pt"    ,&gen_pt_    ,"gen_pt_/f"    );
	tree->Branch("Delta_gen_pt"    ,&Delta_gen_pt_,"Delta_gen_pt_/f"    );
	tree->Branch("isB",&isB_, "isB_/i");
	tree->Branch("isC",&isC_, "isC_/i");
	tree->Branch("isUDS",&isUDS_, "isUDS_/i");
	tree->Branch("isG",&isG_, "isG_/i");

	tree->Branch("npv"    ,&npv_    ,"npv/i"    );

	// jet variables
	tree->Branch("jet_pt", &jet_pt_);
	tree->Branch("jet_eta", &jet_eta_);

	// quark gluon
	tree->Branch("jet_qgl",   &jet_qgl_);
	tree->Branch("jet_ptD",   &jet_ptD_);
	tree->Branch("jet_axis2", &jet_axis2_);
	tree->Branch("jet_mult",  &jet_mult_);


	tree->Branch("gen_pt_Recluster"    ,&gen_pt_Recluster_    ,"gen_pt_Recluster_/f"    );
	tree->Branch("gen_pt_WithNu"    ,&gen_pt_WithNu_    ,"gen_pt_WithNu_/f"    );

}
void ntuple_JetInfo::readEvent(const edm::Event& iEvent){

	iEvent.getByToken(qglToken_, qglHandle);
	iEvent.getByToken(ptDToken_, ptDHandle);
	iEvent.getByToken(axis2Token_, axis2Handle);
	iEvent.getByToken(multToken_, multHandle);


	iEvent.getByToken(genJetMatchReclusterToken_, genJetMatchRecluster);
	iEvent.getByToken(genJetMatchWithNuToken_, genJetMatchWithNu);
}

//use either of these functions

bool ntuple_JetInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
	if(!coll)
		throw std::runtime_error("ntuple_JetInfo::fillBranches: no jet collection");

	/// cuts ///

	// some cuts to contrin training region
	if ( jet.pt() < jetPtMin_ ||  jet.pt() > jetPtMax_ ) return false;                  // apply jet pT cut
	if ( jet.eta() < fabs(jetAbsEtaMin_) ||jet.eta() > fabs(jetAbsEtaMax_) ) return false; // apply jet eta cut

	// often we have way to many gluons that we do not need. This randomply reduces the gluons
	if (gluonReduction_>0 && jet.partonFlavour()==21)
		if(TRandom_.Uniform()>gluonReduction_) return false;

	if(jet.genJet()==NULL)return false;


	//branch fills

	npv_ = vertices()->size();

	const auto jetRef = reco::CandidatePtr(coll->ptrs().at( jetidx));
	jet_qgl_ = (*qglHandle)[jetRef];
	jet_ptD_ = (*ptDHandle)[jetRef];
	jet_axis2_ = (*axis2Handle)[jetRef];
	jet_mult_ = (*multHandle)[jetRef];

	//std::vector<Ptr<pat::Jet> > p= coll->ptrs();


	///re check! FIXME
	isB_= int(abs(jet.partonFlavour())==5);
	isC_= int(abs(jet.partonFlavour())==4);
	isUDS_= int( (abs(jet.partonFlavour())>0) && (abs(jet.partonFlavour())<4 ));
	isG_= int(jet.partonFlavour()==21);

	// check for sumP = 1 if() return false;
	pat::JetCollection h;

	jet_pt_ = jet.pt();
	jet_eta_ = jet.eta();


	gen_pt_ =  jet.genJet()->pt();
	Delta_gen_pt_ =  jet.genJet()->pt()- jet.pt();


	const edm::RefToBase<pat::Jet> patJetRef = coll->refAt(jetidx);
	reco::GenJetRef genjetRecluster = (*genJetMatchRecluster)[patJetRef];

	gen_pt_Recluster_ = 0.;
	if (genjetRecluster.isNonnull() && genjetRecluster.isAvailable()) {
		gen_pt_Recluster_ = genjetRecluster->pt();
	}
	reco::GenJetRef genjetWithNu = (*genJetMatchWithNu)[patJetRef];

	gen_pt_WithNu_ = 0.;
	if (genjetWithNu.isNonnull() && genjetWithNu.isAvailable()) {
		gen_pt_WithNu_ = genjetWithNu->pt();
	}




	return true;
}
