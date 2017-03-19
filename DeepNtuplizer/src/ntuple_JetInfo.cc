/*
 * ntuple_JetInfo.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
        Modifications:George Karathanasis

 */




#include "../interface/ntuple_JetInfo.h"
#include "../interface/helpers.h"
#include <vector>
#include <algorithm>
using namespace std;

void ntuple_JetInfo::getInput(const edm::ParameterSet& iConfig){

	gluonReduction_=(iConfig.getParameter<double>("gluonReduction"));
	jetPtMin_=(iConfig.getParameter<double>("jetPtMin"));
	jetPtMax_=(iConfig.getParameter<double>("jetPtMax"));
	jetAbsEtaMin_=(iConfig.getParameter<double>("jetAbsEtaMin"));
	jetAbsEtaMax_=(iConfig.getParameter<double>("jetAbsEtaMax"));

	vector<string> disc_names = iConfig.getParameter<vector<string> >("bDiscriminators");
	for(auto& name : disc_names) {
		discriminators_[name] = 0.;
	}
}
void ntuple_JetInfo::initBranches(TTree* tree){

	//more general event info, here applied per jet
	addBranch(tree,"npv"    ,&npv_    ,"npv/f"    );
	addBranch(tree,"event_no"    ,&event_no_    ,"npv/i"    );
	addBranch(tree,"jet_no"    ,&jet_no_    ,"npv/i"    );


	// truth labels
	addBranch(tree,"gen_pt"    ,&gen_pt_    ,"gen_pt_/f"    );
	addBranch(tree,"Delta_gen_pt"    ,&Delta_gen_pt_,"Delta_gen_pt_/f"    );
	addBranch(tree,"isB",&isB_, "isB_/i");
	addBranch(tree,"isC",&isC_, "isC_/i");
	addBranch(tree,"isUDS",&isUDS_, "isUDS_/i");
	addBranch(tree,"isG",&isG_, "isG_/i");
	addBranch(tree,"isBlep",&isBlep_, "isBlep_/i");
	addBranch(tree,"isBB",&isBB_, "isBB_/i");
        addBranch(tree,"isBtau",&isBtau_,"isBtau_/i");
        addBranch(tree,"BPt",&BPt_,"BPt_/f");
        	
	// jet variables
	//b=tree->Branch("jet_pt", &jet_pt_);
	addBranch(tree,"jet_pt", &jet_pt_);

	addBranch(tree,"jet_corr_pt", &jet_corr_pt_);
	addBranch(tree,"jet_eta", &jet_eta_);

	// quark gluon
	addBranch(tree,"jet_qgl",   &jet_qgl_);  // qg tagger from jmar
	addBranch(tree,"QG_ptD",   &QG_ptD_);   // momentum fraction per jet constituent
	addBranch(tree,"QG_axis2", &QG_axis2_); // jet shape i.e. gluon are wider than quarks
	addBranch(tree,"QG_mult",  &QG_mult_);  // multiplicity i.e. total num of PFcands reconstructed
	// in the jet


	addBranch(tree,"gen_pt_Recluster"    ,&gen_pt_Recluster_    ,"gen_pt_Recluster_/f"    );
	addBranch(tree,"gen_pt_WithNu"    ,&gen_pt_WithNu_    ,"gen_pt_WithNu_/f"    );
	addBranch(tree,"Delta_gen_pt_Recluster"    ,&Delta_gen_pt_Recluster_    ,"Delta_gen_pt_Recluster_/f"    );
	addBranch(tree,"Delta_gen_pt_WithNu"    ,&Delta_gen_pt_WithNu_    ,"Delta_gen_pt_WithNu_/f"    );

	if(1) // discriminators might need to be filled differently. FIXME
		for(auto& entry : discriminators_) {
			string better_name(entry.first);
			std::replace(better_name.begin(), better_name.end(), ':', '_');
			addBranch(tree,better_name.c_str(), &entry.second, (better_name+"/F").c_str());
		}
}
void ntuple_JetInfo::readEvent(const edm::Event& iEvent){

	iEvent.getByToken(qglToken_, qglHandle);
	iEvent.getByToken(ptDToken_, ptDHandle);
	iEvent.getByToken(axis2Token_, axis2Handle);
	iEvent.getByToken(multToken_, multHandle);

	//		iEvent.getByToken(genParticleToken_, genParticleToken);
	iEvent.getByToken(genJetMatchReclusterToken_, genJetMatchRecluster);
	iEvent.getByToken(genJetMatchWithNuToken_, genJetMatchWithNu);

	//technically a branch fill but per event, therefore here
	event_no_=iEvent.id().event();
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
	for(auto& entry : discriminators_) {
		entry.second = jet.bDiscriminator(entry.first);
	}

	npv_ = vertices()->size();
	jet_no_=jetidx;


	const auto jetRef = reco::CandidatePtr(coll->ptrs().at( jetidx));
	jet_qgl_ = (*qglHandle)[jetRef];
	QG_ptD_ = (*ptDHandle)[jetRef];
	QG_axis2_ = (*axis2Handle)[jetRef];
	QG_mult_ = (*multHandle)[jetRef];

	//std::vector<Ptr<pat::Jet> > p= coll->ptrs();



	pat::JetCollection h;

	jet_pt_ = jet.correctedJet("Uncorrected").pt();
	jet_eta_ = jet.eta();
	jet_corr_pt_ = jet.pt();


	gen_pt_ =  jet.genJet()->pt();
	Delta_gen_pt_ =  jet.genJet()->pt()- jet_pt_;
  
	/*	reco::GenJet * g1;
	//g1=jet.genJet();
        g1=const_cast<reco::GenJet*>(jet.genJet());
        cout<<g1->pt()<<endl;*/
       	//for ( const reco::Candidate & gp:*prun_gen_parts() ){cout<<"skata";}
    	isBlep_=0; 
        isBB_=0; isBtau_=0; 
        isB_=0; isC_=0; isUDS_=0; isG_=0;
	switch(deep_ntuples::jet_flavour(jet,prun_gen_parts())) {
	case deep_ntuples::JetFlavor::Blep:isBlep_=1; break;
        case deep_ntuples::JetFlavor::Btau:isBtau_=1; break;
	case deep_ntuples::JetFlavor::L: isUDS_=1; break;
	case deep_ntuples::JetFlavor::BB:isBB_=1; break; 
	case deep_ntuples::JetFlavor::B: isB_=1; break;
	case deep_ntuples::JetFlavor::C: isC_=1; break;
	case deep_ntuples::JetFlavor::G: isG_=1; break;
	default : break;
	}
	//std::cout<<prun_gen_parts()->size()<<std::endl;
               

	if(!isBtau_ && !isBB_ && !isBlep_ && !isB_ && !isC_ && !isUDS_ && !isG_) return false;
        BPt_=0;
        if (isB_==1 ||isBtau_==1 || isBlep_==1||isBB_==1)
	  BPt_=deep_ntuples::parton_pt(jet,prun_gen_parts());
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

	Delta_gen_pt_Recluster_=gen_pt_Recluster_-jet.pt();
	Delta_gen_pt_WithNu_=gen_pt_WithNu_-jet.pt();


	return true;
}
