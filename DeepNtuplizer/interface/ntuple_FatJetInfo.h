/*
 * ntuple_FatJetInfo.h
 *
 *  Created on: Mar 21, 2017
 *      Author: hqu
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_FATJETINFO_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_FATJETINFO_H_

#include "ntuple_content.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "DeepNTuples/JetAnalysis/interface/FatJetMatching.h"

/*
 * For fatjet specific variables. NOT FILLED for ak4.
 */

class ntuple_FatJetInfo : public ntuple_content{
public:
	ntuple_FatJetInfo() : ntuple_FatJetInfo(0.4) {}
	ntuple_FatJetInfo(double jetR) : ntuple_content(jetR), fjmatch_(jetR>0?jetR:0.8, true) {}
	virtual ~ntuple_FatJetInfo() {}

	void getInput(const edm::ParameterSet& iConfig) override;
	void readEvent(const edm::Event& iEvent) override;

	void initBranches(TTree* tree) override;
	bool fillBranches(const pat::Jet &jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll) override;

	void setGenParticleToken(const edm::EDGetTokenT<reco::GenParticleCollection>& genPartToken) {
		genParticleToken_ = genPartToken;
	}
	void setFatJetToken(const edm::EDGetTokenT<pat::JetCollection>& fatjetToken) {
		fatjetToken_ = fatjetToken;
	}

private:
	double minSoftDropMass_ = 0;
	std::string tagInfoName_ ;
	deep_ntuples::FatJetMatching fjmatch_;
	edm::EDGetTokenT<reco::GenParticleCollection>      genParticleToken_;
	edm::Handle<reco::GenParticleCollection>           genParticleHandle_;
	edm::EDGetTokenT<pat::JetCollection>               fatjetToken_;
	edm::Handle<pat::JetCollection>                    fatjetHandle_;

private:
	// truth labels
	int fj_isLight_ = 0;
	int fj_isW_     = 0;
	int fj_isH_     = 0;
	int fj_isTop_   = 0;

	// gen-matched particle (top/W/etc.)
	float fj_gen_pt_  = 0;
	float fj_gen_eta_ = 0;
	
	double flavour_ = -99;
 	double nbHadrons_  = -99;
 	double ncHadrons_  = -99;

	// fatjet kinematics
	float fj_pt_     = 0;
	float fj_eta_    = 0;
	float fj_phi_    = 0;
	float fj_mass_   = 0;

	// substructure
	float fj_tau1_   = 0;
	float fj_tau2_   = 0;
	float fj_tau3_   = 0;
	float fj_tau21_  = 0;
	float fj_tau32_  = 0;

	// soft drop
	float fj_sdmass_ = 0;
	
	//double-b 
	float fj_doubleb_ =-99;

	//double-b inputs
	float z_ratio_ = -99;
	float trackSipdSig_3_ = -99;
        float trackSipdSig_2_ = -99;
        float trackSipdSig_1_ = -99;
        float trackSipdSig_0_ = -99;
        float trackSipdSig_1_0_ = -99;
        float trackSipdSig_0_0_ = -99;
        float trackSipdSig_1_1_ = -99;
        float trackSipdSig_0_1_ = -99;
        float trackSip2dSigAboveCharm_0_ =  -99;
        float trackSip2dSigAboveBottom_0_ =  -99;
        float trackSip2dSigAboveBottom_1_ =  -99;
        float tau1_trackEtaRel_0_ =  -99;
        float tau1_trackEtaRel_1_ = -99;
        float tau1_trackEtaRel_2_ = -99;
        float tau0_trackEtaRel_0_ = -99;
        float tau0_trackEtaRel_1_ = -99;
        float tau0_trackEtaRel_2_ = -99;
        float tau_vertexMass_0_ = -99;
        float tau_vertexEnergyRatio_0_ = -99;
        float tau_vertexDeltaR_0_ = -99;
        float tau_flightDistance2dSig_0_ = -99;
        float tau_vertexMass_1_ = -99;
        float tau_vertexEnergyRatio_1_ = -99;
        float tau_flightDistance2dSig_1_ =   -99;
        float jetNTracks_ =  -99; 
        float nSV_ =  -99;
	

	// subjets: soft drop gives up to 2 subjets
	float fj_n_sdsubjets_ = 0;

	float fj_sdsj1_pt_    = 0;
	float fj_sdsj1_eta_   = 0;
	float fj_sdsj1_phi_   = 0;
	float fj_sdsj1_mass_  = 0;
	float fj_sdsj1_csv_   = 0;
	float fj_sdsj1_ptD_   = 0;
	float fj_sdsj1_axis1_ = 0;
	float fj_sdsj1_axis2_ = 0;
	float fj_sdsj1_mult_  = 0;

	float fj_sdsj2_pt_    = 0;
	float fj_sdsj2_eta_   = 0;
	float fj_sdsj2_phi_   = 0;
	float fj_sdsj2_mass_  = 0;
	float fj_sdsj2_csv_   = 0;
	float fj_sdsj2_ptD_   = 0;
	float fj_sdsj2_axis1_ = 0;
	float fj_sdsj2_axis2_ = 0;
	float fj_sdsj2_mult_  = 0;

	// some variables used in a baseline tagger
	float fj_ptDR_        = 0;
	float fj_relptdiff_   = 0;
	float fj_sdn2_        = 0;

};

#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_FATJETINFO_H_ */
