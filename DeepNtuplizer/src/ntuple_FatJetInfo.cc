/*
 * ntuple_FatJetInfo.cpp
 *
 *  Created on: Mar 21, 2017
 *      Author: hqu
 */

#include "../interface/ntuple_FatJetInfo.h"
#include "DeepNTuples/JetAnalysis/interface/JetHelper.h"

using namespace deep_ntuples;

void ntuple_FatJetInfo::getInput(const edm::ParameterSet& iConfig) {
	if (jetR() > 0) return;
	minSoftDropMass_ = iConfig.getUntrackedParameter<double>("minSoftDropMass", -1);
}

void ntuple_FatJetInfo::readEvent(const edm::Event& iEvent) {
	if (jetR() > 0) return;

	iEvent.getByToken(genParticleToken_, genParticleHandle_);
	iEvent.getByToken(fatjetToken_, fatjetHandle_);
}

void ntuple_FatJetInfo::initBranches(TTree* tree) {
	if (jetR() > 0) return;

	// truth labels
	addBranch(tree, "fj_isLight",        &fj_isLight_         );
	addBranch(tree, "fj_isW",            &fj_isW_             );
	addBranch(tree, "fj_isTop",          &fj_isTop_           );

	// gen-matched particle (top/W/etc.)
	addBranch(tree, "fj_gen_pt",         &fj_gen_pt_          );
	addBranch(tree, "fj_gen_eta",        &fj_gen_eta_         );

	// fatjet kinematics
	addBranch(tree, "fj_pt",             &fj_pt_              );
	addBranch(tree, "fj_eta",            &fj_eta_             );
	addBranch(tree, "fj_phi",            &fj_phi_             );
	addBranch(tree, "fj_mass",           &fj_mass_            );

	// substructure
	addBranch(tree, "fj_tau1",           &fj_tau1_            );
	addBranch(tree, "fj_tau2",           &fj_tau2_            );
	addBranch(tree, "fj_tau3",           &fj_tau3_            );
	addBranch(tree, "fj_tau21",          &fj_tau21_           );
	addBranch(tree, "fj_tau32",          &fj_tau32_           );

	// soft drop
	addBranch(tree, "fj_sdmass",         &fj_sdmass_          );

	// subjets: soft drop gives up to 2 subjets
	addBranch(tree, "fj_n_sdsubjets",    &fj_n_sdsubjets_     );
	addBranch(tree, "fj_sdsj1_pt",       &fj_sdsj1_pt_        );
	addBranch(tree, "fj_sdsj1_eta",      &fj_sdsj1_eta_       );
	addBranch(tree, "fj_sdsj1_phi",      &fj_sdsj1_phi_       );
	addBranch(tree, "fj_sdsj1_mass",     &fj_sdsj1_mass_      );
	addBranch(tree, "fj_sdsj1_csv",      &fj_sdsj1_csv_       );
	addBranch(tree, "fj_sdsj1_ptD",      &fj_sdsj1_ptD_       );
	addBranch(tree, "fj_sdsj1_axis1",    &fj_sdsj1_axis1_     );
	addBranch(tree, "fj_sdsj1_axis2",    &fj_sdsj1_axis2_     );
	addBranch(tree, "fj_sdsj1_mult",     &fj_sdsj1_mult_      );
	addBranch(tree, "fj_sdsj2_pt",       &fj_sdsj2_pt_        );
	addBranch(tree, "fj_sdsj2_eta",      &fj_sdsj2_eta_       );
	addBranch(tree, "fj_sdsj2_phi",      &fj_sdsj2_phi_       );
	addBranch(tree, "fj_sdsj2_mass",     &fj_sdsj2_mass_      );
	addBranch(tree, "fj_sdsj2_csv",      &fj_sdsj2_csv_       );
	addBranch(tree, "fj_sdsj2_ptD",      &fj_sdsj2_ptD_       );
	addBranch(tree, "fj_sdsj2_axis1",    &fj_sdsj2_axis1_     );
	addBranch(tree, "fj_sdsj2_axis2",    &fj_sdsj2_axis2_     );
	addBranch(tree, "fj_sdsj2_mult",     &fj_sdsj2_mult_      );

	// some variables used in a baseline tagger
	addBranch(tree, "fj_ptDR",           &fj_ptDR_            );
	addBranch(tree, "fj_relptdiff",      &fj_relptdiff_       );
	addBranch(tree, "fj_sdn2",           &fj_sdn2_            );
}

bool ntuple_FatJetInfo::fillBranches(const pat::Jet& jet, const size_t& jetidx, const edm::View<pat::Jet>* coll) {

	if (jetR() > 0) return true;

	const pat::Jet *fj = nullptr;
	for (const auto &j : *fatjetHandle_){
		for (const auto &sj : j.subjets()){
			if (reco::deltaR(jet, *sj)<0.1){
				fj = &j;
				break;
			}
		}
	}
	if (!fj) {
		std::cerr << "Cannot match subjet to fj!" << std::endl;
		fj_pt_ = 0;
		return true;
	}

	// preselection on sdmass and n_sd_subjets
	fj_sdmass_ = fj->userFloat("ak8PFJetsPuppiSoftDropMass");
//	if (fj_sdmass_ < minSoftDropMass_) return false;

	auto sdsubjets = fj->subjets();
	fj_n_sdsubjets_ = sdsubjets.size();

	// truth labels
	auto genmatch = fjmatch_.flavor(fj, *genParticleHandle_);
	fj_isW_     = genmatch.first == FatJetMatching::W;
	fj_isTop_   = genmatch.first == FatJetMatching::Top;
	fj_isLight_ = !(fj_isW_ || fj_isTop_);

	// gen-matched particle (top/W/etc.)
	fj_gen_pt_  = genmatch.second ? genmatch.second->pt()  : -999;
	fj_gen_eta_ = genmatch.second ? genmatch.second->eta() : -999;

	// fatjet kinematics
	fj_pt_     = fj->pt();
	fj_eta_    = fj->eta();
	fj_phi_    = fj->phi();
	fj_mass_   = fj->mass();

	// substructure
	fj_tau1_   = fj->userFloat("NjettinessAK8Puppi:tau1");
	fj_tau2_   = fj->userFloat("NjettinessAK8Puppi:tau2");
	fj_tau3_   = fj->userFloat("NjettinessAK8Puppi:tau3");
	fj_tau21_  = fj_tau1_ > 0 ? fj_tau2_/fj_tau1_ : 1.01;
	fj_tau32_  = fj_tau2_ > 0 ? fj_tau3_/fj_tau2_ : 1.01;

	// subjets: soft drop gives up to 2 subjets
	if (fj_n_sdsubjets_ < 2){
		fj_sdsj1_pt_    = 0;
		fj_sdsj1_eta_   = 0;
		fj_sdsj1_phi_   = 0;
		fj_sdsj1_mass_  = 0;
		fj_sdsj1_csv_   = 0;
		fj_sdsj1_ptD_   = 0;
		fj_sdsj1_axis1_ = 0;
		fj_sdsj1_axis2_ = 0;
		fj_sdsj1_mult_  = 0;

		fj_sdsj2_pt_    = 0;
		fj_sdsj2_eta_   = 0;
		fj_sdsj2_phi_   = 0;
		fj_sdsj2_mass_  = 0;
		fj_sdsj2_csv_   = 0;
		fj_sdsj2_ptD_   = 0;
		fj_sdsj2_axis1_ = 0;
		fj_sdsj2_axis2_ = 0;
		fj_sdsj2_mult_  = 0;

		// some variables used in a baseline tagger
		fj_ptDR_        = 0;
		fj_relptdiff_   = 0;
		fj_sdn2_        = 0;
	} else {
		const auto& sj1 = *sdsubjets.at(0);
		JetHelper jh1(&sj1);
		fj_sdsj1_pt_    = sj1.pt();
		fj_sdsj1_eta_   = sj1.eta();
		fj_sdsj1_phi_   = sj1.phi();
		fj_sdsj1_mass_  = sj1.mass();
		fj_sdsj1_csv_   = sj1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		fj_sdsj1_ptD_   = jh1.ptD();
		fj_sdsj1_axis1_ = jh1.axis1();
		fj_sdsj1_axis2_ = jh1.axis2();
		fj_sdsj1_mult_  = jh1.mult();

		const auto& sj2 = *sdsubjets.at(1);
		JetHelper jh2(&sj2);
		fj_sdsj2_pt_    = sj2.pt();
		fj_sdsj2_eta_   = sj2.eta();
		fj_sdsj2_phi_   = sj2.phi();
		fj_sdsj2_mass_  = sj2.mass();
		fj_sdsj2_csv_   = sj2.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		fj_sdsj2_ptD_   = jh2.ptD();
		fj_sdsj2_axis1_ = jh2.axis1();
		fj_sdsj2_axis2_ = jh2.axis2();
		fj_sdsj2_mult_  = jh2.mult();

		// some variables used in a baseline tagger
		float deltaR = reco::deltaR(sj1, sj2);
		float var_sd_0 = sj2.pt()/(sj1.pt()+sj2.pt());
		fj_ptDR_        = fj->pt() * deltaR;
		fj_relptdiff_   = std::abs(sj1.pt()-sj2.pt())/fj->pt();
		fj_sdn2_        = var_sd_0/std::pow(deltaR,-2);
	}

	return true;
}
