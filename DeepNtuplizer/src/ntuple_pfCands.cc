/*
 * ntuple_pfCands.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_pfCands.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

void ntuple_pfCands::getInput(const edm::ParameterSet& iConfig){

}

void ntuple_pfCands::initBranches(TTree* tree){

	tree->Branch("n_Cpfcand", &n_Cpfcand_,"n_Cpfcand_/i");
	tree->Branch("Cpfcan_pt", &Cpfcan_pt_,"Cpfcan_pt_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_phirel",&Cpfcan_phirel_,"Cpfcan_phirel_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_etarel",&Cpfcan_etarel_,"Cpfcan_etarel_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_puppiw",&Cpfcan_puppiw_,"Cpfcan_puppiw_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_dxy",&Cpfcan_dxy_,"Cpfcan_dxy_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_dz",&Cpfcan_dz_,"Cpfcan_dz_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_VTX_ass",&Cpfcan_VTX_ass_,"Cpfcan_VTX_ass_[n_Cpfcand_]/i");
	tree->Branch("Cpfcan_dptdpt",&Cpfcan_dptdpt_,"Cpfcan_dptdpt_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_detadeta",&Cpfcan_detadeta_,"Cpfcan_detadeta_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_dphidphi",&Cpfcan_dphidphi_,"Cpfcan_dphidphi_[n_Cpfcand_]/f");

	// FIXME gave INFs?
	//  tree->Branch("Cpfcan_dxydxy",&Cpfcan_dxydxy_,"Cpfcan_dxydxy_[n_Cpfcand_]/f");
	//  tree->Branch("Cpfcan_dzdz",&Cpfcan_dzdz_,"Cpfcan_dzdz_[n_Cpfcand_]/f");
	// tree->Branch("Cpfcan_dxydz",&Cpfcan_dxydz_,"Cpfcan_dxydz_[n_Cpfcand_]/f");
	// tree->Branch("Cpfcan_dphidxy",&Cpfcan_dphidxy_,"Cpfcan_dphidxy_[n_Cpfcand_]/f");
	// tree->Branch("Cpfcan_dlambdadz",&Cpfcan_dlambdadz_,"Cpfcan_dlambdadz_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_isMu",&Cpfcan_isMu_,"Cpfcan_isMu_[n_Cpfcand_]/i");
	tree->Branch("Cpfcan_isEl",&Cpfcan_isEl_,"Cpfcan_isEl_[n_Cpfcand_]/i");
	// tree->Branch("Cpfcan_lostInnerHits",&Cpfcan_lostInnerHits_,"Cpfcan_lostInnerHits_[n_Cpfcand_]/i");
	tree->Branch("Cpfcan_chi2",&Cpfcan_chi2_,"Cpfcan_chi2_[n_Cpfcand_]/f");
	tree->Branch("Cpfcan_highPurity",&Cpfcan_highPurity_,"Cpfcan_highPurity_[n_Cpfcand_]/i");

	// did not give integers !!
	//  tree->Branch("Cpfcan_charge",&Cpfcan_charge_,"Cpfcan_charge_[n_Cpfcand_]/i");

	//Neutral Pf candidates
	tree->Branch("n_Npfcand", &n_Npfcand_,"n_Npfcand_/i");
	tree->Branch("Npfcan_pt", &Npfcan_pt_,"Npfcan_pt_[n_Npfcand_]/f");
	tree->Branch("Npfcan_phirel",&Npfcan_phirel_,"Npfcan_phirel_[n_Npfcand_]/f");
	tree->Branch("Npfcan_etarel",&Npfcan_etarel_,"Npfcan_etarel_[n_Npfcand_]/f");
	tree->Branch("Npfcan_isGamma",&Npfcan_isGamma_,"Npfcan_isGamma_[n_Npfcand_]/i");
	tree->Branch("Npfcan_HadFrac",&Npfcan_HadFrac_,"Npfcan_HadFrac_[n_Npfcand_]/f");



}

void ntuple_pfCands::readEvent(const edm::Event& iEvent){

}



//use either of these functions

bool ntuple_pfCands::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

	float etasign = 1.;
	if (jet.eta()<0) etasign =-1.;

	// counts neutral and charged candicates
	n_Cpfcand_ = 0;
	n_Npfcand_ = 0;

	for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++)
	{

		/// This might include more than PF candidates, e.g. Reco muons and could
		/// be double counting. Needs to be checked.!!!!
		///
		/// Split to charged and neutral candidates

		const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
		if(PackedCandidate_->charge()!=0)
		{
			Cpfcan_pt_[n_Cpfcand_] = PackedCandidate_->pt();
			Cpfcan_phirel_[n_Cpfcand_] = reco::deltaPhi(PackedCandidate_->phi(),jet.phi());
			Cpfcan_etarel_[n_Cpfcand_] = etasign*(PackedCandidate_->eta()-jet.eta());
			Cpfcan_dxy_[n_Cpfcand_] = PackedCandidate_->dxy();
			Cpfcan_dz_[n_Cpfcand_] = PackedCandidate_->dz();
			Cpfcan_VTX_ass_[n_Cpfcand_] = PackedCandidate_->pvAssociationQuality();
			Cpfcan_puppiw_[n_Cpfcand_] = PackedCandidate_->puppiWeight();

			const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();
			reco::Track::CovarianceMatrix myCov = PseudoTrack.covariance ();
			//https://github.com/cms-sw/cmssw/blob/CMSSW_9_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h#L394
			Cpfcan_dptdpt_[n_Cpfcand_] = myCov[0][0];
			Cpfcan_detadeta_[n_Cpfcand_]= myCov[1][1];
			Cpfcan_dphidphi_[n_Cpfcand_]= myCov[2][2];
			Cpfcan_dxydxy_[n_Cpfcand_] =  myCov[3][3]; //zero if pvAssociationQuality ==7 ?
			Cpfcan_dzdz_[n_Cpfcand_] =  myCov[4][4]; //zero if pvAssociationQuality ==7 ?
			Cpfcan_dxydz_[n_Cpfcand_] =  myCov[3][4]; //zero if pvAssociationQuality ==7 ?
			Cpfcan_dphidxy_[n_Cpfcand_] =  myCov[2][3]; //zero if pvAssociationQuality ==7 ?
			Cpfcan_dlambdadz_[n_Cpfcand_] =  myCov[1][4]; //zero if pvAssociationQuality ==7 ?

			// TO DO: we can do better than that by including reco::muon informations
			Cpfcan_isMu_[n_Cpfcand_] = 0;
			if(abs(PackedCandidate_->pdgId())==13)    Cpfcan_isMu_[n_Cpfcand_] = 1;

			// TO DO: we can do better than that by including reco::electron informations
			Cpfcan_isEl_[n_Cpfcand_] = 0;
			if(abs(PackedCandidate_->pdgId())==11)    Cpfcan_isEl_[n_Cpfcand_] = 1;

			Cpfcan_charge_[n_Cpfcand_] = PackedCandidate_->charge();
			Cpfcan_lostInnerHits_[n_Cpfcand_] = PackedCandidate_->lostInnerHits();
			Cpfcan_chi2_[n_Cpfcand_] = PseudoTrack.normalizedChi2();
			Cpfcan_highPurity_[n_Cpfcand_] = PseudoTrack.highPurity;
			n_Cpfcand_++;
		}
		else{// neutral candidates
			Npfcan_pt_[n_Npfcand_] = PackedCandidate_->pt();
			Npfcan_phirel_[n_Npfcand_] = reco::deltaPhi(PackedCandidate_->phi(),jet.phi());
			Npfcan_etarel_[n_Npfcand_] = etasign*(PackedCandidate_->eta()-jet.eta());
			Npfcan_isGamma_[n_Npfcand_] = 0;
			if(fabs(PackedCandidate_->pdgId())==22)  Npfcan_isGamma_[n_Npfcand_] = 1;
			Npfcan_HadFrac_[n_Npfcand_] = PackedCandidate_->hcalFraction();
			n_Npfcand_++;
		}

	} // end loop over jet.numberOfDaughters()

	return true; //for making cuts
}
