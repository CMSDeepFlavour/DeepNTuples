/*
 * ntuple_pfCands.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_pfCands.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "../interface/sorting_modules.h"

void ntuple_pfCands::getInput(const edm::ParameterSet& iConfig){

}

void ntuple_pfCands::initBranches(TTree* tree){

	addBranch(tree,"n_Cpfcand", &n_Cpfcand_,"n_Cpfcand_/i");

	addBranch(tree,"nCpfcand", &nCpfcand_,"nCpfcand_/f");

	addBranch(tree,"Cpfcan_pt", &Cpfcan_pt_,"Cpfcan_pt_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_ptrel", &Cpfcan_ptrel_,"Cpfcan_ptrel_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_erel", &Cpfcan_erel_,"Cpfcan_erel_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_phirel",&Cpfcan_phirel_,"Cpfcan_phirel_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_etarel",&Cpfcan_etarel_,"Cpfcan_etarel_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_deltaR",&Cpfcan_deltaR_,"Cpfcan_deltaR_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_puppiw",&Cpfcan_puppiw_,"Cpfcan_puppiw_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_dxy",&Cpfcan_dxy_,"Cpfcan_dxy_[n_Cpfcand_]/f");

	addBranch(tree,"Cpfcan_dxyerr",&Cpfcan_dxyerr_,"Cpfcan_dxyerr_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_dxysig",&Cpfcan_dxysig_,"Cpfcan_dxysig_[n_Cpfcand_]/f");

	addBranch(tree,"Cpfcan_dz",&Cpfcan_dz_,"Cpfcan_dz_[n_Cpfcand_]/f");

	addBranch(tree,"Cpfcan_VTX_ass",&Cpfcan_VTX_ass_,"Cpfcan_VTX_ass_[n_Cpfcand_]/f");

	addBranch(tree,"Cpfcan_fromPV",&Cpfcan_fromPV_,"Cpfcan_fromPV_[n_Cpfcand_]/f");

	addBranch(tree,"Cpfcan_drminsv",&Cpfcan_drminsv_,"Cpfcan_drminsv_[n_Cpfcand_]/f");

//commented ones don't work
	/**///addBranch(tree,"Cpfcan_vertexChi2",&Cpfcan_vertexChi2_,"Cpfcan_vertexChi2_[n_Cpfcand_]/f");
	/**///addBranch(tree,"Cpfcan_vertexNdof",&Cpfcan_vertexNdof_,"Cpfcan_vertexNdof_[n_Cpfcand_]/f");
	/**///addBranch(tree,"Cpfcan_vertexNormalizedChi2",&Cpfcan_vertexNormalizedChi2_,"Cpfcan_vertexNormalizedChi2_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_vertex_rho",&Cpfcan_vertex_rho_,"Cpfcan_vertex_rho_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_vertex_phirel",&Cpfcan_vertex_phirel_,"Cpfcan_vertex_phirel_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_vertex_etarel",&Cpfcan_vertex_etarel_,"Cpfcan_vertex_etarel_[n_Cpfcand_]/f");
	/**///addBranch(tree,"Cpfcan_vertexRef_mass",&Cpfcan_vertexRef_mass_,"Cpfcan_vertexRef_mass_[n_Cpfcand_]/f");

	addBranch(tree,"Cpfcan_dptdpt",&Cpfcan_dptdpt_,"Cpfcan_dptdpt_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_detadeta",&Cpfcan_detadeta_,"Cpfcan_detadeta_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_dphidphi",&Cpfcan_dphidphi_,"Cpfcan_dphidphi_[n_Cpfcand_]/f");


	addBranch(tree,"Cpfcan_dxydxy",&Cpfcan_dxydxy_,"Cpfcan_dxydxy_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_dzdz",&Cpfcan_dzdz_,"Cpfcan_dzdz_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_dxydz",&Cpfcan_dxydz_,"Cpfcan_dxydz_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_dphidxy",&Cpfcan_dphidxy_,"Cpfcan_dphidxy_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_dlambdadz",&Cpfcan_dlambdadz_,"Cpfcan_dlambdadz_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_isMu",&Cpfcan_isMu_,"Cpfcan_isMu_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_isEl",&Cpfcan_isEl_,"Cpfcan_isEl_[n_Cpfcand_]/f");

	//in16 conversion broken
   // addBranch(tree,"Cpfcan_lostInnerHits",&Cpfcan_lostInnerHits_,"Cpfcan_lostInnerHits_[n_Cpfcand_]/i");
	addBranch(tree,"Cpfcan_chi2",&Cpfcan_chi2_,"Cpfcan_chi2_[n_Cpfcand_]/f");
	addBranch(tree,"Cpfcan_quality",&Cpfcan_quality_,"Cpfcan_quality_[n_Cpfcand_]/f");

	// did not give integers !!
	//  addBranch(tree,"Cpfcan_charge",&Cpfcan_charge_,"Cpfcan_charge_[n_Cpfcand_]/i");

	//Neutral Pf candidates
	addBranch(tree,"n_Npfcand", &n_Npfcand_,"n_Npfcand_/i");
	addBranch(tree,"nNpfcand", &nNpfcand_,"nNpfcand/f");

	addBranch(tree,"Npfcan_pt", &Npfcan_pt_,"Npfcan_pt_[n_Npfcand_]/f");
	addBranch(tree,"Npfcan_ptrel", &Npfcan_ptrel_,"Npfcan_ptrel_[n_Npfcand_]/f");
	addBranch(tree,"Npfcan_erel", &Npfcan_erel_,"Npfcan_erel_[n_Npfcand_]/f");

	addBranch(tree,"Npfcan_phirel",&Npfcan_phirel_,"Npfcan_phirel_[n_Npfcand_]/f");
	addBranch(tree,"Npfcan_etarel",&Npfcan_etarel_,"Npfcan_etarel_[n_Npfcand_]/f");
	addBranch(tree,"Npfcan_deltaR",&Npfcan_deltaR_,"Npfcan_deltaR_[n_Npfcand_]/f");
	addBranch(tree,"Npfcan_isGamma",&Npfcan_isGamma_,"Npfcan_isGamma_[n_Npfcand_]/f");
	addBranch(tree,"Npfcan_HadFrac",&Npfcan_HadFrac_,"Npfcan_HadFrac_[n_Npfcand_]/f");
	addBranch(tree,"Npfcan_drminsv",&Npfcan_drminsv_,"Npfcan_drminsv_[n_Npfcand_]/f");


}

void ntuple_pfCands::readEvent(const edm::Event& iEvent){

  iEvent.getByToken(svToken_, secVertices);

}



//use either of these functions

bool ntuple_pfCands::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

	float etasign = 1.;
	if (jet.eta()<0) etasign =-1.;

	std::vector<const pat::PackedCandidate* > pfcands;
	//create collection first, to be able to do some sorting
	for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++){
		const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
		if(PackedCandidate)
			pfcands.push_back(PackedCandidate);
	}



	//sort by dxy significance - many infs and nans - check  - is this the reason for worse performance? FIXME
	std::sort(pfcands.begin(),pfcands.end(),sorting::compareDxyDxyErr<const pat::PackedCandidate*>);


	// counts neutral and charged candicates
	n_Cpfcand_ = 0;
	n_Npfcand_ = 0;

	reco::VertexCompositePtrCandidateCollection cpvtx=*secVertices;


	for (const auto& PackedCandidate_:pfcands){

		if(!PackedCandidate_)continue;

		// get the dr with the closest sv
		float drminpfcandsv_ = mindrsvpfcand(cpvtx,PackedCandidate_); 

		/// This might include more than PF candidates, e.g. Reco muons and could
		/// be double counting. Needs to be checked.!!!!
		///
		/// Split to charged and neutral candidates
		if(PackedCandidate_->charge()!=0 && max_pfcand_>n_Cpfcand_){



			Cpfcan_pt_[n_Cpfcand_] = PackedCandidate_->pt();
			Cpfcan_erel_[n_Cpfcand_] = PackedCandidate_->energy()/jet.energy();
			Cpfcan_phirel_[n_Cpfcand_] = reco::deltaPhi(PackedCandidate_->phi(),jet.phi());
			Cpfcan_etarel_[n_Cpfcand_] = etasign*(PackedCandidate_->eta()-jet.eta());
			Cpfcan_deltaR_[n_Cpfcand_] =reco::deltaR(*PackedCandidate_,jet);
			Cpfcan_dxy_[n_Cpfcand_] = catchInfsAndBound(PackedCandidate_->dxy(),0,-50,50);


			Cpfcan_dxyerr_[n_Cpfcand_]=catchInfs(PackedCandidate_->dxyError(), 5.);

			Cpfcan_dxysig_[n_Cpfcand_]=catchInfsAndBound(PackedCandidate_->dxy()/PackedCandidate_->dxyError(),0.,-5000,5000);


			Cpfcan_dz_[n_Cpfcand_] = PackedCandidate_->dz();
			Cpfcan_VTX_ass_[n_Cpfcand_] = PackedCandidate_->pvAssociationQuality();

			Cpfcan_fromPV_[n_Cpfcand_] = PackedCandidate_->fromPV();

			float tempdontopt=PackedCandidate_->vx();
			tempdontopt++;

			Cpfcan_vertexChi2_[n_Cpfcand_]=PackedCandidate_->vertexChi2();
			Cpfcan_vertexNdof_[n_Cpfcand_]=PackedCandidate_->vertexNdof();
			//divided
			Cpfcan_vertexNormalizedChi2_[n_Cpfcand_]=PackedCandidate_->vertexNormalizedChi2();
			Cpfcan_vertex_rho_[n_Cpfcand_]=catchInfsAndBound(PackedCandidate_->vertex().rho(),0,-1,50);
			Cpfcan_vertex_phirel_[n_Cpfcand_]=reco::deltaPhi(PackedCandidate_->vertex().phi(),jet.phi());
			Cpfcan_vertex_etarel_[n_Cpfcand_]=etasign*(PackedCandidate_->vertex().eta()-jet.eta());
			Cpfcan_vertexRef_mass_[n_Cpfcand_]=PackedCandidate_->vertexRef()->p4().M();


			Cpfcan_puppiw_[n_Cpfcand_] = PackedCandidate_->puppiWeight();

			const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();
			reco::Track::CovarianceMatrix myCov = PseudoTrack.covariance ();
			//https://github.com/cms-sw/cmssw/blob/CMSSW_9_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h#L394

			Cpfcan_dptdpt_[n_Cpfcand_] =    catchInfsAndBound(myCov[0][0],0,-1,1);
			Cpfcan_detadeta_[n_Cpfcand_]=   catchInfsAndBound(myCov[1][1],0,-1,0.01);
			Cpfcan_dphidphi_[n_Cpfcand_]=   catchInfsAndBound(myCov[2][2],0,-1,0.1);

			/*
			 * what makes the most sense here if a track is used in the fit... cerntainly no btag
			 * for now leave it a t zero
			 * infs and nans are set to poor quality
			 */
			Cpfcan_dxydxy_[n_Cpfcand_] =    catchInfsAndBound(myCov[3][3],7.,-1,7); //zero if pvAssociationQuality ==7 ?
			Cpfcan_dzdz_[n_Cpfcand_] =      catchInfsAndBound(myCov[4][4],6.5,-1,6.5); //zero if pvAssociationQuality ==7 ?
			Cpfcan_dxydz_[n_Cpfcand_] =     catchInfsAndBound(myCov[3][4],6.,-6,6); //zero if pvAssociationQuality ==7 ?
			Cpfcan_dphidxy_[n_Cpfcand_] =   catchInfs(myCov[2][3],-0.03); //zero if pvAssociationQuality ==7 ?
			Cpfcan_dlambdadz_[n_Cpfcand_]=  catchInfs(myCov[1][4],-0.03); //zero if pvAssociationQuality ==7 ?


			// TO DO: we can do better than that by including reco::muon informations
			Cpfcan_isMu_[n_Cpfcand_] = 0;
			if(abs(PackedCandidate_->pdgId())==13) {
				Cpfcan_isMu_[n_Cpfcand_] = 1;
			}
			// TO DO: we can do better than that by including reco::electron informations
			Cpfcan_isEl_[n_Cpfcand_] = 0;
			if(abs(PackedCandidate_->pdgId())==11) {
				Cpfcan_isEl_[n_Cpfcand_] = 1;

			}

			Cpfcan_charge_[n_Cpfcand_] = PackedCandidate_->charge();
			Cpfcan_lostInnerHits_[n_Cpfcand_] = catchInfs(PackedCandidate_->lostInnerHits(),2);
			Cpfcan_chi2_[n_Cpfcand_] = catchInfsAndBound(PseudoTrack.normalizedChi2(),300,-1,300);
			//for some reason this returns the quality enum not a mask.
			Cpfcan_quality_[n_Cpfcand_] = PseudoTrack.qualityMask();

			Cpfcan_drminsv_[n_Cpfcand_] = catchInfsAndBound(drminpfcandsv_,5,-1,5);

			n_Cpfcand_++;
		}
		else if(max_pfcand_>n_Npfcand_){// neutral candidates
			Npfcan_pt_[n_Npfcand_] = PackedCandidate_->pt();
			Npfcan_ptrel_[n_Npfcand_] = PackedCandidate_->pt()/jet.pt();
			Npfcan_erel_[n_Npfcand_] = PackedCandidate_->energy()/jet.energy();
			Npfcan_phirel_[n_Npfcand_] = reco::deltaPhi(PackedCandidate_->phi(),jet.phi());
			Npfcan_etarel_[n_Npfcand_] = etasign*(PackedCandidate_->eta()-jet.eta());
			Npfcan_deltaR_[n_Npfcand_] = reco::deltaR(*PackedCandidate_,jet);
			Npfcan_isGamma_[n_Npfcand_] = 0;
			if(fabs(PackedCandidate_->pdgId())==22)  Npfcan_isGamma_[n_Npfcand_] = 1;
			Npfcan_HadFrac_[n_Npfcand_] = PackedCandidate_->hcalFraction();

			Npfcan_drminsv_[n_Npfcand_] = catchInfsAndBound(drminpfcandsv_,5,-1,5);

			n_Npfcand_++;
		}

	} // end loop over jet.numberOfDaughters()

	nCpfcand_=n_Cpfcand_;
	nNpfcand_=n_Npfcand_;

	return true; //for making cuts
}


float ntuple_pfCands::mindrsvpfcand(const std::vector<reco::VertexCompositePtrCandidate> svs, const pat::PackedCandidate* pfcand) {

  float mindr_ = 999.;
  for (unsigned int i0=0; i0<svs.size(); ++i0) {

    float tempdr_ = reco::deltaR(svs[i0],*pfcand);
    if (tempdr_<mindr_) { mindr_ = tempdr_; }

  }

  return mindr_;
}
