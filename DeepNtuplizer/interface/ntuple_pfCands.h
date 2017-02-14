/*
 * ntuple_pfcands.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_

#include "ntuple_content.h"

class ntuple_pfCands: public ntuple_content{
public:

	ntuple_pfCands():ntuple_content(){}

	void getInput(const edm::ParameterSet& iConfig);
	void initBranches(TTree* );
	void readEvent(const edm::Event& iEvent);



	//use either of these functions

	bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);

private:
	unsigned int n_Cpfcand_;
	float  Cpfcan_pt_[100];
	float  Cpfcan_phirel_[100];
	float  Cpfcan_etarel_[100];
	float  Cpfcan_puppiw_[100];
	float   Cpfcan_VTX_ass_[100];

	// covariance
	float  Cpfcan_dz_[100];
	float  Cpfcan_dxy_[100];
	float  Cpfcan_dptdpt_[100];
	float  Cpfcan_detadeta_[100];
	float  Cpfcan_dphidphi_[100];
	float  Cpfcan_dxydxy_[100];
	float  Cpfcan_dzdz_[100];
	float  Cpfcan_dxydz_[100];
	float  Cpfcan_dphidxy_[100];
	float  Cpfcan_dlambdadz_[100];

	// ID, skipped "charged hadron" as that is true if now the other
	// TODO (comment of Markus Stoye) add reco information
	float Cpfcan_isMu_[100]; // pitty that the quality is missing
	float Cpfcan_isEl_[100]; // pitty that the quality is missing
	float Cpfcan_charge_[100];

	// track quality
	float Cpfcan_lostInnerHits_[100];
	float Cpfcan_chi2_[100];
	float Cpfcan_highPurity_[100];

	//Neutral Pf candidates
	int n_Npfcand_;
	float  Npfcan_pt_[100];
	float  Npfcan_phirel_[100];
	float  Npfcan_etarel_[100];
	float  Npfcan_isGamma_[100];
	float  Npfcan_HadFrac_[100];

};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_ */
