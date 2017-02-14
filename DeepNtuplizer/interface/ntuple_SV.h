/*
 * ntuple_SV.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_

#include "ntuple_content.h"

class ntuple_SV: public ntuple_content{
public:

	ntuple_SV();
	~ntuple_SV();

	void getInput(const edm::ParameterSet& iConfig);
	void initBranches(TTree* );
	void readEvent(const edm::Event& iEvent);

	//use either of these functions

	bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);


	void setSVToken(const edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> & t){
		svToken_=t;
	}

private:

	// SV candidates
	int   sv_num_;
	float sv_pt_[100];
	float sv_eta_[100];
	float sv_phi_[100];
	float sv_mass_[100];
	//  float sv_phirel_[100];
	//  float sv_etarel_[100];
	float sv_ntracks_[100];
	float sv_chi2_[100];
	float sv_ndf_[100];
	float sv_dxy_[100];
	float sv_dxyerr_[100];
	float sv_d3d_[100];
	float sv_d3derr_[100];
	float sv_costhetasvpv_[100];
	// TODO (comment of Markus Stoye) add information to PF candidates on which tracks were used !!!
	// this will come later LG

	//tokens to be defined from main analyzer
	edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;

	//helper:
	edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;



	//helper functions:
	Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const ;
	Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const ;
	float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const ;

};



#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_SV_H_ */
