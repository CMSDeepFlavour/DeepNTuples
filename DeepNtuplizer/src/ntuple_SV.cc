/*
 * ntuple_SV.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_SV.h"
// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

ntuple_SV::ntuple_SV():ntuple_content(),sv_num_(0){}
ntuple_SV::~ntuple_SV(){}


void ntuple_SV::getInput(const edm::ParameterSet& iConfig){


}
void ntuple_SV::initBranches(TTree* tree){
	// SV candidates
	tree->Branch("sv_num"         ,&sv_num_         ,"sv_num_/i"                  );
	tree->Branch("sv_pt"          ,&sv_pt_          ,"sv_pt_[sv_num_]/f"          );
	tree->Branch("sv_eta"         ,&sv_eta_         ,"sv_eta_[sv_num_]/f"         );
	tree->Branch("sv_phi"         ,&sv_phi_         ,"sv_phi_[sv_num_]/f"         );
	tree->Branch("sv_mass"        ,&sv_mass_        ,"sv_mass_[sv_num_]/f"        );
	tree->Branch("sv_ntracks"     ,&sv_ntracks_     ,"sv_ntracks_[sv_num_]/i"     );
	tree->Branch("sv_chi2"        ,&sv_chi2_        ,"sv_chi2_[sv_num_]/f"        );
	tree->Branch("sv_ndf"         ,&sv_ndf_         ,"sv_ndf_[sv_num_]/f"         );
	tree->Branch("sv_dxy"         ,&sv_dxy_         ,"sv_dxy_[sv_num_]/f"         );
	tree->Branch("sv_dxyerr"      ,&sv_dxyerr_      ,"sv_dxyerr_[sv_num_]/f"      );
	tree->Branch("sv_d3d"         ,&sv_d3d_         ,"sv_d3d_[sv_num_]/f"         );
	tree->Branch("sv_d3derr"      ,&sv_d3derr_      ,"sv_d3err_[sv_num_]/f"       );
	tree->Branch("sv_costhetasvpv",&sv_costhetasvpv_,"sv_costhetasvpv_[sv_num_]/f");
}


void ntuple_SV::readEvent(const edm::Event& iEvent){

	iEvent.getByToken(svToken_, secVertices);

}



bool ntuple_SV::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

	const reco::Vertex & pv =    vertices()->at(0);
	sv_num_ = 0;
	for (const reco::VertexCompositePtrCandidate &sv : *secVertices) {

		if (reco::deltaR(sv,jet)>0.4) { continue; }

		sv_pt_[sv_num_]           = sv.pt();
		sv_eta_[sv_num_]          = sv.eta();
		sv_phi_[sv_num_]          = sv.phi();
		sv_mass_[sv_num_]         = sv.mass();
		sv_ntracks_[sv_num_]      = sv.numberOfDaughters();
		sv_chi2_[sv_num_]         = sv.vertexChi2();
		sv_ndf_[sv_num_]          = sv.vertexNdof();
		sv_dxy_[sv_num_]          = vertexDxy(sv,pv).value();
		sv_dxyerr_[sv_num_]       = vertexDxy(sv,pv).error();
		sv_d3d_[sv_num_]          = vertexD3d(sv,pv).value();
		sv_d3derr_[sv_num_]       = vertexD3d(sv,pv).error();
		sv_costhetasvpv_[sv_num_] = vertexDdotP(sv,pv); // the pointing angle (i.e. the angle between the sum of the momentum
		// of the tracks in the SV and the flight direction betwen PV and SV)
		sv_num_++;

	} // end of looping over the secondary vertices

	return true;
}









///helpers seldomly touched



Measurement1D ntuple_SV::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const {
	VertexDistanceXY dist;
	reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
	reco::Vertex svtx(svcand.vertex(), csv);
	return dist.distance(svtx, pv);
}

Measurement1D ntuple_SV::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const {
	VertexDistance3D dist;
	reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
	reco::Vertex svtx(svcand.vertex(), csv);
	return dist.distance(svtx, pv);
}

float ntuple_SV::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const {
	reco::Candidate::Vector p = sv.momentum();
	reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
	return p.Unit().Dot(d.Unit());
}

