/*
 * ntuple_bTagVars.h
 *
 *  Created on: 13 Feb 2017
 *      Author: mverzett
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_BTAGVARS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_BTAGVARS_H_

#include "ntuple_content.h"
#include "TRandom3.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

/*
 * For global jet info such as eta, pt, gen info
 */
class ntuple_bTagVars: public ntuple_content{
public:
	ntuple_bTagVars():ntuple_content(){}

	void getInput(const edm::ParameterSet& iConfig);
	void initBranches(TTree* tree);
	void readEvent(const edm::Event& iEvent) {}

	//use either of these functions
	bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);

private:
	template <typename T>
	int dump_vector(reco::TaggingVariableList& from, T* to, reco::btau::TaggingVariableName name) {
		std::vector<T> vals = from.getList(name ,false);
		if(vals.size() > 0) 
			std::copy(vals.begin(), vals.end(), to);
		return vals.size();
	}

	std::string tagInfoName_;

	//jet general
	float trackJetPt_;              // track-based jet transverse momentum
	float jetNTracks_;              // tracks associated to jet
	float jetNSecondaryVertices_;   // number of secondary vertices associated to the jet
	float trackSumJetEtRatio_;      // ratio of track sum transverse energy over jet energy
	float trackSumJetDeltaR_;       // pseudoangular distance between jet axis and track fourvector sum
	float trackSip2dValAboveCharm_; // track 2D signed impact parameter of first track lifting mass above charm
	float trackSip2dSigAboveCharm_; // track 2D signed impact parameter significance of first track lifting mass above charm
	float trackSip3dValAboveCharm_; // track 3D signed impact parameter of first track lifting mass above charm
	float trackSip3dSigAboveCharm_; // track 3D signed impact parameter significance of first track lifting mass above charm
	float vertexCategory_;          // category of secondary vertex (Reco, Pseudo, No)
	//track info
	int   jetNTracksEtaRel_;        // tracks associated to jet for which trackEtaRel is calculated
	int   jetNSelectedTracks_;
	float trackMomentum_[100];    // track momentum
	float trackEta_[100];         // track pseudorapidity
	float trackPhi_[100];         // track polar angle
	float trackPtRel_[100];       // track transverse momentum, relative to the jet axis
	float trackPPar_[100];        // track parallel momentum, along the jet axis
	float trackDeltaR_[100];      // track pseudoangular distance from the jet axis
	float trackPtRatio_[100];     // track transverse momentum, relative to the jet axis, normalized to its energy
	float trackPParRatio_[100];   // track parallel momentum, along the jet axis, normalized to its energy
	float trackSip2dVal_[100];    // track 2D signed impact parameter
	float trackSip2dSig_[100];    // track 2D signed impact parameter significance
	float trackSip3dVal_[100];    // track 3D signed impact parameter
	float trackSip3dSig_[100];    // track 3D signed impact parameter significance
	float trackDecayLenVal_[100]; // track decay length
	float trackDecayLenSig_[100]; // track decay length significance
	float trackJetDistVal_[100];  // minimum track approach distance to jet axis
	float trackJetDistSig_[100];  // minimum track approach distance to jet axis significance
	float trackEtaRel_[100];      // track pseudorapidity, relative to the jet axis
	//SV info
	int   nStoredVertices_;
	float vertexMass_[10];          // mass of track sum at secondary vertex
	float vertexNTracks_[10];       // number of tracks at secondary vertex
	float vertexEnergyRatio_[10];   // ratio of energy at secondary vertex over total energy
	float vertexJetDeltaR_[10];     // pseudoangular distance between jet axis and secondary vertex direction
	float flightDistance2dVal_[10]; // transverse distance between primary and secondary vertex
	float flightDistance2dSig_[10]; // transverse distance significance between primary and secondary vertex
	float flightDistance3dVal_[10]; // distance between primary and secondary vertex
	float flightDistance3dSig_[10]; // distance significance between primary and secondary vertex
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_BTAGVARS_H_ */
