/*
 * ntuple_bTagVars.h
 *
 *  Created on: 13 Feb 2017
 *      Author: mverzett
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_BTAGVARS1_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_BTAGVARS1_H_


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include <iostream>
#include <math.h>
#include <iostream>
#include "TString.h"

#include "TRandom3.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

/*
 * For global jet info such as eta, pt, gen info
 */
class ntuple_bTagVars1{
public:
    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* tree);
    void readEvent(const edm::Event& iEvent) {}
    template <class T>
    void addBranch(TTree* t, const char* name,  T*, const char* leaflist=0);
    //use either of these functions
    bool fillBranches(const reco::ShallowTagInfo tagInfo, const edm::Event& iEvent,  std::vector<float> Disc);

    bool read_;
    std::vector<TString> allbranches_;



    static inline const float& catchInfs(const float& in,const float& replace_value){
      if(in==in){
	if(std::isinf(in))
	  return replace_value;
	else if(in < -1e32 || in > 1e32)
	  return replace_value;
	return in;
      }
      return replace_value;
    }

    static inline float catchInfsAndBound(const float& in,const float& replace_value,
					  const float& lowerbound, const float& upperbound,const float offset=0){
      float withoutinfs=catchInfs(in,replace_value);
      if(withoutinfs+offset<lowerbound) return lowerbound;
      if(withoutinfs+offset>upperbound) return upperbound;
      return withoutinfs;
    }


private:
    template <typename T>
    int dump_vector(reco::TaggingVariableList& from, T* to, reco::btau::TaggingVariableName name, const size_t& max) {
        std::vector<T> vals = from.getList(name ,false);
        size_t size=std::min(vals.size(),max);
        if(size > 0){
            for(size_t i=0;i<size;i++){
                to[i]=catchInfs(vals.at(i),-0.1);
            }
        }
        return size;
    }


    std::string tagInfoName_;

    //jet general
    float OnDeepCSVProbb_;
    float OnCSVProbb_;
    float OnCSVCaloProbb_;
    float OnDeepCSVProbc_;
    float OnDeepCSVProbudsg_;
    float OnDeepCSVCaloProbb_;
    float OnDeepCSVCaloProbc_;
    float OnDeepCSVCaloProbudsg_;


    float OnJetPt_;              // track-based jet transverse momentum
    float OntrackJetPt_;              // track-based jet transverse momentum
    float OnjetNTracks_;              // tracks associated to jet
    float OnjetNSecondaryVertices_;   // number of secondary vertices associated to the jet
    float OntrackSumJetEtRatio_;      // ratio of track sum transverse energy over jet energy
    float OntrackSumJetDeltaR_;       // pseudoangular distance between jet axis and track fourvector sum
    float OntrackSip2dValAboveCharm_; // track 2D signed impact parameter of first track lifting mass above charm
    float OntrackSip2dSigAboveCharm_; // track 2D signed impact parameter significance of first track lifting mass above charm
    float OntrackSip3dValAboveCharm_; // track 3D signed impact parameter of first track lifting mass above charm
    float OntrackSip3dSigAboveCharm_; // track 3D signed impact parameter significance of first track lifting mass above charm
    float OnvertexCategory_;          // category of secondary vertex (Reco, Pseudo, No)
    int   Onn_jetNTracksEtaRel_;        // tracks associated to jet for which trackEtaRel is calculated
    int   Onn_jetNSelectedTracks_;
    int   lumiBlock_;
    int   runNumber_;
    int   eventNumber_;

    float OnjetNTracksEtaRel_;
    float OnjetNSelectedTracks_;


    static constexpr size_t max_OnjetNSelectedTracks_=100;

    float OntrackMomentum_[max_OnjetNSelectedTracks_];    // track momentum
    float OntrackEta_[max_OnjetNSelectedTracks_];         // track pseudorapidity
    float OntrackPhi_[max_OnjetNSelectedTracks_];         // track polar angle
    float OntrackPtRel_[max_OnjetNSelectedTracks_];       // track transverse momentum, relative to the jet axis
    float OntrackPPar_[max_OnjetNSelectedTracks_];        // track parallel momentum, along the jet axis
    float OntrackDeltaR_[max_OnjetNSelectedTracks_];      // track pseudoangular distance from the jet axis
    float OntrackPtRatio_[max_OnjetNSelectedTracks_];     // track transverse momentum, relative to the jet axis, normalized to its energy
    float OntrackPParRatio_[max_OnjetNSelectedTracks_];   // track parallel momentum, along the jet axis, normalized to its energy
    float OntrackSip2dVal_[max_OnjetNSelectedTracks_];    // track 2D signed impact parameter
    float OntrackSip2dSig_[max_OnjetNSelectedTracks_];    // track 2D signed impact parameter significance
    float OntrackSip3dVal_[max_OnjetNSelectedTracks_];    // track 3D signed impact parameter
    float OntrackSip3dSig_[max_OnjetNSelectedTracks_];    // track 3D signed impact parameter significance
    float OntrackDecayLenVal_[max_OnjetNSelectedTracks_]; // track decay length
    float OntrackDecayLenSig_[max_OnjetNSelectedTracks_]; // track decay length significance
    float OntrackJetDistVal_[max_OnjetNSelectedTracks_];  // minimum track approach distance to jet axis
    float OntrackJetDistSig_[max_OnjetNSelectedTracks_];  // minimum track approach distance to jet axis significance
    float OntrackEtaRel_[max_OnjetNSelectedTracks_];      // track pseudorapidity, relative to the jet axis
    //SV info
    int   Onn_StoredVertices_;
    float OnNStoredVertices_;

    static constexpr size_t max_OnnStoredVertices_=10;

    float OnvertexMass_[max_OnnStoredVertices_];          // mass of track sum at secondary vertex
    float OnvertexNTracks_[max_OnnStoredVertices_];       // number of tracks at secondary vertex
    float OnvertexEnergyRatio_[max_OnnStoredVertices_];   // ratio of energy at secondary vertex over total energy
    float OnvertexJetDeltaR_[max_OnnStoredVertices_];     // pseudoangular distance between jet axis and secondary vertex direction
    float OnflightDistance2dVal_[max_OnnStoredVertices_]; // transverse distance between primary and secondary vertex
    float OnflightDistance2dSig_[max_OnnStoredVertices_]; // transverse distance significance between primary and secondary vertex
    float OnflightDistance3dVal_[max_OnnStoredVertices_]; // distance between primary and secondary vertex
    float OnflightDistance3dSig_[max_OnnStoredVertices_]; // distance significance between primary and secondary vertex
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_BTAGVARS_H_ */
