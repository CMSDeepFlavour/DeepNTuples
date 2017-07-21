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
    int   n_jetNTracksEtaRel_;        // tracks associated to jet for which trackEtaRel is calculated
    int   n_jetNSelectedTracks_;

    float jetNTracksEtaRel_;
    float jetNSelectedTracks_;


    static constexpr size_t max_jetNSelectedTracks_=100;

    float trackMomentum_[max_jetNSelectedTracks_];    // track momentum
    float trackEta_[max_jetNSelectedTracks_];         // track pseudorapidity
    float trackPhi_[max_jetNSelectedTracks_];         // track polar angle
    float trackPtRel_[max_jetNSelectedTracks_];       // track transverse momentum, relative to the jet axis
    float trackPPar_[max_jetNSelectedTracks_];        // track parallel momentum, along the jet axis
    float trackDeltaR_[max_jetNSelectedTracks_];      // track pseudoangular distance from the jet axis
    float trackPtRatio_[max_jetNSelectedTracks_];     // track transverse momentum, relative to the jet axis, normalized to its energy
    float trackPParRatio_[max_jetNSelectedTracks_];   // track parallel momentum, along the jet axis, normalized to its energy
    float trackSip2dVal_[max_jetNSelectedTracks_];    // track 2D signed impact parameter
    float trackSip2dSig_[max_jetNSelectedTracks_];    // track 2D signed impact parameter significance
    float trackSip3dVal_[max_jetNSelectedTracks_];    // track 3D signed impact parameter
    float trackSip3dSig_[max_jetNSelectedTracks_];    // track 3D signed impact parameter significance
    float trackDecayLenVal_[max_jetNSelectedTracks_]; // track decay length
    float trackDecayLenSig_[max_jetNSelectedTracks_]; // track decay length significance
    float trackJetDistVal_[max_jetNSelectedTracks_];  // minimum track approach distance to jet axis
    float trackJetDistSig_[max_jetNSelectedTracks_];  // minimum track approach distance to jet axis significance
    float trackEtaRel_[max_jetNSelectedTracks_];      // track pseudorapidity, relative to the jet axis
    //SV info
    int   n_StoredVertices_;
    float NStoredVertices_;

    static constexpr size_t max_nStoredVertices_=10;

    float vertexMass_[max_nStoredVertices_];          // mass of track sum at secondary vertex
    float vertexNTracks_[max_nStoredVertices_];       // number of tracks at secondary vertex
    float vertexEnergyRatio_[max_nStoredVertices_];   // ratio of energy at secondary vertex over total energy
    float vertexJetDeltaR_[max_nStoredVertices_];     // pseudoangular distance between jet axis and secondary vertex direction
    float flightDistance2dVal_[max_nStoredVertices_]; // transverse distance between primary and secondary vertex
    float flightDistance2dSig_[max_nStoredVertices_]; // transverse distance significance between primary and secondary vertex
    float flightDistance3dVal_[max_nStoredVertices_]; // distance between primary and secondary vertex
    float flightDistance3dSig_[max_nStoredVertices_]; // distance significance between primary and secondary vertex
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_BTAGVARS_H_ */
