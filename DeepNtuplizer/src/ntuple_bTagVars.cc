/*
 * ntuple_bTagVars.cc
 *
 *      Author: mverzett
 */

#include "../interface/ntuple_bTagVars.h"
#include <sstream>
#include <algorithm>

using namespace std;

void ntuple_bTagVars::getInput(const edm::ParameterSet& iConfig){
    tagInfoName_=(iConfig.getParameter<string>("tagInfoName"));
}


void ntuple_bTagVars::initBranches(TTree* tree){
  initBranches(tree,"");
}

void ntuple_bTagVars::initBranches(TTree* tree, string On){
    //jet general
    addBranch(tree, (On + "trackJetPt").c_str()            , &trackJetPt_             , (On + "trackJetPt_/F").c_str()            );
    addBranch(tree, (On + "jetNTracks").c_str()            , &jetNTracks_             , (On + "jetNTracks_/F").c_str()            );
    addBranch(tree, (On + "TagVarCSV_jetNSecondaryVertices").c_str() , &jetNSecondaryVertices_  , (On + "jetNSecondaryVertices_/F").c_str() );
    addBranch(tree, (On + "TagVarCSV_trackSumJetEtRatio").c_str()    , &trackSumJetEtRatio_     , (On + "trackSumJetEtRatio_/F").c_str()    );
    addBranch(tree, (On + "TagVarCSV_trackSumJetDeltaR").c_str()     , &trackSumJetDeltaR_      , (On + "trackSumJetDeltaR_/F").c_str()     );
    addBranch(tree, (On + "TagVarCSV_vertexCategory").c_str()        , &vertexCategory_         , (On + "vertexCategory_/F").c_str()        );
    addBranch(tree, (On + "TagVarCSV_trackSip2dValAboveCharm").c_str(), &trackSip2dValAboveCharm_, (On + "trackSip2dValAboveCharm_/F").c_str());
    addBranch(tree, (On + "TagVarCSV_trackSip2dSigAboveCharm").c_str(), &trackSip2dSigAboveCharm_, (On + "trackSip2dSigAboveCharm_/F").c_str());
    addBranch(tree, (On + "TagVarCSV_trackSip3dValAboveCharm").c_str(), &trackSip3dValAboveCharm_, (On + "trackSip3dValAboveCharm_/F").c_str());
    addBranch(tree, (On + "TagVarCSV_trackSip3dSigAboveCharm").c_str(), &trackSip3dSigAboveCharm_, (On + "trackSip3dSigAboveCharm_/F").c_str());
    //track info
    addBranch(tree, (On + "n_TagVarCSV_jetNSelectedTracks").c_str(), &n_jetNSelectedTracks_, (On + "n_jetNSelectedTracks_/i").c_str());
    addBranch(tree, (On + "TagVarCSV_jetNSelectedTracks").c_str(), &jetNSelectedTracks_, (On + "jetNSelectedTracks_/f").c_str());

    addBranch(tree, (On + "TagVarCSVTrk_trackPtRel").c_str()     , &trackPtRel_      , (On + "trackPtRel_[n_jetNSelectedTracks_]/f").c_str()     );
    addBranch(tree, (On + "TagVarCSVTrk_trackDeltaR").c_str()    , &trackDeltaR_     , (On + "trackDeltaR_[n_jetNSelectedTracks_]/f").c_str()    );
    addBranch(tree, (On + "TagVarCSVTrk_trackPtRatio").c_str()   , &trackPtRatio_    , (On + "trackPtRatio_[n_jetNSelectedTracks_]/f").c_str()   );
    addBranch(tree, (On + "TagVarCSVTrk_trackSip3dSig").c_str()  , &trackSip3dSig_   , (On + "trackSip3dSig_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "TagVarCSVTrk_trackSip2dSig").c_str()  , &trackSip2dSig_   , (On + "trackSip2dSig_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "TagVarCSVTrk_trackDecayLenVal").c_str(), &trackDecayLenVal_, (On + "trackDecayLenVal_[n_jetNSelectedTracks_]/f").c_str());
    addBranch(tree, (On + "TagVarCSVTrk_trackJetDistVal").c_str(), &trackJetDistVal_ , (On + "trackJetDistVal_[n_jetNSelectedTracks_]/f").c_str());

    addBranch(tree, (On + "n_TagVarCSV_jetNTracksEtaRel").c_str(), &n_jetNTracksEtaRel_, (On + "n_jetNTracksEtaRel_/i").c_str()               );
    addBranch(tree, (On + "TagVarCSV_jetNTracksEtaRel").c_str(), &jetNTracksEtaRel_, (On + "jetNTracksEtaRel_/f").c_str()               );

    addBranch(tree, (On + "TagVarCSV_trackEtaRel").c_str()    , &trackEtaRel_     , (On + "trackEtaRel_[n_jetNTracksEtaRel_]/f").c_str() );

    addBranch(tree, (On + "trackPParRatio").c_str() , &trackPParRatio_  , (On + "trackPParRatio_[n_jetNSelectedTracks_]/f").c_str() );
    addBranch(tree, (On + "trackSip2dVal").c_str()  , &trackSip2dVal_   , (On + "trackSip2dVal_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "trackSip3dVal").c_str()  , &trackSip3dVal_   , (On + "trackSip3dVal_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "trackMomentum").c_str()  , &trackMomentum_   , (On + "trackMomentum_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "trackEta").c_str()       , &trackEta_        , (On + "trackEta_[n_jetNSelectedTracks_]/f").c_str()       );
    addBranch(tree, (On + "trackPPar").c_str()      , &trackPPar_       , (On + "trackPPar_[n_jetNSelectedTracks_]/f").c_str()      );
    //SV info
    addBranch(tree, (On + "n_StoredVertices").c_str()   , &n_StoredVertices_    , (On + "n_StoredVertices_/i").c_str() );
    addBranch(tree, (On + "NStoredVertices").c_str()   , &NStoredVertices_    , (On + "NStoredVertices_/f").c_str() );

    addBranch(tree, (On + "TagVarCSV_vertexMass").c_str()        , &vertexMass_         , (On + "vertexMass_[n_StoredVertices_]/f").c_str()        );
    addBranch(tree, (On + "TagVarCSV_vertexNTracks").c_str()     , &vertexNTracks_      , (On + "vertexNTracks_[n_StoredVertices_]/f").c_str()     );
    addBranch(tree, (On + "TagVarCSV_vertexEnergyRatio").c_str() , &vertexEnergyRatio_  , (On + "vertexEnergyRatio_[n_StoredVertices_]/f").c_str() );
    addBranch(tree, (On + "TagVarCSV_vertexJetDeltaR").c_str()   , &vertexJetDeltaR_    , (On + "vertexJetDeltaR_[n_StoredVertices_]/f").c_str()   );
    addBranch(tree, (On + "TagVarCSV_flightDistance2dVal").c_str(), &flightDistance2dVal_, (On + "flightDistance2dVal_[n_StoredVertices_]/f").c_str());
    addBranch(tree, (On + "TagVarCSV_flightDistance2dSig").c_str(), &flightDistance2dSig_, (On + "flightDistance2dSig_[n_StoredVertices_]/f").c_str());
    addBranch(tree, (On + "TagVarCSV_flightDistance3dVal").c_str(), &flightDistance3dVal_, (On + "flightDistance3dVal_[n_StoredVertices_]/f").c_str());
    addBranch(tree, (On + "TagVarCSV_flightDistance3dSig").c_str(), &flightDistance3dSig_, (On + "flightDistance3dSig_[n_StoredVertices_]/f").c_str());
}


//use either of these functions
bool ntuple_bTagVars::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> *){
    if(!jet.hasTagInfo(tagInfoName_)) {
        stringstream stream;
        for(auto &lab : jet.tagInfoLabels())
            stream << lab << ", ";
        throw cms::Exception("ValueError") << "There is no tagInfo embedded in the jet labelled: " << tagInfoName_ <<
                ". The available ones are: " << stream.str() << endl;
    }
    const reco::ShallowTagInfo* tagInfo = dynamic_cast<const reco::ShallowTagInfo*>(jet.tagInfo(tagInfoName_)); //to be fixed with new names

    return fillBranches(*tagInfo);
}


bool ntuple_bTagVars::fillBranches(const reco::ShallowTagInfo & tagInfo){


    reco::TaggingVariableList vars = tagInfo.taggingVariables();

    //*******************
    //
    //  jet general
    //
    //*******************
    trackJetPt_                 = vars.get(reco::btau::trackJetPt, -999);
    jetNSecondaryVertices_      = vars.get(reco::btau::jetNSecondaryVertices, -1);
    trackSumJetEtRatio_         = vars.get(reco::btau::trackSumJetEtRatio, -999);
    trackSumJetDeltaR_          = vars.get(reco::btau::trackSumJetDeltaR, -999);
    vertexCategory_             = vars.get(reco::btau::vertexCategory, -999);
    trackSip2dValAboveCharm_    = vars.get(reco::btau::trackSip2dValAboveCharm, -999);
    trackSip2dSigAboveCharm_    = vars.get(reco::btau::trackSip2dSigAboveCharm, -999);
    trackSip3dValAboveCharm_    = vars.get(reco::btau::trackSip3dValAboveCharm, -999);
    trackSip3dSigAboveCharm_    = vars.get(reco::btau::trackSip3dSigAboveCharm, -999);
    jetNTracks_ = vars.get(reco::btau::jetNTracks, -1);
    n_jetNTracksEtaRel_ = vars.get(reco::btau::jetNTracksEtaRel, -1);
    jetNTracksEtaRel_=n_jetNTracksEtaRel_;

    //*******************
    //
    //  track info  //FIXME: check that the order is the same of the charged components!
    //  FIXME: right now there are no default values in the vectors... Is this something we want? or do we want fixed size vectors with zero padding in the end?
    //
    //*******************
    n_jetNSelectedTracks_ = dump_vector(vars, trackMomentum_, reco::btau::trackMomentum,max_jetNSelectedTracks_);
    jetNSelectedTracks_=n_jetNSelectedTracks_;

    dump_vector(vars, trackEta_, reco::btau::trackEta,max_jetNSelectedTracks_)            ;
    dump_vector(vars, trackPtRel_, reco::btau::trackPtRel,max_jetNSelectedTracks_)        ;
    dump_vector(vars, trackPPar_, reco::btau::trackPPar,max_jetNSelectedTracks_)          ;
    dump_vector(vars, trackDeltaR_, reco::btau::trackDeltaR,max_jetNSelectedTracks_)      ;
    dump_vector(vars, trackPtRatio_, reco::btau::trackPtRatio,max_jetNSelectedTracks_)    ;
    dump_vector(vars, trackPParRatio_, reco::btau::trackPParRatio,max_jetNSelectedTracks_);
    dump_vector(vars, trackSip2dVal_, reco::btau::trackSip2dVal,max_jetNSelectedTracks_)  ;
    dump_vector(vars, trackSip2dSig_, reco::btau::trackSip2dSig,max_jetNSelectedTracks_)  ;
    dump_vector(vars, trackSip3dVal_, reco::btau::trackSip3dVal,max_jetNSelectedTracks_)  ;
    dump_vector(vars, trackSip3dSig_, reco::btau::trackSip3dSig,max_jetNSelectedTracks_)  ;
    dump_vector(vars, trackDecayLenVal_, reco::btau::trackDecayLenVal,max_jetNSelectedTracks_);
    dump_vector(vars, trackJetDistVal_, reco::btau::trackJetDistVal,max_jetNSelectedTracks_);
    dump_vector(vars, trackEtaRel_, reco::btau::trackEtaRel,max_jetNSelectedTracks_);

    //*******************
    //
    //  vertex info
    //
    //*******************
    n_StoredVertices_ = dump_vector(vars, vertexMass_, reco::btau::vertexMass,max_nStoredVertices_);
    NStoredVertices_=n_StoredVertices_;

    dump_vector(vars, vertexNTracks_, reco::btau::vertexNTracks,max_nStoredVertices_);
    dump_vector(vars, vertexEnergyRatio_, reco::btau::vertexEnergyRatio,max_nStoredVertices_);
    dump_vector(vars, vertexJetDeltaR_, reco::btau::vertexJetDeltaR,max_nStoredVertices_);
    dump_vector(vars, flightDistance2dVal_, reco::btau::flightDistance2dVal,max_nStoredVertices_);
    dump_vector(vars, flightDistance2dSig_, reco::btau::flightDistance2dSig,max_nStoredVertices_);
    dump_vector(vars, flightDistance3dVal_, reco::btau::flightDistance3dVal,max_nStoredVertices_);
    dump_vector(vars, flightDistance3dSig_, reco::btau::flightDistance3dSig,max_nStoredVertices_);

    return true;
}
