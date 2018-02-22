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
    fillBranches(*tagInfo);

    return true;
}

bool ntuple_bTagVars::fillBranches(const reco::ShallowTagInfo &tagInfo){


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

bool ntuple_bTagVars::Copy(treeReader & Reader,int & jet){

    trackJetPt_                 = Reader.trackJetPt_[jet];
    jetNSecondaryVertices_      = Reader.jetNSecondaryVertices_[jet];
    trackSumJetEtRatio_         = Reader.trackSumJetEtRatio_[jet];
    trackSumJetDeltaR_          = Reader.trackSumJetDeltaR_[jet];
    vertexCategory_             = Reader.vertexCategory_[jet];
    trackSip2dValAboveCharm_    = Reader.trackSip2dValAboveCharm_[jet];
    trackSip2dSigAboveCharm_    = Reader.trackSip2dSigAboveCharm_[jet];
    trackSip3dValAboveCharm_    = Reader.trackSip3dValAboveCharm_[jet];
    trackSip3dSigAboveCharm_    = Reader.trackSip3dSigAboveCharm_[jet];
    jetNTracks_ = Reader.jetNTracks_[jet];
    jetNTracksEtaRel_ = Reader.jetNTracksEtaRel_[jet];
    n_jetNTracksEtaRel_=jetNTracksEtaRel_;
    jetNSelectedTracks_ = Reader.jetNSelectedTracks_[jet];
    n_jetNSelectedTracks_ = jetNSelectedTracks_;
    int TrkIndex = Reader.Jet_nFirstTrkTagVarCSV_[jet];
    for(int z = 0; z<n_jetNSelectedTracks_; z++){
      trackEta_[z] = Reader.trackEta_[TrkIndex];
      trackPtRel_[z] = Reader.trackPtRel_[TrkIndex];
      trackPPar_[z] = Reader.trackPPar_[TrkIndex];
      trackDeltaR_[z] = Reader.trackDeltaR_[TrkIndex];
      trackPtRatio_[z] = Reader.trackPtRatio_[TrkIndex];
      trackPParRatio_[z] = Reader.trackPParRatio_[TrkIndex];
      trackSip2dVal_[z] = Reader.trackSip2dVal_[TrkIndex];
      trackSip2dSig_[z] = Reader.trackSip2dSig_[TrkIndex];
      trackSip3dVal_[z] = Reader.trackSip3dVal_[TrkIndex];
      trackSip3dSig_[z] = Reader.trackSip3dSig_[TrkIndex];
      trackDecayLenVal_[z] = Reader.trackDecayLenVal_[TrkIndex];
      trackJetDistVal_[z] = Reader.trackJetDistVal_[TrkIndex];
      TrkIndex++;
    }
    if((Reader.vertexCategory_[jet] == 0) | (Reader.vertexCategory_[jet] == 1)){
      NStoredVertices_ = 1;
      n_StoredVertices_ = 1;
    }
    else{
      NStoredVertices_ = 0;
      n_StoredVertices_ = 0;    
    }
    if((Reader.vertexCategory_[jet] == 0) | (Reader.vertexCategory_[jet] == 1)){
      vertexMass_[0] = Reader.vertexMass_[jet];
      vertexNTracks_[0] = Reader.vertexNTracks_[jet];
      vertexEnergyRatio_[0] = Reader.vertexEnergyRatio_[jet];
      vertexJetDeltaR_[0] = Reader.vertexJetDeltaR_[jet];
      flightDistance2dVal_[0] = Reader.flightDistance2dVal_[jet];
      flightDistance2dSig_[0] = Reader.flightDistance2dSig_[jet];
      flightDistance3dVal_[0] = Reader.flightDistance3dVal_[jet];
      flightDistance3dSig_[0] = Reader.flightDistance3dSig_[jet];
    }

    int EtaRelIndex = Reader.Jet_nFirstTrkEtaRelTagVarCSV_[jet];
    for(int z = 0; z<n_jetNTracksEtaRel_; z++){
      trackEtaRel_[z] = Reader.trackEtaRel_[EtaRelIndex];
      EtaRelIndex++;
    }
    return true;
}

