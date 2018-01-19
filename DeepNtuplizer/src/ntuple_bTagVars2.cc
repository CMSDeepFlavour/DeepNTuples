/*
 * ntuple_bTagVars.cc
 *
 *      Author: mverzett
 */

#include "../interface/ntuple_bTagVars2.h"
#include <sstream>
#include <algorithm>

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

using namespace std;

void ntuple_bTagVars2::getInput(const edm::ParameterSet& iConfig){
  //tagInfoName_=(iConfig.getParameter<string>("tagInfoName"));
}
template <class T>
void ntuple_bTagVars2::addBranch(TTree* t, const char* name,  T* address, const char* leaflist){
  read_ = false;
  if(read_ ){
    t->SetBranchAddress(name,address);
  }
  else{
    if(leaflist)
      t->Branch(name  ,address  ,leaflist );
    else
      t->Branch(name  ,address);
  }
  //  allbranches_.push_back((TString)name);

}

void ntuple_bTagVars2::initBranches(TTree* tree){

    addBranch(tree,"DeepCSVProbb"             , &DeepCSVProbb_             , "DeepCSVProbb_/F"             );
    addBranch(tree,"DeepCSVProbc"             , &DeepCSVProbc_             , "DeepCSVProbc_/F"             );
    addBranch(tree,"DeepCSVProbudsg"             , &DeepCSVProbudsg_             , "DeepCSVProbudsg_/F"             );
    addBranch(tree,"CSVProbb"             , &CSVProbb_             , "CSVProbb_/F"             );


    addBranch(tree,"JetPt"             , &JetPt_             , "JetPt_/F"             );
    addBranch(tree,"trackJetPt"             , &trackJetPt_             , "trackJetPt_/F"             );
    addBranch(tree,"jetNTracks"             , &jetNTracks_             , "jetNTracks_/F"             );
    addBranch(tree,"TagVarCSV_jetNSecondaryVertices"  , &jetNSecondaryVertices_  , "jetNSecondaryVertices_/F"  );
    addBranch(tree,"TagVarCSV_trackSumJetEtRatio"     , &trackSumJetEtRatio_     , "trackSumJetEtRatio_/F"     );
    addBranch(tree,"TagVarCSV_trackSumJetDeltaR"      , &trackSumJetDeltaR_      , "trackSumJetDeltaR_/F"      );
    addBranch(tree,"TagVarCSV_vertexCategory"         , &vertexCategory_         , "vertexCategory_/F"         );
    addBranch(tree,"TagVarCSV_trackSip2dValAboveCharm", &trackSip2dValAboveCharm_, "trackSip2dValAboveCharm_/F");
    addBranch(tree,"TagVarCSV_trackSip2dSigAboveCharm", &trackSip2dSigAboveCharm_, "trackSip2dSigAboveCharm_/F");
    addBranch(tree,"TagVarCSV_trackSip3dValAboveCharm", &trackSip3dValAboveCharm_, "trackSip3dValAboveCharm_/F");
    addBranch(tree,"TagVarCSV_trackSip3dSigAboveCharm", &trackSip3dSigAboveCharm_, "trackSip3dSigAboveCharm_/F");
    //track info
    addBranch(tree,"n_TagVarCSV_jetNSelectedTracks", &n_jetNSelectedTracks_, "n_jetNSelectedTracks_/i");
    addBranch(tree,"TagVarCSV_jetNSelectedTracks", &jetNSelectedTracks_, "jetNSelectedTracks_/f");

    addBranch(tree,"TagVarCSVTrk_trackPtRel"      , &trackPtRel_      , "trackPtRel_[n_jetNSelectedTracks_]/f"      );
    addBranch(tree,"TagVarCSVTrk_trackDeltaR"     , &trackDeltaR_     , "trackDeltaR_[n_jetNSelectedTracks_]/f"     );
    addBranch(tree,"TagVarCSVTrk_trackPtRatio"    , &trackPtRatio_    , "trackPtRatio_[n_jetNSelectedTracks_]/f"    );
    addBranch(tree,"TagVarCSVTrk_trackSip3dSig"   , &trackSip3dSig_   , "trackSip3dSig_[n_jetNSelectedTracks_]/f"   );
    addBranch(tree,"TagVarCSVTrk_trackSip2dSig"   , &trackSip2dSig_   , "trackSip2dSig_[n_jetNSelectedTracks_]/f"   );
    addBranch(tree,"TagVarCSVTrk_trackDecayLenVal", &trackDecayLenVal_, "trackDecayLenVal_[n_jetNSelectedTracks_]/f");
    addBranch(tree,"TagVarCSVTrk_trackJetDistVal" , &trackJetDistVal_ , "trackJetDistVal_[n_jetNSelectedTracks_]/f" );

    addBranch(tree,"n_TagVarCSV_jetNTracksEtaRel", &n_jetNTracksEtaRel_, "n_jetNTracksEtaRel_/i"                );
    addBranch(tree,"TagVarCSV_jetNTracksEtaRel", &jetNTracksEtaRel_, "jetNTracksEtaRel_/f"                );

    addBranch(tree,"TagVarCSV_trackEtaRel"     , &trackEtaRel_     , "trackEtaRel_[n_jetNTracksEtaRel_]/f"  );

    addBranch(tree,"trackPParRatio"  , &trackPParRatio_  , "trackPParRatio_[n_jetNSelectedTracks_]/f"  );
    addBranch(tree,"trackSip2dVal"   , &trackSip2dVal_   , "trackSip2dVal_[n_jetNSelectedTracks_]/f"   );
    addBranch(tree,"trackSip3dVal"   , &trackSip3dVal_   , "trackSip3dVal_[n_jetNSelectedTracks_]/f"   );
    addBranch(tree,"trackMomentum"   , &trackMomentum_   , "trackMomentum_[n_jetNSelectedTracks_]/f"   );
    addBranch(tree,"trackEta"        , &trackEta_        , "trackEta_[n_jetNSelectedTracks_]/f"        );
    addBranch(tree,"trackPPar"       , &trackPPar_       , "trackPPar_[n_jetNSelectedTracks_]/f"       );
    //SV info
    addBranch(tree,"n_StoredVertices"    , &n_StoredVertices_    , "n_StoredVertices_/i"  );
    addBranch(tree,"NStoredVertices"    , &NStoredVertices_    , "NStoredVertices_/f"  );

    addBranch(tree,"TagVarCSV_vertexMass"         , &vertexMass_         , "vertexMass_[n_StoredVertices_]/f"         );
    addBranch(tree,"TagVarCSV_vertexNTracks"      , &vertexNTracks_      , "vertexNTracks_[n_StoredVertices_]/f"      );
    addBranch(tree,"TagVarCSV_vertexEnergyRatio"  , &vertexEnergyRatio_  , "vertexEnergyRatio_[n_StoredVertices_]/f"  );
    addBranch(tree,"TagVarCSV_vertexJetDeltaR"    , &vertexJetDeltaR_    , "vertexJetDeltaR_[n_StoredVertices_]/f"    );
    addBranch(tree,"TagVarCSV_flightDistance2dVal", &flightDistance2dVal_, "flightDistance2dVal_[n_StoredVertices_]/f");
    addBranch(tree,"TagVarCSV_flightDistance2dSig", &flightDistance2dSig_, "flightDistance2dSig_[n_StoredVertices_]/f");
    addBranch(tree,"TagVarCSV_flightDistance3dVal", &flightDistance3dVal_, "flightDistance3dVal_[n_StoredVertices_]/f");
    addBranch(tree,"TagVarCSV_flightDistance3dSig", &flightDistance3dSig_, "flightDistance3dSig_[n_StoredVertices_]/f");
}

//use either of these functions
bool ntuple_bTagVars2::fillBranches(const reco::ShallowTagInfo tagInfo, const edm::Event& iEvent, std::vector<float> Disc){
 
    reco::TaggingVariableList vars = tagInfo.taggingVariables();

    //*******************
    //
    //  jet general
    //
    //*******************
    DeepCSVProbb_ = Disc.at(0);
    DeepCSVProbc_ = Disc.at(1);
    DeepCSVProbudsg_ = Disc.at(2);
    CSVProbb_ = Disc.at(3);
    JetPt_                 = vars.get(reco::btau::jetPt, -999);
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
