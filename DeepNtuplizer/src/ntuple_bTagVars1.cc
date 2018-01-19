/*
 * ntuple_bTagVars.cc
 *
 *      Author: mverzett
 */

#include "../interface/ntuple_bTagVars1.h"
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

void ntuple_bTagVars1::getInput(const edm::ParameterSet& iConfig){
  //tagInfoName_=(iConfig.getParameter<string>("tagInfoName"));
}
template <class T>
void ntuple_bTagVars1::addBranch(TTree* t, const char* name,  T* address, const char* leaflist){
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

void ntuple_bTagVars1::initBranches(TTree* tree){
    //jet general
    addBranch(tree,"lumiBlock"             , &lumiBlock_             , "lumiBlock_/i"             );
    addBranch(tree,"runNumber"             , &runNumber_             , "runNumber_/i"             );
    addBranch(tree,"eventNumber"             , &eventNumber_             , "eventNumber_/i"             );
    addBranch(tree,"OnDeepCSVProbb"             , &OnDeepCSVProbb_             , "OnDeepCSVProbb_/F"             );
    addBranch(tree,"OnCSVProbb"             , &OnCSVProbb_             , "OnCSVProbb_/F"             );
    addBranch(tree,"OnCSVCaloProbb"             , &OnCSVCaloProbb_             , "OnCSVCaloProbb_/F"             );
    addBranch(tree,"OnDeepCSVProbc"             , &OnDeepCSVProbc_             , "OnDeepCSVProbc_/F"             );
    addBranch(tree,"OnDeepCSVProbudsg"             , &OnDeepCSVProbudsg_             , "OnDeepCSVProbudsg_/F"             );
    addBranch(tree,"OnDeepCSVCaloProbb"             , &OnDeepCSVCaloProbb_             , "OnDeepCSVCaloProbb_/F"             );
    addBranch(tree,"OnDeepCSVCaloProbc"             , &OnDeepCSVCaloProbc_             , "OnDeepCSVCaloProbc_/F"             );
    addBranch(tree,"OnDeepCSVCaloProbudsg"             , &OnDeepCSVCaloProbudsg_             , "OnDeepCSVCaloProbudsg_/F"             );


    addBranch(tree,"OnJetPt"             , &OnJetPt_             , "OnJetPt_/F"             );
    addBranch(tree,"OntrackJetPt"             , &OntrackJetPt_             , "OntrackJetPt_/F"             );
    addBranch(tree,"OnjetNTracks"             , &OnjetNTracks_             , "OnjetNTracks_/F"             );
    addBranch(tree,"OnTagVarCSV_jetNSecondaryVertices"  , &OnjetNSecondaryVertices_  , "OnjetNSecondaryVertices_/F"  );
    addBranch(tree,"OnTagVarCSV_trackSumJetEtRatio"     , &OntrackSumJetEtRatio_     , "OntrackSumJetEtRatio_/F"     );
    addBranch(tree,"OnTagVarCSV_trackSumJetDeltaR"      , &OntrackSumJetDeltaR_      , "OntrackSumJetDeltaR_/F"      );
    addBranch(tree,"OnTagVarCSV_vertexCategory"         , &OnvertexCategory_         , "OnvertexCategory_/F"         );
    addBranch(tree,"OnTagVarCSV_trackSip2dValAboveCharm", &OntrackSip2dValAboveCharm_, "OntrackSip2dValAboveCharm_/F");
    addBranch(tree,"OnTagVarCSV_trackSip2dSigAboveCharm", &OntrackSip2dSigAboveCharm_, "OntrackSip2dSigAboveCharm_/F");
    addBranch(tree,"OnTagVarCSV_trackSip3dValAboveCharm", &OntrackSip3dValAboveCharm_, "OntrackSip3dValAboveCharm_/F");
    addBranch(tree,"OnTagVarCSV_trackSip3dSigAboveCharm", &OntrackSip3dSigAboveCharm_, "OntrackSip3dSigAboveCharm_/F");
    addBranch(tree,"Onn_TagVarCSV_jetNSelectedTracks", &Onn_jetNSelectedTracks_, "Onn_jetNSelectedTracks_/i");
    addBranch(tree,"OnTagVarCSV_jetNSelectedTracks", &OnjetNSelectedTracks_, "OnjetNSelectedTracks_/f");
    addBranch(tree,"OnTagVarCSVTrk_trackPtRel"      , &OntrackPtRel_      , "OntrackPtRel_[Onn_jetNSelectedTracks_]/f"      );
    addBranch(tree,"OnTagVarCSVTrk_trackDeltaR"     , &OntrackDeltaR_     , "OntrackDeltaR_[Onn_jetNSelectedTracks_]/f"     );
    addBranch(tree,"OnTagVarCSVTrk_trackPtRatio"    , &OntrackPtRatio_    , "OntrackPtRatio_[Onn_jetNSelectedTracks_]/f"    );
    addBranch(tree,"OnTagVarCSVTrk_trackSip3dSig"   , &OntrackSip3dSig_   , "OntrackSip3dSig_[Onn_jetNSelectedTracks_]/f"   );
    addBranch(tree,"OnTagVarCSVTrk_trackSip2dSig"   , &OntrackSip2dSig_   , "OntrackSip2dSig_[Onn_jetNSelectedTracks_]/f"   );
    addBranch(tree,"OnTagVarCSVTrk_trackDecayLenVal", &OntrackDecayLenVal_, "OntrackDecayLenVal_[Onn_jetNSelectedTracks_]/f");
    addBranch(tree,"OnTagVarCSVTrk_trackJetDistVal" , &OntrackJetDistVal_ , "OntrackJetDistVal_[Onn_jetNSelectedTracks_]/f" );
    addBranch(tree,"Onn_TagVarCSV_jetNTracksEtaRel", &Onn_jetNTracksEtaRel_, "Onn_jetNTracksEtaRel_/i"                );
    addBranch(tree,"OnTagVarCSV_jetNTracksEtaRel", &OnjetNTracksEtaRel_, "OnjetNTracksEtaRel_/f"                );
    addBranch(tree,"OnTagVarCSV_trackEtaRel"     , &OntrackEtaRel_     , "OntrackEtaRel_[Onn_jetNTracksEtaRel_]/f"  );
    addBranch(tree,"OntrackPParRatio"  , &OntrackPParRatio_  , "OntrackPParRatio_[Onn_jetNSelectedTracks_]/f"  );
    addBranch(tree,"OntrackSip2dVal"   , &OntrackSip2dVal_   , "OntrackSip2dVal_[Onn_jetNSelectedTracks_]/f"   );
    addBranch(tree,"OntrackSip3dVal"   , &OntrackSip3dVal_   , "OntrackSip3dVal_[Onn_jetNSelectedTracks_]/f"   );
    addBranch(tree,"OntrackMomentum"   , &OntrackMomentum_   , "OntrackMomentum_[Onn_jetNSelectedTracks_]/f"   );
    addBranch(tree,"OntrackEta"        , &OntrackEta_        , "OntrackEta_[Onn_jetNSelectedTracks_]/f"        );
    addBranch(tree,"OntrackPPar"       , &OntrackPPar_       , "OntrackPPar_[Onn_jetNSelectedTracks_]/f"       );
    addBranch(tree,"Onn_StoredVertices"    , &Onn_StoredVertices_    , "Onn_StoredVertices_/i"  );
    addBranch(tree,"OnNStoredVertices"    , &OnNStoredVertices_    , "OnNStoredVertices_/f"  );
    addBranch(tree,"OnTagVarCSV_vertexMass"         , &OnvertexMass_         , "OnvertexMass_[Onn_StoredVertices_]/f"         );
    addBranch(tree,"OnTagVarCSV_vertexNTracks"      , &OnvertexNTracks_      , "OnvertexNTracks_[Onn_StoredVertices_]/f"      );
    addBranch(tree,"OnTagVarCSV_vertexEnergyRatio"  , &OnvertexEnergyRatio_  , "OnvertexEnergyRatio_[Onn_StoredVertices_]/f"  );
    addBranch(tree,"OnTagVarCSV_vertexJetDeltaR"    , &OnvertexJetDeltaR_    , "OnvertexJetDeltaR_[Onn_StoredVertices_]/f"    );
    addBranch(tree,"OnTagVarCSV_flightDistance2dVal", &OnflightDistance2dVal_, "OnflightDistance2dVal_[Onn_StoredVertices_]/f");
    addBranch(tree,"OnTagVarCSV_flightDistance2dSig", &OnflightDistance2dSig_, "OnflightDistance2dSig_[Onn_StoredVertices_]/f");
    addBranch(tree,"OnTagVarCSV_flightDistance3dVal", &OnflightDistance3dVal_, "OnflightDistance3dVal_[Onn_StoredVertices_]/f");
    addBranch(tree,"OnTagVarCSV_flightDistance3dSig", &OnflightDistance3dSig_, "OnflightDistance3dSig_[Onn_StoredVertices_]/f");
}

//use either of these functions
bool ntuple_bTagVars1::fillBranches(const reco::ShallowTagInfo tagInfo, const edm::Event& iEvent, std::vector<float> Disc){
 
    reco::TaggingVariableList vars = tagInfo.taggingVariables();

    //*******************
    //
    //  jet general
    //
    //*******************
    lumiBlock_ = iEvent.eventAuxiliary().luminosityBlock();
    runNumber_ = iEvent.eventAuxiliary().run();
    eventNumber_ = iEvent.eventAuxiliary().event();
    OnDeepCSVProbb_          = Disc.at(0);
    OnDeepCSVProbc_          = Disc.at(1);
    OnDeepCSVProbudsg_          = Disc.at(2);
    OnCSVProbb_          = Disc.at(3);
    OnDeepCSVCaloProbb_          = Disc.at(4);
    OnDeepCSVCaloProbc_          = Disc.at(5);
    OnDeepCSVCaloProbudsg_          = Disc.at(6);
    OnCSVCaloProbb_          = Disc.at(7);

    OnJetPt_                 = vars.get(reco::btau::jetPt, -999);
    OntrackJetPt_                 = vars.get(reco::btau::trackJetPt, -999);
    OnjetNSecondaryVertices_      = vars.get(reco::btau::jetNSecondaryVertices, -1);
    OntrackSumJetEtRatio_         = vars.get(reco::btau::trackSumJetEtRatio, -999);
    OntrackSumJetDeltaR_          = vars.get(reco::btau::trackSumJetDeltaR, -999);
    OnvertexCategory_             = vars.get(reco::btau::vertexCategory, -999);
    OntrackSip2dValAboveCharm_    = vars.get(reco::btau::trackSip2dValAboveCharm, -999);
    OntrackSip2dSigAboveCharm_    = vars.get(reco::btau::trackSip2dSigAboveCharm, -999);
    OntrackSip3dValAboveCharm_    = vars.get(reco::btau::trackSip3dValAboveCharm, -999);
    OntrackSip3dSigAboveCharm_    = vars.get(reco::btau::trackSip3dSigAboveCharm, -999);
    OnjetNTracks_ = vars.get(reco::btau::jetNTracks, -1);
    Onn_jetNTracksEtaRel_ = vars.get(reco::btau::jetNTracksEtaRel, -1);
    OnjetNTracksEtaRel_=Onn_jetNTracksEtaRel_;
    Onn_jetNSelectedTracks_ = dump_vector(vars, OntrackMomentum_, reco::btau::trackMomentum,max_OnjetNSelectedTracks_);
    OnjetNSelectedTracks_=Onn_jetNSelectedTracks_;
    
    dump_vector(vars, OntrackEta_, reco::btau::trackEta,max_OnjetNSelectedTracks_)            ;
    dump_vector(vars, OntrackPtRel_, reco::btau::trackPtRel,max_OnjetNSelectedTracks_)        ;
    dump_vector(vars, OntrackPPar_, reco::btau::trackPPar,max_OnjetNSelectedTracks_)          ;
    dump_vector(vars, OntrackDeltaR_, reco::btau::trackDeltaR,max_OnjetNSelectedTracks_)      ;
    dump_vector(vars, OntrackPtRatio_, reco::btau::trackPtRatio,max_OnjetNSelectedTracks_)    ;
    dump_vector(vars, OntrackPParRatio_, reco::btau::trackPParRatio,max_OnjetNSelectedTracks_);
    dump_vector(vars, OntrackSip2dVal_, reco::btau::trackSip2dVal,max_OnjetNSelectedTracks_)  ;
    dump_vector(vars, OntrackSip2dSig_, reco::btau::trackSip2dSig,max_OnjetNSelectedTracks_)  ;
    dump_vector(vars, OntrackSip3dVal_, reco::btau::trackSip3dVal,max_OnjetNSelectedTracks_)  ;
    dump_vector(vars, OntrackSip3dSig_, reco::btau::trackSip3dSig,max_OnjetNSelectedTracks_)  ;
    dump_vector(vars, OntrackDecayLenVal_, reco::btau::trackDecayLenVal,max_OnjetNSelectedTracks_);
    dump_vector(vars, OntrackJetDistVal_, reco::btau::trackJetDistVal,max_OnjetNSelectedTracks_);
    dump_vector(vars, OntrackEtaRel_, reco::btau::trackEtaRel,max_OnjetNSelectedTracks_);

    Onn_StoredVertices_ = dump_vector(vars, OnvertexMass_, reco::btau::vertexMass,max_OnnStoredVertices_);
    OnNStoredVertices_=Onn_StoredVertices_;

    dump_vector(vars, OnvertexNTracks_, reco::btau::vertexNTracks,max_OnnStoredVertices_);
    dump_vector(vars, OnvertexEnergyRatio_, reco::btau::vertexEnergyRatio,max_OnnStoredVertices_);
    dump_vector(vars, OnvertexJetDeltaR_, reco::btau::vertexJetDeltaR,max_OnnStoredVertices_);
    dump_vector(vars, OnflightDistance2dVal_, reco::btau::flightDistance2dVal,max_OnnStoredVertices_);
    dump_vector(vars, OnflightDistance2dSig_, reco::btau::flightDistance2dSig,max_OnnStoredVertices_);
    dump_vector(vars, OnflightDistance3dVal_, reco::btau::flightDistance3dVal,max_OnnStoredVertices_);
    dump_vector(vars, OnflightDistance3dSig_, reco::btau::flightDistance3dSig,max_OnnStoredVertices_);
    return true;
}
