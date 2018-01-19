/*
 * ntuple_bTagVars.cc
 *
 *      Author: mverzett
 */

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
#include "../interface/ntuple_bTagVars1.h"


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
bool ntuple_bTagVars1::fillBranches(const reco::ShallowTagInfo tagInfo, const reco::ShallowTagInfo tagInfo1 , const edm::Event& iEvent, std::vector<float> Disc, std::vector<float> Disc1){
 
    reco::TaggingVariableList vars = tagInfo.taggingVariables();

    reco::TaggingVariableList vars1 = tagInfo1.taggingVariables();
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

    Onn_jetNSelectedTracks_ = dump_vector(vars, OntrackMomentum_, reco::btau::trackMomentum,100);
    OnjetNSelectedTracks_=Onn_jetNSelectedTracks_;
    
    dump_vector(vars, OntrackEta_, reco::btau::trackEta,100)            ;
    dump_vector(vars, OntrackPtRel_, reco::btau::trackPtRel,100)        ;
    dump_vector(vars, OntrackPPar_, reco::btau::trackPPar,100)          ;
    dump_vector(vars, OntrackDeltaR_, reco::btau::trackDeltaR,100)      ;
    dump_vector(vars, OntrackPtRatio_, reco::btau::trackPtRatio,100)    ;
    dump_vector(vars, OntrackPParRatio_, reco::btau::trackPParRatio,100);
    dump_vector(vars, OntrackSip2dVal_, reco::btau::trackSip2dVal,100)  ;
    dump_vector(vars, OntrackSip2dSig_, reco::btau::trackSip2dSig,100)  ;
    dump_vector(vars, OntrackSip3dVal_, reco::btau::trackSip3dVal,100)  ;
    dump_vector(vars, OntrackSip3dSig_, reco::btau::trackSip3dSig,100)  ;
    dump_vector(vars, OntrackDecayLenVal_, reco::btau::trackDecayLenVal,100);
    dump_vector(vars, OntrackJetDistVal_, reco::btau::trackJetDistVal,100);
    dump_vector(vars, OntrackEtaRel_, reco::btau::trackEtaRel,100);

    Onn_StoredVertices_ = dump_vector(vars, OnvertexMass_, reco::btau::vertexMass,10);
    OnNStoredVertices_=Onn_StoredVertices_;

    dump_vector(vars, OnvertexNTracks_, reco::btau::vertexNTracks,10);
    dump_vector(vars, OnvertexEnergyRatio_, reco::btau::vertexEnergyRatio,10);
    dump_vector(vars, OnvertexJetDeltaR_, reco::btau::vertexJetDeltaR,10);
    dump_vector(vars, OnflightDistance2dVal_, reco::btau::flightDistance2dVal,10);
    dump_vector(vars, OnflightDistance2dSig_, reco::btau::flightDistance2dSig,10);
    dump_vector(vars, OnflightDistance3dVal_, reco::btau::flightDistance3dVal,10);
    dump_vector(vars, OnflightDistance3dSig_, reco::btau::flightDistance3dSig,10);

    DeepCSVProbb_ = Disc1.at(0);
    DeepCSVProbc_ = Disc1.at(1);
    DeepCSVProbudsg_ = Disc1.at(2);
    CSVProbb_ = Disc1.at(3);
    JetPt_                 = vars1.get(reco::btau::jetPt, -999);
    trackJetPt_                 = vars1.get(reco::btau::trackJetPt, -999);
    jetNSecondaryVertices_      = vars1.get(reco::btau::jetNSecondaryVertices, -1);
    trackSumJetEtRatio_         = vars1.get(reco::btau::trackSumJetEtRatio, -999);
    trackSumJetDeltaR_          = vars1.get(reco::btau::trackSumJetDeltaR, -999);
    vertexCategory_             = vars1.get(reco::btau::vertexCategory, -999);
    trackSip2dValAboveCharm_    = vars1.get(reco::btau::trackSip2dValAboveCharm, -999);
    trackSip2dSigAboveCharm_    = vars1.get(reco::btau::trackSip2dSigAboveCharm, -999);
    trackSip3dValAboveCharm_    = vars1.get(reco::btau::trackSip3dValAboveCharm, -999);
    trackSip3dSigAboveCharm_    = vars1.get(reco::btau::trackSip3dSigAboveCharm, -999);

    jetNTracks_ = vars1.get(reco::btau::jetNTracks, -1);

    n_jetNTracksEtaRel_ = vars1.get(reco::btau::jetNTracksEtaRel, -1);

    jetNTracksEtaRel_=n_jetNTracksEtaRel_;
    
    n_jetNSelectedTracks_ = dump_vector(vars1, trackMomentum_, reco::btau::trackMomentum,100);

    jetNSelectedTracks_=n_jetNSelectedTracks_;
    
    dump_vector(vars1, trackEta_, reco::btau::trackEta,100);
    dump_vector(vars1, trackPtRel_, reco::btau::trackPtRel,100);
    dump_vector(vars1, trackPPar_, reco::btau::trackPPar,100);
    dump_vector(vars1, trackDeltaR_, reco::btau::trackDeltaR,100);
    dump_vector(vars1, trackPtRatio_, reco::btau::trackPtRatio,100);
    dump_vector(vars1, trackPParRatio_, reco::btau::trackPParRatio,100);
    dump_vector(vars1, trackSip2dVal_, reco::btau::trackSip2dVal,100);
    dump_vector(vars1, trackSip2dSig_, reco::btau::trackSip2dSig,100);
    dump_vector(vars1, trackSip3dVal_, reco::btau::trackSip3dVal,100);
    dump_vector(vars1, trackSip3dSig_, reco::btau::trackSip3dSig,100);
    dump_vector(vars1, trackDecayLenVal_, reco::btau::trackDecayLenVal,100);
    dump_vector(vars1, trackJetDistVal_, reco::btau::trackJetDistVal,100);
    dump_vector(vars1, trackEtaRel_, reco::btau::trackEtaRel,100);

    n_StoredVertices_ = dump_vector(vars1, vertexMass_, reco::btau::vertexMass,10);
    NStoredVertices_=n_StoredVertices_;

    dump_vector(vars1, vertexNTracks_, reco::btau::vertexNTracks,10);
    dump_vector(vars1, vertexEnergyRatio_, reco::btau::vertexEnergyRatio,10);
    dump_vector(vars1, vertexJetDeltaR_, reco::btau::vertexJetDeltaR,10);
    dump_vector(vars1, flightDistance2dVal_, reco::btau::flightDistance2dVal,10);
    dump_vector(vars1, flightDistance2dSig_, reco::btau::flightDistance2dSig,10);
    dump_vector(vars1, flightDistance3dVal_, reco::btau::flightDistance3dVal,10);
    dump_vector(vars1, flightDistance3dSig_, reco::btau::flightDistance3dSig,10);


    return true;
}
