/*
 * TreeReader.cc
 *
 *      Author: ebols
 */

#include "../interface/TreeReader.h"
#include <sstream>
#include <algorithm>

using namespace std;



void TreeReader::initBranches(TTree* tree, bool MCtest){
  initBranches(tree,"", MCtest);
}


void TreeReader::initBranches(TTree* tree, string On, bool MCtest){
    //jet general
    addBranch(tree, (On + "Jet_pt").c_str()            , &Jet_pt_             , (On + "Jet_pt_/F").c_str()            );
    addBranch(tree, (On + "Jet_eta").c_str()            , &Jet_eta_             , (On + "Jet_eta_/F").c_str()            );
    addBranch(tree, (On + "nJet").c_str()            , &nJet_             , (On + "jetNTracks_/F").c_str()            );
    addBranch(tree, (On + "Jet_nFirstSV").c_str()            , &Jet_nFirstSV_             , (On + "jetNTracks_/F").c_str()            );
    addBranch(tree, (On + "Jet_nFirstTrkTagVarCSV").c_str()            , &Jet_nFirstTrkTagVarCSV_             , (On + "trackJetPt_/F").c_str() );
    addBranch(tree, (On + "Jet_nLastTrkTagVarCSV").c_str()            , &Jet_nLastTrkTagVarCSV_             , (On + "trackJetPt_/F").c_str() );
    addBranch(tree, (On + "Jet_nFirstTrkEtaRelTagVarCSV").c_str()            , &Jet_nFirstTrkEtaRelTagVarCSV_             , (On + "trackJetPt_/F").c_str() );
    addBranch(tree, (On + "Jet_nLastTrkEtaRelTagVarCSV").c_str()            , &Jet_nLastTrkEtaRelTagVarCSV_             , (On + "trackJetPt_/F").c_str() );
    addBranch(tree, (On + "TagVarCSV_trackJetPt").c_str()            , &trackJetPt_             , (On + "trackJetPt_/F").c_str()            );
    addBranch(tree, (On + "Jet_ntracks").c_str()            , &jetNTracks_             , (On + "jetNTracks_/F").c_str()            );
    addBranch(tree, (On + "TagVarCSV_jetNSecondaryVertices").c_str() , &jetNSecondaryVertices_  , (On + "jetNSecondaryVertices_/F").c_str() );
    addBranch(tree, (On + "TagVarCSV_trackSumJetEtRatio").c_str()    , &trackSumJetEtRatio_     , (On + "trackSumJetEtRatio_/F").c_str()    );
    addBranch(tree, (On + "TagVarCSV_trackSumJetDeltaR").c_str()     , &trackSumJetDeltaR_      , (On + "trackSumJetDeltaR_/F").c_str()     );
    addBranch(tree, (On + "TagVarCSV_vertexCategory").c_str()        , &vertexCategory_         , (On + "vertexCategory_/F").c_str()        );
    addBranch(tree, (On + "TagVarCSV_trackSip2dValAboveCharm").c_str(), &trackSip2dValAboveCharm_, (On + "trackSip2dValAboveCharm_/F").c_str());
    addBranch(tree, (On + "TagVarCSV_trackSip2dSigAboveCharm").c_str(), &trackSip2dSigAboveCharm_, (On + "trackSip2dSigAboveCharm_/F").c_str());
    addBranch(tree, (On + "TagVarCSV_trackSip3dValAboveCharm").c_str(), &trackSip3dValAboveCharm_, (On + "trackSip3dValAboveCharm_/F").c_str());
    addBranch(tree, (On + "TagVarCSV_trackSip3dSigAboveCharm").c_str(), &trackSip3dSigAboveCharm_, (On + "trackSip3dSigAboveCharm_/F").c_str());
    //track info
    addBranch(tree, (On + "TagVarCSV_jetNTracks").c_str(), &jetNSelectedTracks_, (On + "jetNSelectedTracks_/f").c_str());

    addBranch(tree, (On + "TagVarCSV_trackPtRel").c_str()     , &trackPtRel_      , (On + "trackPtRel_[n_jetNSelectedTracks_]/f").c_str()     );
    addBranch(tree, (On + "TagVarCSV_trackDeltaR").c_str()    , &trackDeltaR_     , (On + "trackDeltaR_[n_jetNSelectedTracks_]/f").c_str()    );
    addBranch(tree, (On + "TagVarCSV_trackPtRatio").c_str()   , &trackPtRatio_    , (On + "trackPtRatio_[n_jetNSelectedTracks_]/f").c_str()   );
    addBranch(tree, (On + "TagVarCSV_trackSip3dSig").c_str()  , &trackSip3dSig_   , (On + "trackSip3dSig_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "TagVarCSV_trackSip2dSig").c_str()  , &trackSip2dSig_   , (On + "trackSip2dSig_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "TagVarCSV_trackDecayLenVal").c_str(), &trackDecayLenVal_, (On + "trackDecayLenVal_[n_jetNSelectedTracks_]/f").c_str());
    addBranch(tree, (On + "TagVarCSV_trackJetDistVal").c_str(), &trackJetDistVal_ , (On + "trackJetDistVal_[n_jetNSelectedTracks_]/f").c_str());

    addBranch(tree, (On + "TagVarCSV_jetNTracksEtaRel").c_str(), &jetNTracksEtaRel_, (On + "jetNTracksEtaRel_/f").c_str()               );

    addBranch(tree, (On + "TagVarCSV_trackEtaRel").c_str()    , &trackEtaRel_     , (On + "trackEtaRel_[n_jetNTracksEtaRel_]/f").c_str() );

    addBranch(tree, (On + "TagVarCSV_trackPParRatio").c_str() , &trackPParRatio_  , (On + "trackPParRatio_[n_jetNSelectedTracks_]/f").c_str() );
    addBranch(tree, (On + "TagVarCSV_trackSip2dVal").c_str()  , &trackSip2dVal_   , (On + "trackSip2dVal_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "TagVarCSV_trackSip3dVal").c_str()  , &trackSip3dVal_   , (On + "trackSip3dVal_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "TagVarCSV_trackMomentum").c_str()  , &trackMomentum_   , (On + "trackMomentum_[n_jetNSelectedTracks_]/f").c_str()  );
    addBranch(tree, (On + "TagVarCSV_trackEta").c_str()       , &trackEta_        , (On + "trackEta_[n_jetNSelectedTracks_]/f").c_str()       );
    addBranch(tree, (On + "TagVarCSV_trackPPar").c_str()      , &trackPPar_       , (On + "trackPPar_[n_jetNSelectedTracks_]/f").c_str()      );
    //SV info
    addBranch(tree, (On + "TagVarCSV_vertexMass").c_str()        , &vertexMass_         , (On + "vertexMass_[n_StoredVertices_]/f").c_str()        );
    addBranch(tree, (On + "TagVarCSV_vertexNTracks").c_str()     , &vertexNTracks_      , (On + "vertexNTracks_[n_StoredVertices_]/f").c_str()     );
    addBranch(tree, (On + "TagVarCSV_vertexEnergyRatio").c_str() , &vertexEnergyRatio_  , (On + "vertexEnergyRatio_[n_StoredVertices_]/f").c_str() );
    addBranch(tree, (On + "TagVarCSV_vertexJetDeltaR").c_str()   , &vertexJetDeltaR_    , (On + "vertexJetDeltaR_[n_StoredVertices_]/f").c_str()   );
    addBranch(tree, (On + "TagVarCSV_flightDistance2dVal").c_str(), &flightDistance2dVal_, (On + "flightDistance2dVal_[n_StoredVertices_]/f").c_str());
    addBranch(tree, (On + "TagVarCSV_flightDistance2dSig").c_str(), &flightDistance2dSig_, (On + "flightDistance2dSig_[n_StoredVertices_]/f").c_str());
    addBranch(tree, (On + "TagVarCSV_flightDistance3dVal").c_str(), &flightDistance3dVal_, (On + "flightDistance3dVal_[n_StoredVertices_]/f").c_str());
    addBranch(tree, (On + "TagVarCSV_flightDistance3dSig").c_str(), &flightDistance3dSig_, (On + "flightDistance3dSig_[n_StoredVertices_]/f").c_str());


    if(MCtest){
      addBranch(tree, "Jet_hadronFlavour", &Jet_hadronFlavour, "Jet_hadronFlavour");
    }

    addBranch(tree, "Jet_DeepCSVBDisc", &Jet_DeepCSVBDisc, "&b_Jet_DeepCSVBDisc");
    addBranch(tree, "Jet_DeepCSVBDiscN", &Jet_DeepCSVBDiscN, "&b_Jet_DeepCSVBDiscN");
    addBranch(tree, "Jet_DeepCSVBDiscP", &Jet_DeepCSVBDiscP, "&b_Jet_DeepCSVBDiscP");
    addBranch(tree, "Jet_DeepCSVCvsLDisc", &Jet_DeepCSVCvsLDisc, "&b_Jet_DeepCSVCvsLDisc");
    addBranch(tree, "Jet_DeepCSVCvsLDiscN", &Jet_DeepCSVCvsLDiscN, "&b_Jet_DeepCSVCvsLDiscN");
    addBranch(tree, "Jet_DeepCSVCvsLDiscP", &Jet_DeepCSVCvsLDiscP, "&b_Jet_DeepCSVCvsLDiscP");
    addBranch(tree, "Jet_DeepCSVCvsBDisc", &Jet_DeepCSVCvsBDisc, "&b_Jet_DeepCSVCvsBDisc");
    addBranch(tree, "Jet_DeepCSVCvsBDiscN", &Jet_DeepCSVCvsBDiscN, "&b_Jet_DeepCSVCvsBDiscN");
    addBranch(tree, "Jet_DeepCSVCvsBDiscP", &Jet_DeepCSVCvsBDiscP, "&b_Jet_DeepCSVCvsBDiscP");
    addBranch(tree, "Jet_DeepCSVb", &Jet_DeepCSVb, "&b_Jet_DeepCSVb");
    addBranch(tree, "Jet_DeepCSVc", &Jet_DeepCSVc, "&b_Jet_DeepCSVc");
    addBranch(tree, "Jet_DeepCSVl", &Jet_DeepCSVl, "&b_Jet_DeepCSVl");
    addBranch(tree, "Jet_DeepCSVbb", &Jet_DeepCSVbb, "&b_Jet_DeepCSVbb");
    addBranch(tree, "Jet_DeepCSVcc", &Jet_DeepCSVcc, "&b_Jet_DeepCSVcc");
    addBranch(tree, "Jet_DeepCSVbN", &Jet_DeepCSVbN, "&b_Jet_DeepCSVbN");
    addBranch(tree, "Jet_DeepCSVcN", &Jet_DeepCSVcN, "&b_Jet_DeepCSVcN");
    addBranch(tree, "Jet_DeepCSVlN", &Jet_DeepCSVlN, "&b_Jet_DeepCSVlN");
    addBranch(tree, "Jet_DeepCSVbbN", &Jet_DeepCSVbbN, "&b_Jet_DeepCSVbbN");
    addBranch(tree, "Jet_DeepCSVccN", &Jet_DeepCSVccN, "&b_Jet_DeepCSVccN");
    addBranch(tree, "Jet_DeepCSVbP", &Jet_DeepCSVbP, "&b_Jet_DeepCSVbP");
    addBranch(tree, "Jet_DeepCSVcP", &Jet_DeepCSVcP, "&b_Jet_DeepCSVcP");
    addBranch(tree, "Jet_DeepCSVlP", &Jet_DeepCSVlP, "&b_Jet_DeepCSVlP");
    addBranch(tree, "Jet_DeepCSVbbP", &Jet_DeepCSVbbP, "&b_Jet_DeepCSVbbP");
    addBranch(tree, "Jet_DeepCSVccP", &Jet_DeepCSVccP, "&b_Jet_DeepCSVccP");

    addBranch(tree, "nBitTrigger", &nBitTrigger, "&b_nBitTrigger");
    addBranch(tree, "BitTrigger", &BitTrigger, "&b_BitTrigger");
    addBranch(tree, "nPV", &nPV, "&b_nPV");
    addBranch(tree, "pthat", &pthat, "&b_pthat");


}

