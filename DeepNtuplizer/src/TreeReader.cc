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
    addBranch(tree, (On + "Jet_pt").c_str()            , &Jet_pt_);
    addBranch(tree, (On + "Jet_eta").c_str()            , &Jet_eta_);
    addBranch(tree, (On + "nJet").c_str()            , &nJet_);
    addBranch(tree, (On + "Jet_nFirstSV").c_str()            , &Jet_nFirstSV_);
    addBranch(tree, (On + "Jet_nFirstTrkTagVarCSV").c_str()            , &Jet_nFirstTrkTagVarCSV_ );
    addBranch(tree, (On + "Jet_nLastTrkTagVarCSV").c_str()            , &Jet_nLastTrkTagVarCSV_  );
    addBranch(tree, (On + "Jet_nFirstTrkEtaRelTagVarCSV").c_str()            , &Jet_nFirstTrkEtaRelTagVarCSV_ );
    addBranch(tree, (On + "Jet_nLastTrkEtaRelTagVarCSV").c_str()            , &Jet_nLastTrkEtaRelTagVarCSV_ );
    addBranch(tree, (On + "TagVarCSV_trackJetPt").c_str()            , &trackJetPt_);
    addBranch(tree, (On + "Jet_ntracks").c_str()            , &jetNTracks_ );
    addBranch(tree, (On + "TagVarCSV_jetNSecondaryVertices").c_str() , &jetNSecondaryVertices_);
    addBranch(tree, (On + "TagVarCSV_trackSumJetEtRatio").c_str()    , &trackSumJetEtRatio_);
    addBranch(tree, (On + "TagVarCSV_trackSumJetDeltaR").c_str()     , &trackSumJetDeltaR_);
    addBranch(tree, (On + "TagVarCSV_vertexCategory").c_str()        , &vertexCategory_ );
    addBranch(tree, (On + "TagVarCSV_trackSip2dValAboveCharm").c_str(), &trackSip2dValAboveCharm_);
    addBranch(tree, (On + "TagVarCSV_trackSip2dSigAboveCharm").c_str(), &trackSip2dSigAboveCharm_);
    addBranch(tree, (On + "TagVarCSV_trackSip3dValAboveCharm").c_str(), &trackSip3dValAboveCharm_);
    addBranch(tree, (On + "TagVarCSV_trackSip3dSigAboveCharm").c_str(), &trackSip3dSigAboveCharm_);
    //track info
    addBranch(tree, (On + "TagVarCSV_jetNTracks").c_str(), &jetNSelectedTracks_);

    addBranch(tree, (On + "TagVarCSV_trackPtRel").c_str()     , &trackPtRel_      );
    addBranch(tree, (On + "TagVarCSV_trackDeltaR").c_str()    , &trackDeltaR_     );
    addBranch(tree, (On + "TagVarCSV_trackPtRatio").c_str()   , &trackPtRatio_    );
    addBranch(tree, (On + "TagVarCSV_trackSip3dSig").c_str()  , &trackSip3dSig_   );
    addBranch(tree, (On + "TagVarCSV_trackSip2dSig").c_str()  , &trackSip2dSig_   );
    addBranch(tree, (On + "TagVarCSV_trackDecayLenVal").c_str(), &trackDecayLenVal_ );
    addBranch(tree, (On + "TagVarCSV_trackJetDistVal").c_str(), &trackJetDistVal_ );

    addBranch(tree, (On + "TagVarCSV_jetNTracksEtaRel").c_str(), &jetNTracksEtaRel_ );

    addBranch(tree, (On + "TagVarCSV_trackEtaRel").c_str()    , &trackEtaRel_      );

    addBranch(tree, (On + "TagVarCSV_trackPParRatio").c_str() , &trackPParRatio_   );
    addBranch(tree, (On + "TagVarCSV_trackSip2dVal").c_str()  , &trackSip2dVal_    );
    addBranch(tree, (On + "TagVarCSV_trackSip3dVal").c_str()  , &trackSip3dVal_    );
    addBranch(tree, (On + "TagVarCSV_trackMomentum").c_str()  , &trackMomentum_    );
    addBranch(tree, (On + "TagVarCSV_trackEta").c_str()       , &trackEta_         );
    addBranch(tree, (On + "TagVarCSV_trackPPar").c_str()      , &trackPPar_        );
    //SV info
    addBranch(tree, (On + "TagVarCSV_vertexMass").c_str()        , &vertexMass_          );
    addBranch(tree, (On + "TagVarCSV_vertexNTracks").c_str()     , &vertexNTracks_       );
    addBranch(tree, (On + "TagVarCSV_vertexEnergyRatio").c_str() , &vertexEnergyRatio_   );
    addBranch(tree, (On + "TagVarCSV_vertexJetDeltaR").c_str()   , &vertexJetDeltaR_     );
    addBranch(tree, (On + "TagVarCSV_flightDistance2dVal").c_str(), &flightDistance2dVal_);
    addBranch(tree, (On + "TagVarCSV_flightDistance2dSig").c_str(), &flightDistance2dSig_);
    addBranch(tree, (On + "TagVarCSV_flightDistance3dVal").c_str(), &flightDistance3dVal_);
    addBranch(tree, (On + "TagVarCSV_flightDistance3dSig").c_str(), &flightDistance3dSig_);


    if(MCtest){
      addBranch(tree, "Jet_hadronFlavour", &Jet_hadronFlavour);
    }

    addBranch(tree, "Jet_DeepCSVBDisc", &Jet_DeepCSVBDisc);
    addBranch(tree, "Jet_DeepCSVBDiscN", &Jet_DeepCSVBDiscN);
    addBranch(tree, "Jet_DeepCSVBDiscP", &Jet_DeepCSVBDiscP);
    addBranch(tree, "Jet_DeepCSVCvsLDisc", &Jet_DeepCSVCvsLDisc);
    addBranch(tree, "Jet_DeepCSVCvsLDiscN", &Jet_DeepCSVCvsLDiscN);
    addBranch(tree, "Jet_DeepCSVCvsLDiscP", &Jet_DeepCSVCvsLDiscP);
    addBranch(tree, "Jet_DeepCSVCvsBDisc", &Jet_DeepCSVCvsBDisc);
    addBranch(tree, "Jet_DeepCSVCvsBDiscN", &Jet_DeepCSVCvsBDiscN);
    addBranch(tree, "Jet_DeepCSVCvsBDiscP", &Jet_DeepCSVCvsBDiscP);
    addBranch(tree, "Jet_DeepCSVb", &Jet_DeepCSVb);
    addBranch(tree, "Jet_DeepCSVc", &Jet_DeepCSVc);
    addBranch(tree, "Jet_DeepCSVl", &Jet_DeepCSVl);
    addBranch(tree, "Jet_DeepCSVbb", &Jet_DeepCSVbb);
    addBranch(tree, "Jet_DeepCSVcc", &Jet_DeepCSVcc);
    addBranch(tree, "Jet_DeepCSVbN", &Jet_DeepCSVbN);
    addBranch(tree, "Jet_DeepCSVcN", &Jet_DeepCSVcN);
    addBranch(tree, "Jet_DeepCSVlN", &Jet_DeepCSVlN);
    addBranch(tree, "Jet_DeepCSVbbN", &Jet_DeepCSVbbN);
    addBranch(tree, "Jet_DeepCSVccN", &Jet_DeepCSVccN);
    addBranch(tree, "Jet_DeepCSVbP", &Jet_DeepCSVbP);
    addBranch(tree, "Jet_DeepCSVcP", &Jet_DeepCSVcP);
    addBranch(tree, "Jet_DeepCSVlP", &Jet_DeepCSVlP);
    addBranch(tree, "Jet_DeepCSVbbP", &Jet_DeepCSVbbP);
    addBranch(tree, "Jet_DeepCSVccP", &Jet_DeepCSVccP);

    addBranch(tree, "nBitTrigger", &nBitTrigger);
    addBranch(tree, "BitTrigger", &BitTrigger);
    addBranch(tree, "nPV", &nPV);
    addBranch(tree, "pthat", &pthat);


}

