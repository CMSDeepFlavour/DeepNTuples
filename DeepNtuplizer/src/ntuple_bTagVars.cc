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
	//jet general
	tree->Branch("trackJetPt"             , &trackJetPt_             , "trackJetPt_/F"             );
	tree->Branch("jetNTracks"             , &jetNTracks_             , "jetNTracks_/F"             );
	tree->Branch("TagVarCSV_jetNSecondaryVertices"  , &jetNSecondaryVertices_  , "jetNSecondaryVertices_/F"  );
	tree->Branch("TagVarCSV_trackSumJetEtRatio"     , &trackSumJetEtRatio_     , "trackSumJetEtRatio_/F"     );
	tree->Branch("TagVarCSV_trackSumJetDeltaR"      , &trackSumJetDeltaR_      , "trackSumJetDeltaR_/F"      );
	tree->Branch("TagVarCSV_vertexCategory"         , &vertexCategory_         , "vertexCategory_/F"         );
	tree->Branch("TagVarCSV_trackSip2dValAboveCharm", &trackSip2dValAboveCharm_, "trackSip2dValAboveCharm_/F");
	tree->Branch("TagVarCSV_trackSip2dSigAboveCharm", &trackSip2dSigAboveCharm_, "trackSip2dSigAboveCharm_/F");
	tree->Branch("TagVarCSV_trackSip3dValAboveCharm", &trackSip3dValAboveCharm_, "trackSip3dValAboveCharm_/F");
	tree->Branch("TagVarCSV_trackSip3dSigAboveCharm", &trackSip3dSigAboveCharm_, "trackSip3dSigAboveCharm_/F");
	//track info
	tree->Branch("TagVarCSV_jetNSelectedTracks", &jetNSelectedTracks_, "jetNSelectedTracks_/i");
	tree->Branch("TagVarCSVTrk_trackPtRel"      , &trackPtRel_      , "trackPtRel_[jetNSelectedTracks_]/f"      );
	tree->Branch("TagVarCSVTrk_trackDeltaR"     , &trackDeltaR_     , "trackDeltaR_[jetNSelectedTracks_]/f"     );
	tree->Branch("TagVarCSVTrk_trackPtRatio"    , &trackPtRatio_    , "trackPtRatio_[jetNSelectedTracks_]/f"    );
	tree->Branch("TagVarCSVTrk_trackSip3dSig"   , &trackSip3dSig_   , "trackSip3dSig_[jetNSelectedTracks_]/f"   );
	tree->Branch("TagVarCSVTrk_trackSip2dSig"   , &trackSip2dSig_   , "trackSip2dSig_[jetNSelectedTracks_]/f"   );
	tree->Branch("TagVarCSVTrk_trackDecayLenVal", &trackDecayLenVal_, "trackDecayLenVal_[jetNSelectedTracks_]/f");
	tree->Branch("TagVarCSVTrk_trackJetDistVal" , &trackJetDistVal_ , "trackJetDistVal_[jetNSelectedTracks_]/f" );
	tree->Branch("TagVarCSV_jetNTracksEtaRel", &jetNTracksEtaRel_, "jetNTracksEtaRel_/i"                );
	tree->Branch("TagVarCSV_trackEtaRel"     , &trackEtaRel_     , "trackEtaRel_[jetNTracksEtaRel_]/f"  );
	tree->Branch("trackPParRatio"  , &trackPParRatio_  , "trackPParRatio_[jetNSelectedTracks_]/f"  );
	tree->Branch("trackSip2dVal"   , &trackSip2dVal_   , "trackSip2dVal_[jetNSelectedTracks_]/f"   );
	tree->Branch("trackSip3dVal"   , &trackSip3dVal_   , "trackSip3dVal_[jetNSelectedTracks_]/f"   );
	tree->Branch("trackMomentum"   , &trackMomentum_   , "trackMomentum_[jetNSelectedTracks_]/f"   );
	tree->Branch("trackEta"        , &trackEta_        , "trackEta_[jetNSelectedTracks_]/f"        );
	tree->Branch("trackPPar"       , &trackPPar_       , "trackPPar_[jetNSelectedTracks_]/f"       );
	//SV info
	tree->Branch("n_StoredVertices"    , &nStoredVertices_    , "nStoredVertices_/i"  );
	tree->Branch("TagVarCSV_vertexMass"         , &vertexMass_         , "vertexMass_[nStoredVertices_]/f"         );
	tree->Branch("TagVarCSV_vertexNTracks"      , &vertexNTracks_      , "vertexNTracks_[nStoredVertices_]/f"      );
	tree->Branch("TagVarCSV_vertexEnergyRatio"  , &vertexEnergyRatio_  , "vertexEnergyRatio_[nStoredVertices_]/f"  );
	tree->Branch("TagVarCSV_vertexJetDeltaR"    , &vertexJetDeltaR_    , "vertexJetDeltaR_[nStoredVertices_]/f"    );
	tree->Branch("TagVarCSV_flightDistance2dVal", &flightDistance2dVal_, "flightDistance2dVal_[nStoredVertices_]/f");
	tree->Branch("TagVarCSV_flightDistance2dSig", &flightDistance2dSig_, "flightDistance2dSig_[nStoredVertices_]/f");
	tree->Branch("TagVarCSV_flightDistance3dVal", &flightDistance3dVal_, "flightDistance3dVal_[nStoredVertices_]/f");
	tree->Branch("TagVarCSV_flightDistance3dSig", &flightDistance3dSig_, "flightDistance3dSig_[nStoredVertices_]/f");	
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
	reco::TaggingVariableList vars = tagInfo->taggingVariables();
        
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
	jetNTracksEtaRel_ = vars.get(reco::btau::jetNTracksEtaRel, -1);
        
	//*******************
	//
	//  track info  //FIXME: check that the order is the same of the charged components!
	//  FIXME: right now there are no default values in the vectors... Is this something we want? or do we want fixed size vectors with zero padding in the end?
	//
	//*******************
	jetNSelectedTracks_ = dump_vector(vars, trackMomentum_, reco::btau::trackMomentum);
	dump_vector(vars, trackEta_, reco::btau::trackEta)            ;
	dump_vector(vars, trackPtRel_, reco::btau::trackPtRel)        ;
	dump_vector(vars, trackPPar_, reco::btau::trackPPar)          ;
	dump_vector(vars, trackDeltaR_, reco::btau::trackDeltaR)      ;
	dump_vector(vars, trackPtRatio_, reco::btau::trackPtRatio)    ;
	dump_vector(vars, trackPParRatio_, reco::btau::trackPParRatio);
	dump_vector(vars, trackSip2dVal_, reco::btau::trackSip2dVal)  ;
	dump_vector(vars, trackSip2dSig_, reco::btau::trackSip2dSig)  ;
	dump_vector(vars, trackSip3dVal_, reco::btau::trackSip3dVal)  ;
	dump_vector(vars, trackSip3dSig_, reco::btau::trackSip3dSig)  ;
	dump_vector(vars, trackDecayLenVal_, reco::btau::trackDecayLenVal);
	dump_vector(vars, trackJetDistVal_, reco::btau::trackJetDistVal);
	dump_vector(vars, trackEtaRel_, reco::btau::trackEtaRel);

	//*******************
	//
	//  vertex info
	//
	//*******************
	nStoredVertices_ = dump_vector(vars, vertexMass_, reco::btau::vertexMass);
	dump_vector(vars, vertexNTracks_, reco::btau::vertexNTracks);
	dump_vector(vars, vertexEnergyRatio_, reco::btau::vertexEnergyRatio);
	dump_vector(vars, vertexJetDeltaR_, reco::btau::vertexJetDeltaR);
	dump_vector(vars, flightDistance2dVal_, reco::btau::flightDistance2dVal);
	dump_vector(vars, flightDistance2dSig_, reco::btau::flightDistance2dSig);
	dump_vector(vars, flightDistance3dVal_, reco::btau::flightDistance3dVal);
	dump_vector(vars, flightDistance3dSig_, reco::btau::flightDistance3dSig);

	return true;
}
