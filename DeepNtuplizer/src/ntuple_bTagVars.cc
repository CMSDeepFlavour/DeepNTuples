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
	tree->Branch("jetNTracksEtaRel"       , &jetNTracksEtaRel_       , "jetNTracksEtaRel_/F"       );
	tree->Branch("jetNSecondaryVertices"  , &jetNSecondaryVertices_  , "jetNSecondaryVertices_/F"  );
	tree->Branch("trackSumJetEtRatio"     , &trackSumJetEtRatio_     , "trackSumJetEtRatio_/F"     );
	tree->Branch("trackSumJetDeltaR"      , &trackSumJetDeltaR_      , "trackSumJetDeltaR_/F"      );
	tree->Branch("vertexCategory"         , &vertexCategory_         , "vertexCategory_/F"         );
	tree->Branch("trackSip2dValAboveCharm", &trackSip2dValAboveCharm_, "trackSip2dValAboveCharm_/F");
	tree->Branch("trackSip2dSigAboveCharm", &trackSip2dSigAboveCharm_, "trackSip2dSigAboveCharm_/F");
	tree->Branch("trackSip3dValAboveCharm", &trackSip3dValAboveCharm_, "trackSip3dValAboveCharm_/F");
	tree->Branch("trackSip3dSigAboveCharm", &trackSip3dSigAboveCharm_, "trackSip3dSigAboveCharm_/F");
	//track info
	tree->Branch("nStoredTracks"   , &nStoredTracks_   , "nStoredTracks_/i");
	tree->Branch("trackMomentum"   , &trackMomentum_   , "trackMomentum_[nStoredTracks_]/f"   );
	tree->Branch("trackEta"        , &trackEta_        , "trackEta_[nStoredTracks_]/f"        );
	tree->Branch("trackPtRel"      , &trackPtRel_      , "trackPtRel_[nStoredTracks_]/f"      );
	tree->Branch("trackPPar"       , &trackPPar_       , "trackPPar_[nStoredTracks_]/f"       );
	tree->Branch("trackDeltaR"     , &trackDeltaR_     , "trackDeltaR_[nStoredTracks_]/f"     );
	tree->Branch("trackPtRatio"    , &trackPtRatio_    , "trackPtRatio_[nStoredTracks_]/f"    );
	tree->Branch("trackPParRatio"  , &trackPParRatio_  , "trackPParRatio_[nStoredTracks_]/f"  );
	tree->Branch("trackSip2dVal"   , &trackSip2dVal_   , "trackSip2dVal_[nStoredTracks_]/f"   );
	tree->Branch("trackSip2dSig"   , &trackSip2dSig_   , "trackSip2dSig_[nStoredTracks_]/f"   );
	tree->Branch("trackSip3dVal"   , &trackSip3dVal_   , "trackSip3dVal_[nStoredTracks_]/f"   );
	tree->Branch("trackSip3dSig"   , &trackSip3dSig_   , "trackSip3dSig_[nStoredTracks_]/f"   );
	tree->Branch("trackDecayLenVal", &trackDecayLenVal_, "trackDecayLenVal_[nStoredTracks_]/f");
	tree->Branch("trackJetDistVal" , &trackJetDistVal_ , "trackJetDistVal_[nStoredTracks_]/f" );
	tree->Branch("trackEtaRel"     , &trackEtaRel_     , "trackEtaRel_[nStoredTracks_]/f"     );
	//SV info
	tree->Branch("nStoredVertices"    , &nStoredVertices_    , "nStoredVertices_/i"  );
	tree->Branch("vertexMass"         , &vertexMass_         , "vertexMass_[nStoredVertices_]/f"         );
	tree->Branch("vertexNTracks"      , &vertexNTracks_      , "vertexNTracks_[nStoredVertices_]/f"      );
	tree->Branch("vertexEnergyRatio"  , &vertexEnergyRatio_  , "vertexEnergyRatio_[nStoredVertices_]/f"  );
	tree->Branch("vertexJetDeltaR"    , &vertexJetDeltaR_    , "vertexJetDeltaR_[nStoredVertices_]/f"    );
	tree->Branch("flightDistance2dVal", &flightDistance2dVal_, "flightDistance2dVal_[nStoredVertices_]/f");
	tree->Branch("flightDistance2dSig", &flightDistance2dSig_, "flightDistance2dSig_[nStoredVertices_]/f");
	tree->Branch("flightDistance3dVal", &flightDistance3dVal_, "flightDistance3dVal_[nStoredVertices_]/f");
	tree->Branch("flightDistance3dSig", &flightDistance3dSig_, "flightDistance3dSig_[nStoredVertices_]/f");	
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
	nStoredTracks_ = 100;
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackMomentum_, reco::btau::trackMomentum)  );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackEta_, reco::btau::trackEta)            );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackPtRel_, reco::btau::trackPtRel)        );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackPPar_, reco::btau::trackPPar)          );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackDeltaR_, reco::btau::trackDeltaR)      );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackPtRatio_, reco::btau::trackPtRatio)    );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackPParRatio_, reco::btau::trackPParRatio));
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackSip2dVal_, reco::btau::trackSip2dVal)  );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackSip2dSig_, reco::btau::trackSip2dSig)  );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackSip3dVal_, reco::btau::trackSip3dVal)  );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackSip3dSig_, reco::btau::trackSip3dSig)  );
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackDecayLenVal_, reco::btau::trackDecayLenVal));
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackJetDistVal_, reco::btau::trackJetDistVal));
	nStoredTracks_ = min(nStoredTracks_ , dump_vector(vars, trackEtaRel_, reco::btau::trackEtaRel));

	//*******************
	//
	//  vertex info
	//
	//*******************
	dump_vector(vars, vertexMass_, reco::btau::vertexMass);
	dump_vector(vars, vertexNTracks_, reco::btau::vertexNTracks);
	dump_vector(vars, vertexEnergyRatio_, reco::btau::vertexEnergyRatio);
	dump_vector(vars, vertexJetDeltaR_, reco::btau::vertexJetDeltaR);
	dump_vector(vars, flightDistance2dVal_, reco::btau::flightDistance2dVal);
	dump_vector(vars, flightDistance2dSig_, reco::btau::flightDistance2dSig);
	dump_vector(vars, flightDistance3dVal_, reco::btau::flightDistance3dVal);
	dump_vector(vars, flightDistance3dSig_, reco::btau::flightDistance3dSig);

	return true;
}
