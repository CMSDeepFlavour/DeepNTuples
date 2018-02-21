/*
 * ntuple_eventInfo.cc
 *
 *  Created on: 20 Feb 2018
 *      Author: dwalter
 */

#include "../interface/ntuple_eventInfo.h"


#include <vector>



void ntuple_eventInfo::getInput(const edm::ParameterSet& iConfig){

    useLHEWeights_ = (iConfig.getParameter<bool>("useLHEWeights"));

}
void ntuple_eventInfo::initBranches(TTree* tree){

    addBranch(tree,"event_weight", &event_weight_);
}


void ntuple_eventInfo::readEvent(const edm::Event& iEvent){

    if(!iEvent.isRealData()){

        iEvent.getByToken(lheToken_, lheInfo);
    }

    if(!iEvent.isRealData()){
        for (auto const& v : *pupInfo()) {
            int bx = v.getBunchCrossing();
            if (bx == 0) {
                ntrueInt_ = v.getTrueNumInteractions();
            }
        }
        double lheWeight = 1.;

        if(useLHEWeights_){
            lheWeight = lheInfo->weights()[0].wgt/std::abs(lheInfo->weights()[0].wgt);
        }

        double pupWeight = 0;
        if(ntrueInt_ < pupWeights.size()){
            pupWeight = pupWeights.at(ntrueInt_);
        }

        event_weight_ = lheWeight * pupWeight;
    }
    else{
        event_weight_ = 1.;
    }

}

//use either of these functions

bool ntuple_eventInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
    return true;
}
