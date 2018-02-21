/*
 * ntuple_eventInfo.cc
 *
 *  Created on: 20 Feb 2018
 *      Author: dwalter
 */

#include "../interface/ntuple_eventInfo.h"



void ntuple_eventInfo::getInput(const edm::ParameterSet& iConfig){

    useLHEWeights_ = (iConfig.getParameter<bool>("useLHEWeights"));
    pupDataDir_ = (iConfig.getParameter<std::string>("pileupDataDir"));
    pupMCDir_ = (iConfig.getParameter<std::string>("pileupMCDir"));

    if(pupDataDir_=="" || pupMCDir_==""){
        std::cout<<"no pileup histograms, proceed without pileup reweighting. \n";
    }
    else{
        pupMCFile = new TFile(pupMCDir_.c_str());
        pupDataFile = new TFile(pupDataDir_.c_str());

        pupMCHist = (TH1F*)pupMCFile->Get("pileup");
        pupDataHist = (TH1F*)pupDataFile->Get("pileup");

        pupMCHist->Scale(1./pupMCHist->Integral());
        pupDataHist->Scale(1./pupDataHist->Integral());

        pupWeights.push_back(1.0);
        for(int bin = 1; bin < pupMCHist->GetNbinsX() + 1; bin++){
            pupWeights.push_back(pupDataHist->GetBinContent(bin)/pupMCHist->GetBinContent(bin));
        }
    }



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

        double pupWeight = 1.;
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
