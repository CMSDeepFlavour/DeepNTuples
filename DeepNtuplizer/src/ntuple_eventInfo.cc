/*
 * ntuple_eventInfo.cc
 *
 *  Created on: 20 Feb 2018
 *      Author: dwalter
 */

#include "../interface/ntuple_eventInfo.h"



void ntuple_eventInfo::getInput(const edm::ParameterSet& iConfig){

    useLHEWeights_ = (iConfig.getParameter<bool>("useLHEWeights"));

    pupDataDir_ = (iConfig.getParameter<std::string>("pileupData"));
    pupMCDir_ = (iConfig.getParameter<std::string>("pileupMC"));
    sfMuonIdDir_ = (iConfig.getParameter<std::string>("sfMuonId"));
    sfMuonIdName_ = (iConfig.getParameter<std::string>("sfMuonIdName"));


    if(pupDataDir_=="" || pupMCDir_==""){
        std::cout<<"no pileup histograms, proceed without pileup reweighting. \n";
    }
    else{
        TFile *pupMCFile = new TFile(pupMCDir_.c_str());
        TFile *pupDataFile = new TFile(pupDataDir_.c_str());

        TH1F *pupMCHist = (TH1F*)pupMCFile->Get("pileup");
        TH1F *pupDataHist = (TH1F*)pupDataFile->Get("pileup");

        pupMCHist->Scale(1./pupMCHist->Integral());
        pupDataHist->Scale(1./pupDataHist->Integral());

        pupWeights.push_back(1.0);
        for(int bin = 1; bin < pupMCHist->GetNbinsX() + 1; bin++){
            pupWeights.push_back(pupDataHist->GetBinContent(bin)/pupMCHist->GetBinContent(bin));
        }
    }

    if(sfMuonIdDir_ == ""){
        std::cout<<"no muon id scalefactor histogram, proceed without muon id scalefactor"<<std::endl;
    }
    else{
        TFile *sfMuonIdFile = new TFile(sfMuonIdDir_.c_str());

        sfMuonIdHist = (TH2F*)sfMuonIdFile->Get(sfMuonIdName_.c_str());

        sfMuonIdHist_xaxis = sfMuonIdHist->GetXaxis();
        sfMuonIdHist_yaxis = sfMuonIdHist->GetYaxis();
    }
}


void ntuple_eventInfo::initBranches(TTree* tree){

    addBranch(tree,"event_weight", &event_weight_);
}


void ntuple_eventInfo::readEvent(const edm::Event& iEvent){

    if(!iEvent.isRealData()){

        iEvent.getByToken(lheToken_, lheInfo);
        iEvent.getByToken(muonToken_, muons);

        ///// pileup weights
        for (auto const& v : *pupInfo()) {
            int bx = v.getBunchCrossing();
            if (bx == 0) {
                ntrueInt_ = v.getTrueNumInteractions();
            }
        }
        double pupWeight = 1.;
        if(ntrueInt_ < pupWeights.size()){
            pupWeight = pupWeights.at(ntrueInt_);
        }

        ///// lhe weights
        double lheWeight = 1.;
        if(useLHEWeights_){
            lheWeight = lheInfo->weights()[0].wgt/std::abs(lheInfo->weights()[0].wgt);
        }

        /////scalefactors
        // Muon scalefactors
        double leadingMuon_pt = 0.;
        double leadingMuon_eta = 0.;
        for (size_t i = 0; i < muons->size(); ++i){
            const auto & muon = (*muons).at(i);
            if(muon.pt() > leadingMuon_pt){
                leadingMuon_pt = muon.pt();
                leadingMuon_eta = muon.eta();
            }
        }
        // Muon ID
        double muonIdSf = 1.;
        if(sfMuonIdDir_ != ""){
            int binx = sfMuonIdHist_xaxis->FindBin(std::abs(leadingMuon_eta));
            int biny = sfMuonIdHist_yaxis->FindBin(leadingMuon_pt);
            if(leadingMuon_pt > 120.)   //dont take overflow bin, but the last one
                biny -= 1;
            muonIdSf = sfMuonIdHist->GetBinContent(binx, biny);
        }

        std::cout<<"leading muon has pt "<<leadingMuon_pt<<" and eta "<<leadingMuon_eta<<" muon id scalefactor = "<<muonIdSf<<std::endl;



        event_weight_ = lheWeight * pupWeight * muonIdSf;
    }
    else{
        event_weight_ = 1.;
    }

}

//use either of these functions

bool ntuple_eventInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
    return true;
}
