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
    sfMuonIdName_ = (iConfig.getParameter<std::string>("sfMuonIdHist"));
    sfMuonIsoDir_ = (iConfig.getParameter<std::string>("sfMuonIso"));
    sfMuonIsoName_ = (iConfig.getParameter<std::string>("sfMuonIsoHist"));
    sfMuonTrackingDir_ = (iConfig.getParameter<std::string>("sfMuonTracking"));
    sfMuonTrackingName_ = (iConfig.getParameter<std::string>("sfMuonTrackingHist"));
    sfElIdAndIsoDir_ = (iConfig.getParameter<std::string>("sfElIdAndIso"));
    sfElIdAndIsoName_ = (iConfig.getParameter<std::string>("sfElIdAndIsoHist"));

    if(pupDataDir_=="" || pupMCDir_=="")
        std::cout<<"no pileup histograms, proceed without pileup reweighting. \n";
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

    if(sfMuonIdDir_ == "")
        std::cout<<"no muon id scalefactor histogram, proceed without muon id scalefactor"<<std::endl;
    else{
        TFile *sfMuonIdFile = new TFile(sfMuonIdDir_.c_str());

        sfMuonIdHist = (TH2F*)sfMuonIdFile->Get(sfMuonIdName_.c_str());

        sfMuonIdHist_xaxis = sfMuonIdHist->GetXaxis();
        sfMuonIdHist_yaxis = sfMuonIdHist->GetYaxis();
    }

    if(sfMuonIsoDir_ == "")
        std::cout<<"no muon iso scalefactor histogram, proceed without muon id scalefactor"<<std::endl;
    else{
        TFile *sfMuonIsoFile = new TFile(sfMuonIsoDir_.c_str());

        sfMuonIsoHist = (TH2F*)sfMuonIsoFile->Get(sfMuonIsoName_.c_str());

        sfMuonIsoHist_xaxis = sfMuonIsoHist->GetXaxis();
        sfMuonIsoHist_yaxis = sfMuonIsoHist->GetYaxis();
    }

    if(sfElIdAndIsoDir_ == ""){
        std::cout<<"no electron id and iso scalefactor histogram, proceed without electron id and iso scalefactor"<<std::endl;
    }
    else{
        TFile *sfElIdAndIsoFile = new TFile(sfElIdAndIsoDir_.c_str());

        sfElIdAndIsoHist = (TH2F*)sfElIdAndIsoFile->Get(sfElIdAndIsoName_.c_str());

        sfElIdAndIsoHist_xaxis = sfElIdAndIsoHist->GetXaxis();
        sfElIdAndIsoHist_yaxis = sfElIdAndIsoHist->GetYaxis();
    }

    if(sfMuonTrackingDir_ == ""){
        std::cout<<"no electron id and iso scalefactor histogram, proceed without electron id and iso scalefactor"<<std::endl;
    }
    else{
        TFile *sfMuonTrackingFile = new TFile(sfMuonTrackingDir_.c_str());

        sfMuonTrackingTGraph = (TGraphAsymmErrors*)sfMuonTrackingFile->Get(sfMuonTrackingName_.c_str());

        readHistoFromGraph(sfMuonTrackingTGraph, &sfMuonTrackingHist, "sfMuonTrackingHist");

        sfMuonTrackingHist_axis = sfMuonTrackingHist->GetXaxis();
    }
}


void ntuple_eventInfo::initBranches(TTree* tree){

    addBranch(tree,"event_weight", &event_weight_);
}


void ntuple_eventInfo::readEvent(const edm::Event& iEvent){

    if(!iEvent.isRealData()){

        iEvent.getByToken(lheToken_, lheInfo);
        iEvent.getByToken(muonToken_, muons);
        iEvent.getByToken(electronToken_, electrons);

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
        // Muon ISO
        double muonIsoSf = 1.;
        if(sfMuonIsoDir_ != ""){
            int binx = sfMuonIsoHist_xaxis->FindBin(std::abs(leadingMuon_eta));
            int biny = sfMuonIsoHist_yaxis->FindBin(leadingMuon_pt);
            if(leadingMuon_pt > 120.)   //dont take overflow bin, but the last one
                biny -= 1;
            muonIsoSf = sfMuonIsoHist->GetBinContent(binx, biny);
        }
        //Muon tracking
        double muonTrackingSf = 1.;
        if(sfMuonTrackingDir_ != ""){
            int binx = sfMuonTrackingHist_axis->FindBin(std::abs(leadingMuon_eta));
            muonTrackingSf = sfMuonTrackingHist->GetBinContent(binx);
        }

        // Electron scalefactors
        double leadingElectron_pt = 0.;
        double leadingElectron_sueta = 0.;
        for (size_t i = 0; i < electrons->size(); ++i){
            const auto & electron = (*electrons).at(i);
            if(electron.pt() > leadingElectron_pt){
                leadingElectron_pt = electron.pt();
                leadingElectron_sueta = electron.superCluster()->eta();
            }
        }
        //Electron ID and ISO
        double elIdAndIsoSf = 1.;
        if(sfElIdAndIsoDir_ != ""){
            int binx = sfElIdAndIsoHist_xaxis->FindBin(std::abs(leadingElectron_sueta));
            int biny = sfElIdAndIsoHist_yaxis->FindBin(leadingElectron_pt);
            if(leadingElectron_pt > 120.)
                biny -= 1;
            elIdAndIsoSf = sfElIdAndIsoHist->GetBinContent(binx, biny);
        }


        event_weight_ = lheWeight * pupWeight * muonIdSf * muonIsoSf * muonTrackingSf * elIdAndIsoSf;
    }
    else{
        event_weight_ = 1.;
    }

}

//use either of these functions

bool ntuple_eventInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
    return true;
}

void readHistoFromGraph(TGraphAsymmErrors* graph, TH1D** h, TString name)
{
    const int npoints = graph->GetN();
    const double* x_centers = graph->GetX();
    const double* y_centers = graph->GetY();
    double x_lows[npoints];
    double x_highs[npoints];
    double y_lows[npoints];
    double y_highs[npoints];
    for(int i=0; i<npoints; ++i)
    {
        x_lows[i] = graph->GetErrorXlow(i);
        x_highs[i] = graph->GetErrorXhigh(i);
        y_lows[i] = graph->GetErrorYlow(i);
        y_highs[i] = graph->GetErrorYhigh(i);
    }

    double x_edges[npoints+1];
    for(int i=0; i<npoints; ++i)
    {
        x_edges[i] = x_centers[i] - x_lows[i];
    }
    x_edges[npoints] = x_centers[npoints-1] + x_highs[npoints-1];


    *h = new TH1D(name, name, npoints, x_edges);
    (*h)->SetDirectory(0); // without this the histo will get deleted when a currently open TFile is closed

    for(int i=0; i<npoints; ++i)
    {
//     cout << i << ":" << y_centers[i] << "+-" << max(y_lows[i], y_highs[i]) << endl;
        (*h)->SetBinContent(i+1, y_centers[i]);
        (*h)->SetBinError(i+1, std::max(y_lows[i], y_highs[i]));
//     cout << (*h)->GetBinError(i) << endl;
    }
}
