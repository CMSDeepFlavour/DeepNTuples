/*
 * ntuple_eventInfo.cc
 *
 *  Created on: 20 Feb 2018
 *      Author: dwalter
 */

#include "../interface/ntuple_eventInfo.h"



void ntuple_eventInfo::getInput(const edm::ParameterSet& iConfig){

    useLHEWeights_ = (iConfig.getParameter<bool>("useLHEWeights"));

    periods = (iConfig.getParameter<std::vector<std::string>>("periods"));
    lumis = (iConfig.getParameter<std::vector<double>>("lumis"));

    pupDataDir_ = (iConfig.getParameter<std::string>("pileupData"));
    pupMCDir_ = (iConfig.getParameter<std::string>("pileupMC"));

    sfTrigger_mu_Dir_ = (iConfig.getParameter<std::vector<std::string>>("sfTrigger_mu"));
    sfTrigger_mu_Name_ = (iConfig.getParameter<std::vector<std::string>>("sfTrigger_mu_Hist"));

    sfTrigger_emu_Dir_ = (iConfig.getParameter<std::vector<std::string>>("sfTrigger_emu"));
    sfTrigger_emu_Name_ = (iConfig.getParameter<std::vector<std::string>>("sfTrigger_emu_Hist"));
    sfMuonId_Dir_ = (iConfig.getParameter<std::vector<std::string>>("sfMuonId"));
    sfMuonId_Name_ = (iConfig.getParameter<std::vector<std::string>>("sfMuonId_Hist"));
    sfMuonIso_Dir_ = (iConfig.getParameter<std::vector<std::string>>("sfMuonIso"));
    sfMuonIso_Name_ = (iConfig.getParameter<std::vector<std::string>>("sfMuonIso_Hist"));
    sfElIdAndIso_Dir_ = (iConfig.getParameter<std::vector<std::string>>("sfElIdAndIso"));
    sfElIdAndIso_Name_ = (iConfig.getParameter<std::vector<std::string>>("sfElIdAndIso_Hist"));

    sfMuonTracking_Dir_ = (iConfig.getParameter<std::vector<std::string>>("sfMuonTracking"));
    sfMuonTracking_Name_ = (iConfig.getParameter<std::vector<std::string>>("sfMuonTracking_Hist"));

    triggers = (iConfig.getParameter<std::vector<std::string>>("triggers"));


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

    initializeScalefactor(sfTrigger_mu_Dir_, sfTrigger_mu_Name_, &sfTrigger_mu_Hist, periods);
    initializeScalefactor(sfTrigger_emu_Dir_, sfTrigger_emu_Name_, &sfTrigger_emu_Hist, periods);
    initializeScalefactor(sfMuonId_Dir_, sfMuonId_Name_, &sfMuonId_Hist, periods);
    initializeScalefactor(sfMuonIso_Dir_, sfMuonIso_Name_, &sfMuonIso_Hist, periods);
    initializeScalefactor(sfElIdAndIso_Dir_, sfElIdAndIso_Name_, &sfElIdAndIso_Hist, periods);
    initializeScalefactor(sfMuonTracking_Dir_, sfMuonTracking_Name_, &sfMuonTracking_Hist, periods);

}


void ntuple_eventInfo::initBranches(TTree* tree){

    addBranch(tree,"event_weight", &event_weight_,"event_weight_/f"    );
}


void ntuple_eventInfo::readEvent(const edm::Event& iEvent){
    TriggerInfo triggerInfo(iEvent,triggerToken_);

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

        ///// scalefactors

        double leadingMuon_pt = 0.;
        double leadingMuon_eta = 0.;
        for (size_t i = 0; i < muons->size(); ++i){
            const auto & muon = (*muons).at(i);
            if(muon.pt() > leadingMuon_pt){
                leadingMuon_pt = muon.pt();
                leadingMuon_eta = muon.eta();
            }
        }
        double leadingElectron_pt = 0.;
        double leadingElectron_sueta = 0.;
        double leadingElectron_eta = 0.;
        for (size_t i = 0; i < electrons->size(); ++i){
            const auto & electron = (*electrons).at(i);
            if(electron.pt() > leadingElectron_pt){
                leadingElectron_pt = electron.pt();
                leadingElectron_sueta = electron.superCluster()->eta();
                leadingElectron_eta = electron.eta();
            }
        }


        event_weight_ = 0.;
        for(unsigned int i=0; i< periods.size(); i++){

            double isf = 1.;
            bool itriggered = 1;

            isf *= getScalefactor(std::abs(leadingMuon_eta),        leadingMuon_pt,             sfTrigger_mu_Hist, i);
            isf *= getScalefactor(std::abs(leadingElectron_eta),    std::abs(leadingMuon_eta),  sfTrigger_emu_Hist, i);
            isf *= getScalefactor(std::abs(leadingMuon_eta),        leadingMuon_pt,             sfMuonId_Hist, i);
            isf *= getScalefactor(std::abs(leadingMuon_eta),        leadingMuon_pt,             sfMuonIso_Hist, i);
            isf *= getScalefactor(std::abs(leadingElectron_sueta),  leadingElectron_pt,         sfElIdAndIso_Hist, i);
            isf *= getScalefactor(std::abs(leadingMuon_eta),                                    sfMuonTracking_Hist, i);

            if(triggers.size() != 0){
                itriggered = triggerInfo.IsAnyTriggered(triggers[i]);
                //std::cout<<"istriggered "<<itriggered<<std::endl;
            }

            //std::cout<<getScalefactor(std::abs(leadingMuon_eta), leadingMuon_pt, sfTrigger_mu_Hist, i)<<std::endl;
            //std::cout<<getScalefactor(std::abs(leadingElectron_eta), std::abs(leadingMuon_eta), sfTrigger_emu_Hist, i)<<std::endl;
            //std::cout<<getScalefactor(std::abs(leadingMuon_eta), leadingMuon_pt, sfMuonId_Hist, i)<<std::endl;
            //std::cout<<getScalefactor(std::abs(leadingMuon_eta), leadingMuon_pt, sfMuonIso_Hist, i)<<std::endl;
            //std::cout<<getScalefactor(std::abs(leadingElectron_sueta), leadingElectron_pt, sfElIdAndIso_Hist, i)<<std::endl;
            //std::cout<<getScalefactor(std::abs(leadingMuon_eta), sfMuonTracking_Hist, i)<<std::endl;

            //std::cout<<"period "<<periods[i]<<" has sf = "<<isf<<std::endl;

            event_weight_ += lumis[i] * isf * itriggered;

            //std::cout<<"period "<<periods[i]<<" has part weight = "<<event_weight_<<std::endl;

        }
        event_weight_ *= (float)(lheWeight * pupWeight);

        //std::cout<<"total eventweight = "<<event_weight_<<std::endl;
    }
    else{
        event_weight_ = 1.;
    }

}

bool ntuple_eventInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
    return true;
}

double getScalefactor(double x, double y, std::vector<TH2F*> hists, unsigned int period){

    if(hists.size()==0) return 1.;
    else if(hists.size()==1){
        int binx = hists[0]->GetXaxis()->FindBin(x);
        int biny = hists[0]->GetYaxis()->FindBin(y);
        if(x > hists[0]->GetXaxis()->GetXmax())     //dont take overflow bins, but the last ones
            binx -= 1;
        if(y > hists[0]->GetYaxis()->GetXmax())
            biny -= 1;
        return(hists[0]->GetBinContent(binx, biny));
    }
    else if(hists.size() > period){
        int binx = hists[period]->GetXaxis()->FindBin(x);
        int biny = hists[period]->GetYaxis()->FindBin(y);
        if(x > hists[period]->GetXaxis()->GetXmax())     //dont take overflow bins, but the last ones
            binx -= 1;
        if(y > hists[period]->GetYaxis()->GetXmax())
            biny -= 1;
        return(hists[period]->GetBinContent(binx, biny));
    }
    else std::cout<<"ERROR: something went wrong in getScalefactor function"<<std::endl;

    return 0.;
}

double getScalefactor(double x, std::vector<TH1D*> hists, unsigned int period){

    if(hists.size()==0) return 1.;
    else if(hists.size()==1){
        int bin = hists[0]->GetXaxis()->FindBin(x);
        if(x > hists[0]->GetXaxis()->GetXmax())     //dont take overflow bins, but the last ones
            bin -= 1;
        return(hists[0]->GetBinContent(bin));
    }
    else if(hists.size() > period){
        int bin = hists[period]->GetXaxis()->FindBin(x);
        if(x > hists[period]->GetXaxis()->GetXmax())     //dont take overflow bins, but the last ones
            bin -= 1;
        return(hists[period]->GetBinContent(bin));
    }
    else std::cout<<"ERROR: something went wrong in getScalefactor function"<<std::endl;

    return 0.;
}

void initializeScalefactor(std::vector<std::string> dirs, std::vector<std::string> names, std::vector<TH2F*> *sfHists, std::vector<std::string> periods){
    if(dirs.size()==0)
        std::cout<<"no scalefactor "<<std::endl;
    else if(dirs.size()==periods.size()){
        for(unsigned int i = 0; i < periods.size(); i++){
            std::cout<<"initialize a scalefactor histogram for period "<<periods[i] <<std::endl;
            TFile *file = new TFile(dirs[i].c_str());
            sfHists->push_back((TH2F*)file->Get(names[i].c_str()));
        }
    }
    else if(dirs.size()==1){
        std::cout<<"initialize one scalefactor histogram for all periods " <<std::endl;
        TFile *file = new TFile(dirs[0].c_str());
        sfHists->push_back((TH2F*)file->Get(names[0].c_str()));
    }
    else{
        std::cout<<"ERROR: something went wrong by initialization of a scalefactor " <<std::endl;
    }
}

void initializeScalefactor(std::vector<std::string> dirs, std::vector<std::string> names, std::vector<TH1D*> *sfHists, std::vector<std::string> periods){
    if(dirs.size()==0)
        std::cout<<"no scalefactor "<<std::endl;
    else if(dirs.size()==periods.size()){
        for(unsigned int i = 0; i < periods.size(); i++){
            std::cout<<"initialize a scalefactor histogram for period "<<periods[i] <<std::endl;
            TGraphAsymmErrors* sfTGraph;
            TH1D* hist;

            TFile *file = new TFile(dirs[i].c_str());

            sfTGraph = (TGraphAsymmErrors*)file->Get(names[i].c_str());
            readHistoFromGraph(sfTGraph, &hist, "sfMuonTrackingHist");
            sfHists->push_back(hist);
        }
    }
    else if(dirs.size()==1){
        std::cout<<"initialize one scalefactor histogram for all periods " <<std::endl;
        TGraphAsymmErrors* sfTGraph;
        TH1D* hist;

        TFile *file = new TFile(dirs[0].c_str());

        sfTGraph = (TGraphAsymmErrors*)file->Get(names[0].c_str());
        readHistoFromGraph(sfTGraph, &hist, "sfMuonTrackingHist");
        sfHists->push_back(hist);
    }
    else{
        std::cout<<"ERROR: something went wrong by initialization of a scalefactor " <<std::endl;
    }
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
