/*
 * ntuple_eventInfo.h
 *
 *  Created on: 20 Feb 2018
 *      Author: dwalter
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_

#include "ntuple_content.h"

#include <vector>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "TGraphAsymmErrors.h"
#include <TObjString.h>



#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

void readHistoFromGraph(TGraphAsymmErrors* graph, TH1D** h, TString name);

/*
 * For MC weights such as pileup, lhe, ... later: lepton scalefactors
 */
class ntuple_eventInfo: public ntuple_content{
public:
    ntuple_eventInfo():ntuple_content(),
    isData_(false)
{}

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);

    void setLHEToken(edm::EDGetTokenT<LHEEventProduct> lheToken) {
        lheToken_ = lheToken;
    }
    void setMuonsToken(edm::EDGetTokenT<edm::View<pat::Muon> > muonToken){
        muonToken_ = muonToken;
    }
    void setElectronsToken(edm::EDGetTokenT<edm::View<pat::Electron> > electronToken){
        electronToken_ = electronToken;
    }


private:

    edm::EDGetTokenT<LHEEventProduct> lheToken_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muonToken_;
    edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_;

    edm::Handle<LHEEventProduct> lheInfo;
    edm::Handle<edm::View<pat::Muon> > muons;
    edm::Handle<edm::View<pat::Electron> > electrons;


    bool isData_;
    bool useLHEWeights_;

    std::string pupDataDir_;
    std::string pupMCDir_;
    std::string sfMuonTriggerDir_;
    std::string sfMuonTriggerName_;
    std::string sfMuonIdDir_;
    std::string sfMuonIdName_;
    std::string sfMuonIsoDir_;
    std::string sfMuonIsoName_;
    std::string sfMuonTrackingDir_;
    std::string sfMuonTrackingName_;
    std::string sfElIdAndIsoDir_;
    std::string sfElIdAndIsoName_;

    TH2F *sfMuonTriggerHist;
    TH2F *sfMuonIdHist;
    TH2F *sfMuonIsoHist;
    TH2F *sfElIdAndIsoHist;
    TGraphAsymmErrors* sfMuonTrackingTGraph;
    TH1D *sfMuonTrackingHist;

    TAxis *sfMuonTriggerHist_xaxis;
    TAxis *sfMuonTriggerHist_yaxis;
    TAxis *sfMuonIdHist_xaxis;
    TAxis *sfMuonIdHist_yaxis;
    TAxis *sfMuonIsoHist_xaxis;
    TAxis *sfMuonIsoHist_yaxis;
    TAxis *sfElIdAndIsoHist_xaxis;
    TAxis *sfElIdAndIsoHist_yaxis;
    TAxis *sfMuonTrackingHist_axis;

    // global variables

    float ntrueInt_;
    std::vector<double> pupWeights;

    /////////branches
    float event_weight_;


};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_ */
