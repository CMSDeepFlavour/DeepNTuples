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
void initializeScalefactor(std::vector<std::string> dirs, std::vector<std::string> names, std::vector<TH2F*>* sfHists, std::vector<std::string> periods);
void initializeScalefactor(std::vector<std::string> dirs, std::vector<std::string> names, std::vector<TH1D*>* sfHists, std::vector<std::string> periods);
double getScalefactor(double x, double y, std::vector<TH2F*> hists, unsigned int period);
double getScalefactor(double x, std::vector<TH1D*> hists, unsigned int period);


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

    std::vector<double> lumis;
    std::vector<std::string> periods;

    std::string pupDataDir_;
    std::string pupMCDir_;

    //std::string sfTrigger_mu_BCDEF_Dir_;
    //std::string sfTrigger_mu_BCDEF_Name_;
    //std::string sfTrigger_mu_GH_Dir_;
    //std::string sfTrigger_mu_GH_Name_;

    std::vector<std::string> sfTrigger_mu_Dir_;
    std::vector<std::string> sfTrigger_mu_Name_;
    std::vector<std::string> sfTrigger_emu_Dir_;
    std::vector<std::string> sfTrigger_emu_Name_;
    std::vector<std::string> sfMuonId_Dir_;
    std::vector<std::string> sfMuonId_Name_;
    std::vector<std::string> sfMuonIso_Dir_;
    std::vector<std::string> sfMuonIso_Name_;
    std::vector<std::string> sfElIdAndIso_Dir_;
    std::vector<std::string> sfElIdAndIso_Name_;
    std::vector<std::string> sfMuonTracking_Dir_;
    std::vector<std::string> sfMuonTracking_Name_;

    std::vector<TH2F*> sfTrigger_mu_Hist;
    std::vector<TH2F*> sfTrigger_emu_Hist;
    std::vector<TH2F*> sfMuonId_Hist;
    std::vector<TH2F*> sfMuonIso_Hist;
    std::vector<TH2F*> sfElIdAndIso_Hist;

    std::vector<TH1D*> sfMuonTracking_Hist;


    // global variables

    float ntrueInt_;
    std::vector<double> pupWeights;

    /////////branches
    float event_weight_;



};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_ */
