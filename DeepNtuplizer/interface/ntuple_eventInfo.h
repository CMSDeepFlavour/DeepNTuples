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
#include <TObjString.h>



#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

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
    void setMuonsToken(edm::EDGetTokenT<edm::View<pat::Muon>> muonToken){
        muonToken_ = muonToken;
    }


private:

    edm::EDGetTokenT<LHEEventProduct> lheToken_;
    edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;

    edm::Handle<LHEEventProduct> lheInfo;
    edm::Handle<edm::View<pat::Muon>> muons;

    bool isData_;
    bool useLHEWeights_;

    std::string pupDataDir_;
    std::string pupMCDir_;
    std::string sfMuonIdDir_;
    std::string sfMuonIdName_;

    TH2F *sfMuonIdHist;

    TAxis *sfMuonIdHist_xaxis;
    TAxis *sfMuonIdHist_yaxis;

    // global variables

    float ntrueInt_;

    std::vector<double> pupWeights;

    /////////branches
    float event_weight_;


};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_ */
