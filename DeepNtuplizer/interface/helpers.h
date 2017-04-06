/*
 *      Author: mverzett
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

namespace deep_ntuples {
    std::vector<std::size_t> jet_muonsIds(const pat::Jet& jet, const std::vector<pat::Muon>& event_muons); 
    enum JetFlavor {UNDEFINED, G, L, C, B, BB, LeptonicB, LeptonicB_C};
    JetFlavor jet_flavour(const pat::Jet& jet, std::vector<reco::GenParticle> neutrinosLepB, std::vector<reco::GenParticle> neutrinosLepB_C);
}

#endif //DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
