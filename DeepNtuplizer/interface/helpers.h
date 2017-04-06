/*
 *      Author: mverzett
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

namespace deep_ntuples {
    enum JetFlavor {UNDEFINED, G, UD, S, C, B, BB, LeptonicB, LeptonicB_C};
    JetFlavor jet_flavour(const pat::Jet& jet, std::vector<reco::GenParticle> neutrinosLepB, std::vector<reco::GenParticle> neutrinosLepB_C, bool usePhysForLightAndUndefined=false);
    std::vector<std::size_t> jet_muonsIds(const pat::Jet& jet, const std::vector<pat::Muon>& event_muons); 
}

#endif //DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
