#include "../interface/helpers.h"
#include <vector>

namespace deep_ntuples {
    std::vector<std::size_t> jet_muonsIds(const pat::Jet& jet, const std::vector<pat::Muon>& event_muons) {
        std::vector <std::size_t> muonIds;
        for (std::size_t i = 0; i < event_muons.size(); i++) {
            const auto & muon = event_muons.at(i);
            if(reco::deltaR(muon.eta(),muon.phi(),jet.eta(),jet.phi()) < 0.4) muonIds.emplace_back(i);
        }
        return muonIds;
    }
    JetFlavor jet_flavour(const pat::Jet& jet, std::vector<reco::GenParticle> neutrinosLepB, std::vector<reco::GenParticle> neutrinosLepB_C) {
        int hflav = abs(jet.hadronFlavour());
        int pflav = abs(jet.partonFlavour());
        size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
        size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

        if(hflav == 5) { //B jet
            if(nbs > 1) return JetFlavor::BB;
            else if(nbs == 1) {
                for (std::vector<reco::GenParticle>::iterator it = neutrinosLepB.begin(); it != neutrinosLepB.end(); ++it){
                    if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4) {
                        return JetFlavor::LeptonicB; 
                    }
                }
                for (std::vector<reco::GenParticle>::iterator it = neutrinosLepB_C.begin(); it != neutrinosLepB_C.end(); ++it){
                    if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4) {
                        return JetFlavor::LeptonicB_C; 
                    }
                }
                return JetFlavor::B; 
            }
            else return JetFlavor::UNDEFINED;           
        }
        else if(hflav == 4) { //C jet
            return JetFlavor::C;
        }
        else { //not a heavy jet
            if(std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs) {
                return JetFlavor::UNDEFINED;
            }
            else {
                if(pflav == 21) return JetFlavor::G;
                else return JetFlavor::L;
            }
        }
    }
}
