#include "../interface/helpers.h"

namespace deep_ntuples {
JetFlavor jet_flavour(const pat::Jet& jet, std::vector<reco::GenParticle> neutrinosLepB, std::vector<reco::GenParticle> neutrinosLepB_C, bool usePhysForLightAndUndefined) {
    int hflav = abs(jet.hadronFlavour());
    int pflav = abs(jet.partonFlavour());
    int physflav = 0;
    if(jet.genParton()) physflav=abs(jet.genParton()->pdgId());
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
        else {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return JetFlavor::G;
                else if(physflav == 3) return JetFlavor::S;
                else if(physflav == 2 || physflav ==1) return JetFlavor::UD;
                else return JetFlavor::UNDEFINED;
            }
            else return JetFlavor::UNDEFINED;
        }
    }
    else if(hflav == 4) { //C jet
        return JetFlavor::C;
    }
    else { //not a heavy jet
        if(std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs) {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return JetFlavor::G;
                else if(physflav == 3) return JetFlavor::S;
                else if(physflav == 2 || physflav ==1) return JetFlavor::UD;
                else return JetFlavor::UNDEFINED;
            }
            else return JetFlavor::UNDEFINED;
        }
        else if(usePhysForLightAndUndefined){
            if(physflav == 21) return JetFlavor::G;
            else if(physflav == 3) return JetFlavor::S;
            else if(physflav == 2 || physflav ==1) return JetFlavor::UD;
            else return JetFlavor::UNDEFINED;
        }
        else {
            if(pflav == 21) return JetFlavor::G;
            else if(pflav == 3) return JetFlavor::S;
            else if(pflav == 2 || pflav ==1) return JetFlavor::UD;
            else return JetFlavor::UNDEFINED;
        }
    }
}
}
