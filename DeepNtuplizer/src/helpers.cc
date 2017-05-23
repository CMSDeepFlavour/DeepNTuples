#include "../interface/helpers.h"
#include <vector>

namespace deep_ntuples {

    std::vector<std::size_t> jet_muonsIds(const pat::Jet& jet, const std::vector<pat::Muon>& event_muons) {
        std::vector <std::size_t> muonsIds;
        for (std::size_t i = 0; i < event_muons.size(); i++) {
            const auto & muon = event_muons.at(i);
            if(reco::deltaR(muon.eta(),muon.phi(),jet.eta(),jet.phi()) < 0.4) muonsIds.emplace_back(i);
        }
        return muonsIds;
    }

    std::vector<std::size_t> jet_electronsIds(const pat::Jet& jet, const std::vector<pat::Electron>& event_electrons) {
        std::vector <std::size_t> electronsIds;
        for (std::size_t i = 0; i < event_electrons.size(); i++) {
            const auto & electron = event_electrons.at(i);
            if(reco::deltaR(electron.eta(),electron.phi(),jet.eta(),jet.phi()) < 0.4) electronsIds.emplace_back(i);
        }
        return electronsIds;
    }

    JetFlavor jet_flavour(const pat::Jet& jet, std::vector<reco::GenParticle> gToBB, std::vector<reco::GenParticle> gToCC, std::vector<reco::GenParticle> neutrinosLepB, std::vector<reco::GenParticle> neutrinosLepB_C, bool usePhysForLightAndUndefined) { 
        int hflav = abs(jet.hadronFlavour());
        int pflav = abs(jet.partonFlavour());
        int physflav = 0;
        if(jet.genParton()) physflav=abs(jet.genParton()->pdgId());
        std::size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
        std::size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

        unsigned int nbFromGSP(0);
        for (reco::GenParticle p : gToBB) {
          double dr(reco::deltaR(jet, p));
          if (dr < 0.4) ++nbFromGSP;
        }

        unsigned int ncFromGSP(0);
        for (reco::GenParticle p : gToCC) {
          double dr(reco::deltaR(jet, p));
          if (dr < 0.4) ++ncFromGSP;
        }

        std::cout << " jet pt = " << jet.pt() << " hfl = " << hflav << " pfl = " << pflav << " genpart = " << physflav 
          << " nbFromGSP = " << nbFromGSP << " ncFromGSP = " << ncFromGSP
          << " nBhadrons " << nbs << " nCHadrons " << ncs << std::endl;

        if(hflav == 5) { //B jet
            if(nbs > 1) {
              if (nbFromGSP > 0) return JetFlavor::GBB; 
              else return JetFlavor::BB;
            }
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
            if (ncs > 1) {
              if (ncFromGSP > 0) return JetFlavor::GCC;
              else return JetFlavor::CC;
            }
            else return JetFlavor::C;
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
