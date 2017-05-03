/*
 * FatJetMatching.cc
 *
 *  Created on: Feb 1, 2017
 *      Author: hqu
 */

#include "DeepNTuples/JetAnalysis/interface/FatJetMatching.h"

#include <unordered_set>
#include "TString.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace deep_ntuples;

std::pair<FatJetMatching::FatJetFlavor, const reco::GenParticle*> FatJetMatching::flavor(const pat::Jet* jet,
    const reco::GenParticleCollection& genParticles) {

  processed_.clear();

  if (debug_) {
    std::cout << "\n=======\nJet (energy, pT, eta, phi) = "
        << jet->energy() << ", " << jet->pt() << ", " << jet->eta() << ", " << jet->phi()
        << std::endl << std::endl;
    printGenInfoHeader();
    for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
      printGenParticleInfo(&genParticles[ipart], ipart);
    }
  }

  for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
    const auto *gp = &genParticles[ipart];

    if (processed_.count(gp)) continue;
    processed_.insert(gp);

    auto pdgid = std::abs(gp->pdgId());
    if (pdgid == ParticleID::p_t){
      // top
      auto top = getFinal(gp);
      // find the W and test if it's hadronic
      const reco::GenParticle *w_from_top = nullptr, *b_from_top = nullptr;
      for (const auto &dau : top->daughterRefVector()){
        if (std::abs(dau->pdgId()) == ParticleID::p_Wplus){
          w_from_top = getFinal(&(*dau));
        }else if (std::abs(dau->pdgId()) <= ParticleID::p_b){
          // ! use <= p_b ! -- can also have charms etc.
          b_from_top = getFinal(&(*dau));
        }
      }
      if (!w_from_top || !b_from_top) throw std::logic_error("[FatJetMatching::flavor] Cannot find b or W from top decay: "+std::to_string(ipart));
      if (isHadronic(w_from_top)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "top: "; printGenParticleInfo(top, -1);
          cout << "b:   "; printGenParticleInfo(b_from_top, -1);
          cout << "W:   "; printGenParticleInfo(w_from_top, -1);
        }
        if (!requiresQuarksContained_) {
          double dr_top = reco::deltaR(jet->p4(), top->p4());
          if (debug_){
            using namespace std;
            cout << "deltaR(jet, top): " << dr_top << endl;
          }
          if (dr_top < jetR_) return std::make_pair(FatJetFlavor::Top, top);
        }
        else{
          double dr_wdaus = maxDeltaRToDaughterQuarks(jet, w_from_top);
          double dr_b     = reco::deltaR(jet->p4(), b_from_top->p4());
          if (debug_){
            using namespace std;
            cout << "deltaR(jet, b)     : " << dr_b << endl;
            cout << "deltaR(jet, w daus): " << dr_wdaus << endl;
          }
          if (dr_wdaus < jetR_ && dr_b < jetR_) return std::make_pair(FatJetFlavor::Top, top);
          if (dr_wdaus < jetR_ && dr_b > jetR_) return std::make_pair(FatJetFlavor::W, w_from_top);
        }
      }
    }else if (pdgid == ParticleID::p_Wplus){
      // W: not from top, or top not in jet cone
      auto w = getFinal(gp);
      if (isHadronic(w)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "W:   "; printGenParticleInfo(w, -1);
        }
        if (!requiresQuarksContained_){
          if (reco::deltaR(jet->p4(), w->p4()) < jetR_) return std::make_pair(FatJetFlavor::W, w);
        }
        else{
          double dr_wdaus = maxDeltaRToDaughterQuarks(jet, w);
          if (debug_){
            using namespace std;
            cout << "deltaR(jet, w daus): " << dr_wdaus << endl;
          }
          if (dr_wdaus < jetR_) return std::make_pair(FatJetFlavor::W, w);
        }
      }
    }else if (pdgid == ParticleID::p_Z0) {
      // Z
        auto z = getFinal(gp);
        if (isHadronic(z)) {
          if (debug_){
            using namespace std;
            cout << "jet: " << jet->polarP4() << endl;
            cout << "Z:   "; printGenParticleInfo(z, -1);
          }
          if (!requiresQuarksContained_){
            if (reco::deltaR(jet->p4(), z->p4()) < jetR_) return std::make_pair(FatJetFlavor::Z, z);
          }
          else{
            double dr_zdaus = maxDeltaRToDaughterQuarks(jet, z);
            if (debug_){
              using namespace std;
              cout << "deltaR(jet, Z daus): " << dr_zdaus << endl;
            }
            if (dr_zdaus < jetR_) return std::make_pair(FatJetFlavor::Z, z);
          }
        }
    }else if (pdgid == ParticleID::p_h0) {
      // Higgs
        auto h = getFinal(gp);
        if (isHadronic(h)) {
          if (debug_){
            using namespace std;
            cout << "jet: " << jet->polarP4() << endl;
            cout << "H:   "; printGenParticleInfo(h, -1);
          }
          if (!requiresQuarksContained_){
            if (reco::deltaR(jet->p4(), h->p4()) < jetR_) return std::make_pair(FatJetFlavor::H, h);
          }
          else{
            double dr_hdaus = maxDeltaRToDaughterQuarks(jet, h);
            if (debug_){
              using namespace std;
              cout << "deltaR(jet, h daus): " << dr_hdaus << endl;
            }
            if (dr_hdaus < jetR_) return std::make_pair(FatJetFlavor::H, h);
          }
        }
    }else {
    	// ?
    }
  }

  if (genParticles.size() != processed_.size())
    throw std::logic_error("[FatJetMatching::flavor] Not all genParticles are processed!");

  return std::make_pair(FatJetFlavor::Default, nullptr);

}

void FatJetMatching::printGenInfoHeader() const {
  using namespace std;
  cout    << right << setw(6) << "#" << " " << setw(10) << "pdgId"
      << "  " << "Chg" << "  " << setw(10) << "Mass" << "  " << setw(48) << " Momentum"
      << left << "  " << setw(10) << "Mothers" << " " << setw(30) << "Daughters" << endl;
}

void FatJetMatching::printGenParticleInfo(const reco::GenParticle* genParticle, const int idx) const {
  using namespace std;
  cout  << right << setw(3) << genParticle->status();
  cout  << right << setw(3) << idx << " " << setw(10) << genParticle->pdgId() << "  ";
  cout  << right << "  " << setw(3) << genParticle->charge() << "  " << TString::Format("%10.3g", genParticle->mass() < 1e-5 ? 0 : genParticle->mass());
  cout  << left << setw(50) << TString::Format("  (E=%6.4g pT=%6.4g eta=%7.3g phi=%7.3g)", genParticle->energy(), genParticle->pt(), genParticle->eta(), genParticle->phi());

  TString                     mothers;
  for (unsigned int iMom = 0; iMom < genParticle->numberOfMothers(); ++iMom) {
    if (mothers.Length())     mothers        += ",";
    mothers   += genParticle->motherRef(iMom).key();
  }
  cout << "  " << setw(10) << mothers;
  TString                     daughters;
  for (unsigned int iDau = 0; iDau < genParticle->numberOfDaughters(); ++iDau) {
    if (daughters.Length())   daughters      += ",";
    daughters += genParticle->daughterRef(iDau).key();
  }
  cout << " " << setw(30) << daughters << endl;
}

const reco::GenParticle* FatJetMatching::getFinal(const reco::GenParticle* particle) {
  // will mark intermediate particles as processed
  if (!particle) return nullptr;
  processed_.insert(particle);
  const reco::GenParticle *final = particle;

  while (final->numberOfDaughters()) {
    const reco::GenParticle *chain = nullptr;
    for (unsigned idau = 0; idau < final->numberOfDaughters(); ++idau){
      if (final->daughter(idau)->pdgId() == particle->pdgId()) {
        chain = dynamic_cast<const reco::GenParticle*>(final->daughter(idau));
        processed_.insert(chain);
        break;
      }
    }
    if (!chain) break;
    final = chain;
  }
  return final;
}

bool FatJetMatching::isHadronic(const reco::GenParticle* particle) const {
  // particle needs to be the final version before decay
  if (!particle) throw std::invalid_argument("[FatJetMatching::isHadronic()] Null particle!");
  for(const auto &dau : particle->daughterRefVector()){
    auto pdgid = std::abs(dau->pdgId());
    if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b) return true;
  }
  return false;
}

double FatJetMatching::maxDeltaRToDaughterQuarks(const pat::Jet* jet, const reco::GenParticle* mother) const {
  // mother particle needs to be the final version before decay
  double maxDeltaR = -1;
  for (const auto &q : mother->daughterRefVector()){
    if (std::abs(q->pdgId()) > ParticleID::p_b) continue;
    double deltaR = reco::deltaR(q->p4(), jet->p4());
    if (deltaR > maxDeltaR) maxDeltaR = deltaR;
  }
  return maxDeltaR > 0 ? maxDeltaR : 1e9;
}
