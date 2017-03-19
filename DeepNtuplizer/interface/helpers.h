/*
 *      Author: mverzett,George Karathanasis
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"


namespace deep_ntuples {
  enum JetFlavor {UNDEFINED, G, L, C,Blep,Btau, B, BB};
  JetFlavor jet_flavour(const pat::Jet& jet,  const reco::GenParticleCollection* PgenC );
  float parton_pt(const pat::Jet& jet,  const reco::GenParticleCollection* PgenC );
}

#endif //DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
