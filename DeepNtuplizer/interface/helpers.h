/*
 *      Author: mverzett
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace deep_ntuples {
	enum JetFlavor {UNDEFINED, G, L, C, CC, B, BB};
	JetFlavor jet_flavour(const pat::Jet& jet);
}

#endif //DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_HELPERS_H_
