#include "../interface/helpers.h"

namespace deep_ntuples {
	JetFlavor jet_flavour(const pat::Jet& jet) {
		int hflav = abs(jet.hadronFlavour());
		int pflav = abs(jet.partonFlavour());
		size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
		size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

		if(hflav == 5) { //B jet
			if(nbs > 1) return JetFlavor::BB;
			else if(nbs == 1) return JetFlavor::B;
			else return JetFlavor::UNDEFINED;			
		}
		else if(hflav == 4) { //C jet
			if(ncs > 1) return JetFlavor::CC;
			else if(ncs == 1) return JetFlavor::C;
			else return JetFlavor::UNDEFINED;			
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
