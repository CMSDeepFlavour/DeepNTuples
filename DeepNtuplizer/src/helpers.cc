/* Mauro Verzeti, George Karathanasis */

#include "../interface/helpers.h"


namespace deep_ntuples {
  JetFlavor jet_flavour(const pat::Jet& jet, const reco::GenParticleCollection* PgenC) {
    int hflav = abs(jet.hadronFlavour());
    int pflav = abs(jet.partonFlavour());
    size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
    size_t ncs = jet.jetFlavourInfo().getcHadrons().size();
    //std::cout<<PgenC->size()<<std::endl;
       
		
    //const reco::GenParticle & genC= prun_gen_parts()->at(0);
    std::vector<int> tau;
   std::vector < reco::GenParticle> neutrinos;
   //for (unsigned int i=0; i<PgenC->size(); i++ )
   for (const reco::Candidate &genC : *PgenC )
    {
      const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
      //std::cout<<gen.pt()<<std::endl; 
     //for (std::vector<reco::GenParticle>::iterator it = PgenC->at(i).begin(); it !=  PgenC->at(i).end(); ++it)
      //for (unsigned int j=0; j<PgenC->at(i)->size(); j++ )
      //std::cout<<PgenC->at(i)<<std::endl;
     
	if(abs(gen.pdgId())==12||abs(gen.pdgId())==14||abs(gen.pdgId())==16){
	 if((gen.fromHardProcessDecayed()+ gen.fromHardProcessFinalState()+ gen.isHardProcess()+gen.fromHardProcessBeforeFSR()) != 0 ) continue;
         const reco::GenParticle* mother =  static_cast< const reco::GenParticle*> (gen.mother());
         if (mother==NULL) continue;
         
	 if(abs(mother->pdgId())==15&& (abs(mother->mother()->mother()->pdgId())==24 || abs(mother->mother()->pdgId())==24)) /*continue;*/ {tau.push_back(1);  neutrinos.push_back(gen);}
	 if((abs(mother->pdgId())>500&&abs(mother->pdgId())<600)||(abs(mother->pdgId())>5000&&abs(mother->pdgId())<6000)) {neutrinos.push_back(gen); tau.push_back(0);}
	}
    }
   unsigned int tau_counter=0;  
   for (std::vector<reco::GenParticle>::iterator it = neutrinos.begin(); it != neutrinos.end(); ++it){
     if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4 && hflav==5)        // return JetFlavor::Blep; 
      {if(tau[tau_counter]==1) return JetFlavor::Btau;
     else  return JetFlavor::Blep;  }
     tau_counter+=1;


      }

	
		if(hflav == 5) { //B jet
			if(nbs > 1) return JetFlavor::BB;
			else if(nbs == 1) return JetFlavor::B;
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

  float  parton_pt(const pat::Jet& jet, const reco::GenParticleCollection* PgenC) {
   std::vector < reco::GenParticle> bhm;
   //for (unsigned int i=0; i<PgenC->size(); i++ )
   int hflav = abs(jet.hadronFlavour());
   for (const reco::Candidate &genC : *PgenC )
    {
      const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);          
	 if((abs(gen.pdgId())>500&&abs(gen.pdgId())<600)||(abs(gen.pdgId())>5000&&abs(gen.pdgId())<6000)) bhm.push_back(gen); 
    }
    
     
   for (std::vector<reco::GenParticle>::iterator it = bhm.begin(); it != bhm.end(); ++it){
     if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4 && hflav==5)       return it->pt();     
   }
   return 0;
 }



}
