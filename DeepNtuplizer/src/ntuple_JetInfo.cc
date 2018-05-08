/*
 * ntuple_JetInfo.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */




#include "../interface/ntuple_JetInfo.h"
#include "../interface/helpers.h"
//#include "../interface/leptonsInJet.h"
#include <vector>
#include <algorithm>
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

void ntuple_JetInfo::getInput(const edm::ParameterSet& iConfig){

    gluonReduction_=(iConfig.getParameter<double>("gluonReduction"));
    jetPtMin_=(iConfig.getParameter<double>("jetPtMin"));
    jetPtMax_=(iConfig.getParameter<double>("jetPtMax"));
    jetAbsEtaMin_=(iConfig.getParameter<double>("jetAbsEtaMin"));
    jetAbsEtaMax_=(iConfig.getParameter<double>("jetAbsEtaMax"));

    vector<string> disc_names = iConfig.getParameter<vector<string> >("bDiscriminators");
    for(auto& name : disc_names) {
        discriminators_[name] = 0.;
    }
}
void ntuple_JetInfo::initBranches(TTree* tree){

    //more general event info, here applied per jet
    addBranch(tree,"npv"    ,&npv_    ,"npv/f"    );
    addBranch(tree,"rho", &rho_, "rho/f");
    addBranch(tree,"ntrueInt",&ntrueInt_,"ntrueInt/f");
    addBranch(tree,"event_no"    ,&event_no_    ,"event_no/i"    );
    addBranch(tree,"jet_no"    ,&jet_no_    ,"jet_no/i"    );


    // truth labels
    addBranch(tree,"gen_pt"    ,&gen_pt_    ,"gen_pt_/f"    );
    addBranch(tree,"Delta_gen_pt"    ,&Delta_gen_pt_,"Delta_gen_pt_/f"    );
    addBranch(tree,"isB",&isB_, "isB_/i");
    addBranch(tree,"isGBB",&isGBB_, "isGBB_/i");
    addBranch(tree,"isBB",&isBB_, "isBB_/i");
    addBranch(tree,"isLeptonicB",&isLeptonicB_, "isLeptonicB_/i");
    addBranch(tree,"isLeptonicB_C",&isLeptonicB_C_, "isLeptonicB_C_/i");
    addBranch(tree,"isC",&isC_, "isC_/i");
    addBranch(tree,"isGCC",&isGCC_, "isGCC_/i");
    addBranch(tree,"isCC",&isCC_, "isCC_/i");
//    addBranch(tree,"isTau",&isTau_, "isTau_/i");
    addBranch(tree,"isUD",&isUD_, "isUD_/i");
    addBranch(tree,"isS",&isS_, "isS_/i");
    addBranch(tree,"isG",&isG_, "isG_/i");
    addBranch(tree,"isUndefined",&isUndefined_, "isUndefined_/i");
    addBranch(tree,"genDecay",&genDecay_, "genDecay_/f"); //dxy corresponds to the distance the Bhadron traveled

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    addBranch(tree,"isPhysB",&isPhysB_, "isPhysB_/i");
    addBranch(tree,"isPhysGBB",&isPhysGBB_, "isPhysGBB_/i");
    addBranch(tree,"isPhysBB",&isPhysBB_, "isPhysBB_/i");
    addBranch(tree,"isPhysLeptonicB",&isPhysLeptonicB_, "isPhysLeptonicB_/i");
    addBranch(tree,"isPhysLeptonicB_C",&isPhysLeptonicB_C_, "isPhysLeptonicB_C_/i");
    addBranch(tree,"isPhysC",&isPhysC_, "isPhysC_/i");
    addBranch(tree,"isPhysGCC",&isPhysGCC_, "isPhysGCC_/i");
    addBranch(tree,"isPhysCC",&isPhysCC_, "isPhysCC_/i");
//    addBranch(tree,"isPhysTau",&isPhysTau_, "isPhysTau_/i");
    addBranch(tree,"isPhysUD",&isPhysUD_, "isPhysUD_/i");
    addBranch(tree,"isPhysS",&isPhysS_, "isPhysS_/i");
    addBranch(tree,"isPhysG",&isPhysG_, "isPhysG_/i");
    addBranch(tree,"isPhysUndefined",&isPhysUndefined_, "isPhysUndefined_/i");

    // jet variables
    //b=tree->Branch("jet_pt", &jet_pt_);
    addBranch(tree,"jet_pt", &jet_pt_);
    addBranch(tree,"jet_corr_pt", &jet_corr_pt_);
    addBranch(tree,"jet_eta", &jet_eta_);
    addBranch(tree,"jet_phi", &jet_phi_);
    addBranch(tree,"jet_mass", &jet_mass_);
    addBranch(tree,"jet_energy", &jet_energy_);


    //jet id
    addBranch(tree,"jet_looseId", &jet_looseId_);

    // quark gluon
    addBranch(tree,"jet_qgl",   &jet_qgl_);  // qg tagger from jmar
    addBranch(tree,"QG_ptD",   &QG_ptD_);   // momentum fraction per jet constituent
    addBranch(tree,"QG_axis2", &QG_axis2_); // jet shape i.e. gluon are wider than quarks
    addBranch(tree,"QG_mult",  &QG_mult_);  // multiplicity i.e. total num of PFcands reconstructed

    // yutas quark-gluon info
    addBranch(tree,"y_multiplicity"    ,&y_multiplicity_,"y_multiplicity_/f"    );
    addBranch(tree,"y_charged_multiplicity"    ,&y_charged_multiplicity_,"y_charged_multiplicity_/f"    );
    addBranch(tree,"y_neutral_multiplicity"    ,&y_neutral_multiplicity_,"y_neutral_multiplicity_/f"    );
    addBranch(tree,"y_ptD"    ,&y_ptD_,"y_ptD_/f"    );
    addBranch(tree,"y_axis1"    ,&y_axis1_,"y_axis1_/f"    );
    addBranch(tree,"y_axis2"    ,&y_axis2_,"y_axis2_/f"    );
    addBranch(tree,"y_pt_dr_log"    ,&y_pt_dr_log_,"y_pt_dr_log_/f"    );


    // in the jet

    addBranch(tree,"muons_number", &muons_number_, "muons_number_/i");
    addBranch(tree,"electrons_number", &electrons_number_, "electrons_number_/i");

    addBranch(tree,"muons_isLooseMuon", &muons_isLooseMuon_, "muons_isLooseMuon_[muons_number_]/i");
    addBranch(tree,"muons_isTightMuon", &muons_isTightMuon_, "muons_isTightMuon_[muons_number_]/i");
    addBranch(tree,"muons_isSoftMuon", &muons_isSoftMuon_, "muons_isSoftMuon_[muons_number_]/i");
    addBranch(tree,"muons_isHighPtMuon", &muons_isHighPtMuon_, "muons_isHighPtMuon_[muons_number_]/i");
    addBranch(tree,"muons_pt", &muons_pt_, "muons_pt_[muons_number_]/f");
    addBranch(tree,"muons_relEta", &muons_relEta_, "muons_relEta_[muons_number_]/f");
    addBranch(tree,"muons_relPhi", &muons_relPhi_, "muons_relPhi_[muons_number_]/f");
    addBranch(tree,"muons_energy", &muons_energy_, "muons_energy_[muons_number_]/f");
    addBranch(tree,"electrons_pt", &electrons_pt_, "electrons_pt_[electrons_number_]/f");
    addBranch(tree,"electrons_relEta", &electrons_relEta_, "electrons_relEta_[electrons_number_]/f");
    addBranch(tree,"electrons_relPhi", &electrons_relPhi_, "electrons_relPhi_[electrons_number_]/f");
    addBranch(tree,"electrons_energy", &electrons_energy_, "electrons_energy_[electrons_number_]/f");


    addBranch(tree,"gen_pt_Recluster"    ,&gen_pt_Recluster_    ,"gen_pt_Recluster_/f"    );
    addBranch(tree,"gen_pt_WithNu"    ,&gen_pt_WithNu_    ,"gen_pt_WithNu_/f"    );
    addBranch(tree,"Delta_gen_pt_Recluster"    ,&Delta_gen_pt_Recluster_    ,"Delta_gen_pt_Recluster_/f"    );
    addBranch(tree,"Delta_gen_pt_WithNu"    ,&Delta_gen_pt_WithNu_    ,"Delta_gen_pt_WithNu_/f"    );

    if(1) // discriminators might need to be filled differently. FIXME
        for(auto& entry : discriminators_) {
            string better_name(entry.first);
            std::replace(better_name.begin(), better_name.end(), ':', '_');
            addBranch(tree,better_name.c_str(), &entry.second, (better_name+"/F").c_str());
        }
}


void ntuple_JetInfo::readEvent(const edm::Event& iEvent){

    iEvent.getByToken(qglToken_, qglHandle);
    iEvent.getByToken(ptDToken_, ptDHandle);
    iEvent.getByToken(axis2Token_, axis2Handle);
    iEvent.getByToken(multToken_, multHandle);

    iEvent.getByToken(genJetMatchReclusterToken_, genJetMatchRecluster);
    iEvent.getByToken(genJetMatchWithNuToken_, genJetMatchWithNu);

    iEvent.getByToken(genParticlesToken_, genParticlesHandle);


    iEvent.getByToken(muonsToken_, muonsHandle);
    iEvent.getByToken(electronsToken_, electronsHandle);

    event_no_=iEvent.id().event();

    //presumably this whole part can be removed!


    neutrinosLepB.clear();
    neutrinosLepB_C.clear();
    gToBB.clear();
    gToCC.clear();
    alltaus_.clear();
    Bhadron_.clear();
    Bhadron_daughter_.clear();

    //std::cout << " start search for a b in this event "<<std::endl;
 for (const reco::Candidate &genC : *genParticlesHandle)
   {
     const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
     
     if((abs(gen.pdgId())>500&&abs(gen.pdgId())<600)||(abs(gen.pdgId())>5000&&abs(gen.pdgId())<6000)) {

       //std::cout<<gen.end_vertex()<<endl;

       Bhadron_.push_back(gen);
       if(gen.numberOfDaughters()>0){
     
	 if( (abs(gen.daughter(0)->pdgId())>500&&abs(gen.daughter(0)->pdgId())<600)||(abs(gen.daughter(0)->pdgId())>5000&&abs(gen.daughter(0)->pdgId())<6000))
	   {
	     if(gen.daughter(0)->numberOfDaughters()>0)
	       {
		
		 const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*(gen.daughter(0)->daughter(0)));
		 
		 if(daughter_.vx()!=gen.vx())
		   { 
		     Bhadron_daughter_.push_back(daughter_);
		   }
                 else Bhadron_daughter_.push_back(gen);
		 //  std::cout << "only b daughters " << endl;
		 // }
	       }
	     else  Bhadron_daughter_.push_back(gen);
	     
	   }
	 else{
	   //  std::cout<<gen.daughter(0)->vx()<< " oh  " <<gen.vx()<<" "<<gen.pt() <<" "<<  gen.daughter(0)->pdgId() <<std::endl; 
	  
	   const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*gen.daughter(0));
	   Bhadron_daughter_.push_back(daughter_);
	 }

       }// if daughter is there
       else {
	 
	 //std::cout << " lonly B hadron, has NO daughter??? "<<std::endl;
	 Bhadron_daughter_.push_back(gen);
       }
     }
   }

 for (const reco::Candidate &genC : *genParticlesHandle) {
        const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
        if(abs(gen.pdgId())==12||abs(gen.pdgId())==14||abs(gen.pdgId())==16) {
            const reco::GenParticle* mother =  static_cast< const reco::GenParticle*> (gen.mother());
            if(mother!=NULL) {
                if((abs(mother->pdgId())>500&&abs(mother->pdgId())<600)||(abs(mother->pdgId())>5000&&abs(mother->pdgId())<6000)) {
                    neutrinosLepB.emplace_back(gen);
                }
                if((abs(mother->pdgId())>400&&abs(mother->pdgId())<500)||(abs(mother->pdgId())>4000&&abs(mother->pdgId())<5000)) {
                    neutrinosLepB_C.emplace_back(gen);
                }
            }
            else {
                std::cout << "No mother" << std::endl;
            }
        }

        int id(std::abs(gen.pdgId())); 
        int status(gen.status());

        if (id == 21 && status >= 21 && status <= 59) { //// Pythia8 hard scatter, ISR, or FSR
            if ( gen.numberOfDaughters() == 2 ) {
                const reco::Candidate* d0 = gen.daughter(0);
                const reco::Candidate* d1 = gen.daughter(1);
                if ( std::abs(d0->pdgId()) == 5 && std::abs(d1->pdgId()) == 5
                        && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4) gToBB.push_back(gen) ;
                if ( std::abs(d0->pdgId()) == 4 && std::abs(d1->pdgId()) == 4
                        && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4) gToCC.push_back(gen) ;
            }
        }

        if(id == 15 && false){
            alltaus_.push_back(gen);
        }

    }
    //technically a branch fill but per event, therefore here
}

//use either of these functions

bool ntuple_JetInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
    if(!coll)
        throw std::runtime_error("ntuple_JetInfo::fillBranches: no jet collection");

    /// cuts ///
    bool returnval=true;

    // some cuts to contrin training region
    if ( jet.pt() < jetPtMin_ ||  jet.pt() > jetPtMax_ ) returnval=false;                  // apply jet pT cut
    if ( fabs(jet.eta()) < jetAbsEtaMin_ || fabs(jet.eta()) > jetAbsEtaMax_ ) returnval=false; // apply jet eta cut


    // often we have way to many gluons that we do not need. This randomply reduces the gluons
    if (gluonReduction_>0 && jet.partonFlavour()==21)
        if(TRandom_.Uniform()>gluonReduction_) returnval=false;

    if(jet.genJet()==NULL)returnval=false;


    //branch fills
    for(auto& entry : discriminators_) {
        entry.second = catchInfs(jet.bDiscriminator(entry.first),-0.1);
    }

    npv_ = vertices()->size();

    for (auto const& v : *pupInfo()) {
        int bx = v.getBunchCrossing();
        if (bx == 0) {
            ntrueInt_ = v.getTrueNumInteractions();
        }
    }
    rho_ = rhoInfo()[0];


    jet_no_=jetidx;

    const auto jetRef = reco::CandidatePtr(coll->ptrs().at( jetidx));

    jet_qgl_ = (*qglHandle)[jetRef];
    QG_ptD_ = (*ptDHandle)[jetRef];
    QG_axis2_ = (*axis2Handle)[jetRef];
    QG_mult_ = (*multHandle)[jetRef];


    //std::vector<Ptr<pat::Jet> > p= coll->ptrs();

    isB_=0; isGBB_=0; isBB_=0; isC_=0; isGCC_=0; isCC_=0; isUD_=0;isTau_=0;
    isS_=0; isG_=0, isLeptonicB_=0, isLeptonicB_C_=0, isUndefined_=0;
    auto muIds = deep_ntuples::jet_muonsIds(jet,*muonsHandle);
    auto elecIds = deep_ntuples::jet_electronsIds(jet,*electronsHandle);

    muons_number_ = muIds.size();
    electrons_number_ = elecIds.size();

    float etasign = 1.;
    if (jet.eta()<0) etasign = -1.;

    for(std::size_t i=0; i<max_num_lept; i++) {
        if (i < muIds.size()) {
            const auto & muon = (*muonsHandle).at(muIds.at(i));
            muons_isLooseMuon_[i] = muon.isLooseMuon();
            muons_isTightMuon_[i] = muon.isTightMuon(vertices()->at(0));
            muons_isSoftMuon_[i] = muon.isSoftMuon(vertices()->at(0));
            muons_isHighPtMuon_[i] = muon.isHighPtMuon(vertices()->at(0));
            muons_pt_[i] = muon.pt();
            muons_relEta_[i] = etasign*(muon.eta()-jet.eta());
            muons_relPhi_[i] = reco::deltaPhi(muon.phi(),jet.phi());
            muons_energy_[i] = muon.energy()/jet.energy();
        }
        if (i < elecIds.size()) {
            const auto & electron = (*electronsHandle).at(elecIds.at(i));
            electrons_pt_[i] = electron.pt();
            electrons_relEta_[i] = etasign*(electron.eta()-jet.eta());
            electrons_relPhi_[i] = reco::deltaPhi(electron.phi(),jet.phi());
            electrons_energy_[i] = electron.energy()/jet.energy();
        }
    }

    //// Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
    if(jet.genJet()!=NULL){
        switch(deep_ntuples::jet_flavour(jet, gToBB, gToCC, neutrinosLepB, neutrinosLepB_C, alltaus_)) {
        case deep_ntuples::JetFlavor::B:  isB_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB: isLeptonicB_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB_C: isLeptonicB_C_=1; break;
        case deep_ntuples::JetFlavor::GBB: isGBB_=1; break;
        case deep_ntuples::JetFlavor::BB: isBB_=1; break;
        case deep_ntuples::JetFlavor::C:  isC_=1; break;
        case deep_ntuples::JetFlavor::GCC: isGCC_=1; break;
        case deep_ntuples::JetFlavor::CC: isCC_=1; break;
        case deep_ntuples::JetFlavor::TAU: isTau_=1;break;
        case deep_ntuples::JetFlavor::G:  isG_=1; break;
        case deep_ntuples::JetFlavor::UD: isUD_=1; break;
        case deep_ntuples::JetFlavor::S:  isS_=1; break;
        default : isUndefined_=1; break;
        }
    }

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    //// Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
    isPhysB_=0; isPhysBB_=0; isPhysGBB_=0; isPhysC_=0; isPhysCC_=0;
    isPhysGCC_=0; isPhysUD_=0; isPhysS_=0; isPhysG_=0, isPhysLeptonicB_=0, isPhysLeptonicB_C_=0, isPhysUndefined_=0;
    isPhysTau_=0;
    if(jet.genJet()!=NULL){
        switch(deep_ntuples::jet_flavour(jet, gToBB, gToCC, neutrinosLepB, neutrinosLepB_C, alltaus_,true)) {
        case deep_ntuples::JetFlavor::UD: isPhysUD_=1; break;
        case deep_ntuples::JetFlavor::S:  isPhysS_=1; break;
        case deep_ntuples::JetFlavor::B:  isPhysB_=1; break;
        case deep_ntuples::JetFlavor::BB: isPhysBB_=1; break;
        case deep_ntuples::JetFlavor::GBB: isPhysGBB_=1; break;
        case deep_ntuples::JetFlavor::C:  isPhysC_=1; break;
        case deep_ntuples::JetFlavor::CC: isPhysCC_=1; break;
        case deep_ntuples::JetFlavor::GCC: isPhysGCC_=1; break;
        case deep_ntuples::JetFlavor::TAU: isPhysTau_=1;break;
        case deep_ntuples::JetFlavor::G:  isPhysG_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB: isPhysLeptonicB_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB_C: isPhysLeptonicB_C_=1; break;
        default : isPhysUndefined_=1; break;
        }
    }

    if(!jet.genJet()){//for data
        isUndefined_=1;isPhysUndefined_=1;
    }

    if(isUndefined_ && isPhysUndefined_) returnval=false; //skip event, if neither standard flavor definition nor physics definition fallback define a "proper flavor"

    pat::JetCollection h;

    jet_pt_ = jet.correctedJet("Uncorrected").pt();
    jet_eta_ = jet.eta();
    jet_phi_ = jet.phi();
    jet_corr_pt_ = jet.pt();
    jet_mass_ = jet.mass();
    jet_energy_ = jet.energy();

    genDecay_ = -1.;

    reco::GenParticleRefVector Bhadrons_in_jet = jet.jetFlavourInfo().getbHadrons();

    if (Bhadrons_in_jet.size() > 0){ 

        for (unsigned int idx=0; idx<Bhadron_.size(); ++idx){

            reco::GenParticle bhad = Bhadron_[idx];

            bool bhad_is_in_jet = false;

            for (reco::GenParticleRefVector::const_iterator bhad_in_jet = Bhadrons_in_jet.begin(); bhad_in_jet!=Bhadrons_in_jet.end(); ++bhad_in_jet) {

                //check if bhad is identical to bhad_in_jet
                if ( (*bhad_in_jet)->pt() == bhad.pt() && (*bhad_in_jet)->eta() == bhad.eta()
                        && (*bhad_in_jet)->phi() == bhad.phi() && (*bhad_in_jet)->pdgId() == bhad.pdgId())              
                    bhad_is_in_jet = true;
            }
            if (bhad_is_in_jet){

                if (Bhadron_daughter_[idx].vx()!=bhad.vx()){
                    
                    float vx = Bhadron_daughter_[idx].vx() - bhad.vx();
                    float vy = Bhadron_daughter_[idx].vy() - bhad.vy();

                    float dxy = sqrt(vx*vx+vy*vy);
                    if (dxy > genDecay_)
                        genDecay_= dxy;
                }
                else if (genDecay_ < 0) 
                    genDecay_ = -0.1;
            }
        }
    }



    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    try{
        float NHF  = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF  = jet.chargedHadronEnergyFraction();
        //float MUF  = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
        float NumNeutralParticles =jet.neutralMultiplicity();
        float CHM      = jet.chargedMultiplicity();

        jet_looseId_ = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(jet_eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jet_eta_)>2.4) && abs(jet_eta_)<=2.7) ||
                (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(jet_eta_)>2.7 && abs(jet_eta_)<=3.0 ) ||
                (NEMF<0.90 && NumNeutralParticles>10 && abs(jet_eta_)>3.0 );
    }catch(const cms::Exception &e){
        jet_looseId_ = 1;
    }


    gen_pt_ =  0;
    Delta_gen_pt_ =  0;
    gen_pt_Recluster_=0;
    gen_pt_WithNu_=0;
    Delta_gen_pt_Recluster_=0;
    Delta_gen_pt_WithNu_=0;

    if(jet.genJet()){
        gen_pt_ =  jet.genJet()->pt();
        Delta_gen_pt_ =  jet.genJet()->pt()- jet_pt_;

        const edm::RefToBase<pat::Jet> patJetRef = coll->refAt(jetidx);
        reco::GenJetRef genjetRecluster = (*genJetMatchRecluster)[patJetRef];

        gen_pt_Recluster_ = 0.;
        if (genjetRecluster.isNonnull() && genjetRecluster.isAvailable()) {
            gen_pt_Recluster_ = genjetRecluster->pt();
        }
        reco::GenJetRef genjetWithNu = (*genJetMatchWithNu)[patJetRef];

        gen_pt_WithNu_ = 0.;
        if (genjetWithNu.isNonnull() && genjetWithNu.isAvailable()) {
            gen_pt_WithNu_ = genjetWithNu->pt();
        }

        Delta_gen_pt_Recluster_=gen_pt_Recluster_-jet.pt();
        Delta_gen_pt_WithNu_=gen_pt_WithNu_-jet.pt();
    }


    auto qgtuple=yuta::calcVariables(&jet);
    //(multiplicity, charged_multiplicity, neutral_multiplicity, ptD, axis1, axis2, pt_dr_log);

    y_multiplicity_=std::get<0>(qgtuple);
    y_charged_multiplicity_=std::get<1>(qgtuple);
    y_neutral_multiplicity_=std::get<2>(qgtuple);
    y_ptD_    =  std::get<3>(qgtuple);
    y_axis1_  =  std::get<4>(qgtuple);
    y_axis2_  =  std::get<5>(qgtuple);
    y_pt_dr_log_=std::get<6>(qgtuple);

   



    return returnval;
}
