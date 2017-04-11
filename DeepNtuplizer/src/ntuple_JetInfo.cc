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
    addBranch(tree,"event_no"    ,&event_no_    ,"event_no/i"    );
    addBranch(tree,"jet_no"    ,&jet_no_    ,"jet_no/i"    );


    // truth labels
    addBranch(tree,"gen_pt"    ,&gen_pt_    ,"gen_pt_/f"    );
    addBranch(tree,"Delta_gen_pt"    ,&Delta_gen_pt_,"Delta_gen_pt_/f"    );
    addBranch(tree,"isB",&isB_, "isB_/i");
    addBranch(tree,"isBB",&isBB_, "isBB_/i");
    addBranch(tree,"isLeptonicB",&isLeptonicB_, "isLeptonicB_/i");
    addBranch(tree,"isLeptonicB_C",&isLeptonicB_C_, "isLeptonicB_C_/i");
    addBranch(tree,"isC",&isC_, "isC_/i");
    addBranch(tree,"isUD",&isUD_, "isUD_/i");
    addBranch(tree,"isS",&isS_, "isS_/i");
    addBranch(tree,"isG",&isG_, "isG_/i");
    addBranch(tree,"isUndefined",&isUndefined_, "isUndefined_/i");

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    addBranch(tree,"isPhysB",&isPhysB_, "isPhysB_/i");
    addBranch(tree,"isPhysBB",&isPhysBB_, "isPhysBB_/i");
    addBranch(tree,"isPhysLeptonicB",&isPhysLeptonicB_, "isPhysLeptonicB_/i");
    addBranch(tree,"isPhysLeptonicB_C",&isPhysLeptonicB_C_, "isPhysLeptonicB_C_/i");
    addBranch(tree,"isPhysC",&isPhysC_, "isPhysC_/i");
    addBranch(tree,"isPhysUD",&isPhysUD_, "isPhysUD_/i");
    addBranch(tree,"isPhysS",&isPhysS_, "isPhysS_/i");
    addBranch(tree,"isPhysG",&isPhysG_, "isPhysG_/i");
    addBranch(tree,"isPhysUndefined",&isPhysUndefined_, "isPhysUndefined_/i");

    // jet variables
    //b=tree->Branch("jet_pt", &jet_pt_);
    addBranch(tree,"jet_pt", &jet_pt_);

    addBranch(tree,"jet_corr_pt", &jet_corr_pt_);
    addBranch(tree,"jet_eta", &jet_eta_);


    //jet id
    addBranch(tree,"jet_looseId", &jet_looseId_);

    // quark gluon
    addBranch(tree,"jet_qgl",   &jet_qgl_);  // qg tagger from jmar
    addBranch(tree,"QG_ptD",   &QG_ptD_);   // momentum fraction per jet constituent
    addBranch(tree,"QG_axis2", &QG_axis2_); // jet shape i.e. gluon are wider than quarks
    addBranch(tree,"QG_mult",  &QG_mult_);  // multiplicity i.e. total num of PFcands reconstructed
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

    neutrinosLepB.clear();
    neutrinosLepB_C.clear();

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
    }
    //technically a branch fill but per event, therefore here
    event_no_=iEvent.id().event();
}

//use either of these functions

bool ntuple_JetInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
    if(!coll)
        throw std::runtime_error("ntuple_JetInfo::fillBranches: no jet collection");

    /// cuts ///

    // some cuts to contrin training region
    if ( jet.pt() < jetPtMin_ ||  jet.pt() > jetPtMax_ ) return false;                  // apply jet pT cut
    if ( jet.eta() < fabs(jetAbsEtaMin_) ||jet.eta() > fabs(jetAbsEtaMax_) ) return false; // apply jet eta cut

    // often we have way to many gluons that we do not need. This randomply reduces the gluons
    if (gluonReduction_>0 && jet.partonFlavour()==21)
        if(TRandom_.Uniform()>gluonReduction_) return false;

    if(jet.genJet()==NULL)return false;


    //branch fills
    for(auto& entry : discriminators_) {
        entry.second = catchInfs(jet.bDiscriminator(entry.first),-0.1);
    }

    npv_ = vertices()->size();
    jet_no_=jetidx;

    const auto jetRef = reco::CandidatePtr(coll->ptrs().at( jetidx));
    jet_qgl_ = (*qglHandle)[jetRef];
    QG_ptD_ = (*ptDHandle)[jetRef];
    QG_axis2_ = (*axis2Handle)[jetRef];
    QG_mult_ = (*multHandle)[jetRef];

    //std::vector<Ptr<pat::Jet> > p= coll->ptrs();

    isB_=0; isBB_=0; isC_=0; isUD_=0; isS_=0; isG_=0, isLeptonicB_=0, isLeptonicB_C_=0, isUndefined_=0;
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

    switch(deep_ntuples::jet_flavour(jet, neutrinosLepB, neutrinosLepB_C)) {
    case deep_ntuples::JetFlavor::UD: isUD_=1; break;
    case deep_ntuples::JetFlavor::S:  isS_=1; break;
    case deep_ntuples::JetFlavor::B:  isB_=1; break;
    case deep_ntuples::JetFlavor::BB: isBB_=1; break;
    case deep_ntuples::JetFlavor::C:  isC_=1; break;
    case deep_ntuples::JetFlavor::G:  isG_=1; break;
    case deep_ntuples::JetFlavor::LeptonicB: isLeptonicB_=1; break;                                 
    case deep_ntuples::JetFlavor::LeptonicB_C: isLeptonicB_C_=1; break; 
    default : isUndefined_=1; break;
    }

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    isPhysB_=0; isPhysBB_=0; isPhysC_=0; isPhysUD_=0; isPhysS_=0; isPhysG_=0, isPhysLeptonicB_=0, isPhysLeptonicB_C_=0, isPhysUndefined_=0;
    switch(deep_ntuples::jet_flavour(jet, neutrinosLepB, neutrinosLepB_C, true)) {
    case deep_ntuples::JetFlavor::UD: isPhysUD_=1; break;
    case deep_ntuples::JetFlavor::S:  isPhysS_=1; break;
    case deep_ntuples::JetFlavor::B:  isPhysB_=1; break;
    case deep_ntuples::JetFlavor::BB: isPhysBB_=1; break;
    case deep_ntuples::JetFlavor::C:  isPhysC_=1; break;
    case deep_ntuples::JetFlavor::G:  isPhysG_=1; break;
    case deep_ntuples::JetFlavor::LeptonicB: isPhysLeptonicB_=1; break;                                 
    case deep_ntuples::JetFlavor::LeptonicB_C: isPhysLeptonicB_C_=1; break; 
    default : isPhysUndefined_=1; break;
    }

    if(isUndefined_ && isPhysUndefined_) return false; //skip event, if neither standard flavor definition nor physics definition fallback define a "proper flavor"

    pat::JetCollection h;

    jet_pt_ = jet.correctedJet("Uncorrected").pt();
    jet_eta_ = jet.eta();
    jet_corr_pt_ = jet.pt();


    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
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


    return true;
}
