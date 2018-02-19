/*
 * ntuple_JetInfo.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_

#include "ntuple_content.h"
#include "TRandom3.h"
#include <map>
#include <string>

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

/*
 * For global jet info such as eta, pt, gen info
 */
class ntuple_JetInfo: public ntuple_content{
public:
    ntuple_JetInfo():ntuple_content(),
    gluonReduction_(0),
    useherwcompat_matching_(false),
    isherwig_(false),
    isData_(false)
{}

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);

    void setAxis2Token(edm::EDGetTokenT<edm::ValueMap<float> > axis2Token) {
        axis2Token_ = axis2Token;
    }

    void setMultToken(edm::EDGetTokenT<edm::ValueMap<int> > multToken) {
        multToken_ = multToken;
    }

    void setPtDToken(edm::EDGetTokenT<edm::ValueMap<float> > ptDToken) {
        ptDToken_ = ptDToken;
    }

    void setQglToken(edm::EDGetTokenT<edm::ValueMap<float> > qglToken) {
        qglToken_ = qglToken;
    }

    void setGenJetMatchReclusterToken(
            edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchReclusterToken) {
        genJetMatchReclusterToken_ = genJetMatchReclusterToken;
    }

    void setGenJetMatchWithNuToken(
            edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchWithNuToken) {
        genJetMatchWithNuToken_ = genJetMatchWithNuToken;
    }

    void setGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken) {
        genParticlesToken_ = genParticlesToken;
    }

    void setMuonsToken(edm::EDGetTokenT<pat::MuonCollection> muonsToken) {
        muonsToken_ = muonsToken;
    }

    void setElectronsToken(edm::EDGetTokenT<pat::ElectronCollection> electronsToken) {
        electronsToken_ = electronsToken;
    }

    void setLHEToken(edm::EDGetTokenT<LHEEventProduct> lheToken) {
        lheToken_ = lheToken;
    }

    void setUseHerwigCompatibleMatching(const bool use){
        useherwcompat_matching_=use;
    }
    void setIsHerwig(const bool use){
        isherwig_=use;
    }
    void setIsData(const bool use){
        isData_=use;
    }

    //private:

    double                    jetPtMin_;
    double                    jetPtMax_;
    double                    jetAbsEtaMin_;
    double                    jetAbsEtaMax_;

    //Quark gluon likelihood
    edm::EDGetTokenT<edm::ValueMap<float>>   qglToken_;
    edm::EDGetTokenT<edm::ValueMap<float>>   ptDToken_;
    edm::EDGetTokenT<edm::ValueMap<float>>   axis2Token_;
    edm::EDGetTokenT<edm::ValueMap<int>>     multToken_;

    edm::Handle<edm::ValueMap<float>> qglHandle;
    edm::Handle<edm::ValueMap<float>> ptDHandle;
    edm::Handle<edm::ValueMap<float>> axis2Handle;
    edm::Handle<edm::ValueMap<int>> multHandle;


    edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchReclusterToken_;
    edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchWithNuToken_;

    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

    edm::EDGetTokenT<pat::MuonCollection> muonsToken_;       
    edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;

    edm::EDGetTokenT<LHEEventProduct> lheToken_;

    edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchRecluster;
    edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchWithNu;

    edm::Handle<reco::GenParticleCollection> genParticlesHandle;

    edm::Handle<pat::MuonCollection> muonsHandle;
    edm::Handle<pat::ElectronCollection> electronsHandle;

    edm::Handle<LHEEventProduct> lheInfo;


    TRandom3 TRandom_;
    float gluonReduction_;

    std::vector <reco::GenParticle> neutrinosLepB;
    std::vector <reco::GenParticle> neutrinosLepB_C;

    std::vector<reco::GenParticle> gToBB;
    std::vector<reco::GenParticle> gToCC;
    std::vector<reco::GenParticle> alltaus_;

    std::vector<reco::GenParticle> Bhadron_;
    std::vector<reco::GenParticle> Bhadron_daughter_;



    bool useherwcompat_matching_;
    bool isherwig_;
    bool isData_;

    /////////branches

    // labels (MC truth)
    // regressions pt, Deta, Dphi
    float gen_pt_;
    float Delta_gen_pt_;
    //classification
    int isB_;
    int isGBB_;
    int isBB_;
    int isC_;
    int isGCC_;
    int isCC_;
    int isUD_;
    int isS_;
    int isG_;
    int isUndefined_;
    float genDecay_;
    int isLeptonicB_;
    int isLeptonicB_C_;
    int isTau_;

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    int isPhysB_;
    int isPhysGBB_;
    int isPhysBB_;
    int isPhysC_;
    int isPhysGCC_;
    int isPhysCC_;
    int isPhysUD_;
    int isPhysS_;
    int isPhysG_;
    int isPhysUndefined_;
    int isPhysLeptonicB_;
    int isPhysLeptonicB_C_;
    int isPhysTau_;
    int isRealData_;

    // global variables
    float npv_;
    float ntrueInt_;
    float rho_;
    unsigned int event_no_;
    unsigned int jet_no_;

    // jet variables
    float jet_pt_;
    float jet_corr_pt_;
    float  jet_eta_;
    float  jet_phi_;
    float  jet_mass_;
    float  jet_energy_;

    // variables for event weights
    bool useLHEWeights_;
    std::vector<double> pupWeights = {0., 0.9059208563049067, 0.4689946658078284, 0.4646384395209148, 0.35969952879169353, 0.20581098696265324, 0.1272975292012222, 0.11484070280016173, 0.18997631713506571, 0.2134442228038988, 0.2643260644534947, 0.38290419689123334, 0.5233231080095909, 0.6531355038202962, 0.7413472652417477, 0.8226986246368242, 0.909166618129074, 0.9754246515040906, 1.0116680456990728, 1.0419003275306946, 1.0762295787102802, 1.1034517536607453, 1.1240057591005752, 1.1262609366258989, 1.1198925605086503, 1.1078563848872773, 1.083911370342871, 1.0496801331833463, 1.0169474777740053, 0.9813073284102061, 0.9426647112682524, 0.9233559974737714, 0.894377685151131, 0.884518640484559, 0.8804996845439541, 0.888865624721422, 0.9107720191678217, 0.953824079957025, 1.0112115792549745, 1.0947266350670766, 1.1906625100966544, 1.3370686699103784, 1.510754710039529, 1.7098870624229212, 1.9391935863130325, 2.1529480239752585, 2.396292845558282, 2.7067102383965915, 2.9081731387455574, 3.126536829527805, 3.263323444233829, 3.3766746819738014, 3.3825183331846596, 3.332122539734733, 3.309580249670872, 2.962530063365449, 2.609107759179973, 2.4216997456337266, 2.228050917243899, 1.963846134469284, 1.844961576497657, 1.431817556288543, 1.5742387193813208, 1.731514171818786, 2.21989485741883, 2.8198818386713915, 3.8909891214054464, 5.4911613943458395, 8.610687700958811, 12.634799153629881, 20.32733145948569};
    float crossSection_;
    float luminosity_;
    float efficiency_;
    float event_weight_;

    float jet_looseId_;

    // quark/gluon
    float jet_qgl_;
    float QG_ptD_;
    float QG_axis2_;
    float QG_mult_;


    float y_multiplicity_;
    float y_charged_multiplicity_;
    float y_neutral_multiplicity_;
    float y_ptD_;
    float y_axis1_;
    float y_axis2_;
    float y_pt_dr_log_;

    static constexpr std::size_t max_num_lept = 5;
    int muons_isLooseMuon_[max_num_lept];
    int muons_isTightMuon_[max_num_lept];
    int muons_isSoftMuon_[max_num_lept];
    int muons_isHighPtMuon_[max_num_lept]; 
    float muons_pt_[max_num_lept]; 
    float muons_relEta_[max_num_lept]; 
    float muons_relPhi_[max_num_lept]; 
    float muons_energy_[max_num_lept]; 
    float electrons_pt_[max_num_lept]; 
    float electrons_relEta_[max_num_lept]; 
    float electrons_relPhi_[max_num_lept]; 
    float electrons_energy_[max_num_lept];

    int muons_number_ = 0;
    int electrons_number_ = 0;

    float gen_pt_Recluster_;
    float gen_pt_WithNu_;
    float Delta_gen_pt_Recluster_;
    float Delta_gen_pt_WithNu_;
    std::map<std::string, float> discriminators_;
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_ */
