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

/*
 * For global jet info such as eta, pt, gen info
 */
class ntuple_JetInfo: public ntuple_content{
public:
    ntuple_JetInfo():ntuple_content(),gluonReduction_(0){}

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

    void setGenParticlesToken(
            edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken) {
        genParticlesToken_ = genParticlesToken;
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

    edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchRecluster;
    edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchWithNu;

    edm::Handle<reco::GenParticleCollection> genParticlesHandle;

    TRandom3 TRandom_;
    float gluonReduction_;

    std::vector < reco::GenParticle> neutrinosLepB;
    std::vector < reco::GenParticle> neutrinosLepB_C;

    // labels (MC truth)
    // regressions pt, Deta, Dphi
    float gen_pt_;
    float Delta_gen_pt_;
    //classification
    int isB_;
    int isBB_;
    int isC_;
    int isUD_;
    int isS_;
    int isG_;
    int isUndefined_;
    int isLeptonicB_;
    int isLeptonicB_C_;

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    int isPhysB_;
    int isPhysBB_;
    int isPhysC_;
    int isPhysUD_;
    int isPhysS_;
    int isPhysG_;
    int isPhysUndefined_;
    int isPhysLeptonicB_;
    int isPhysLeptonicB_C_;

    // global variables
    float npv_;
    unsigned int event_no_;
    unsigned int jet_no_;

    // jet variables
    float jet_pt_;
    float jet_corr_pt_;
    float  jet_eta_;

    float jet_looseId_;
    // quark/gluon
    float jet_qgl_;
    float QG_ptD_;
    float QG_axis2_;
    float QG_mult_;

    float gen_pt_Recluster_;
    float gen_pt_WithNu_;
    float Delta_gen_pt_Recluster_;
    float Delta_gen_pt_WithNu_;
    std::map<std::string, float> discriminators_;
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_ */
