/*
 * ntuple_pfcands.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_

#include "ntuple_content.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class ntuple_pfCands: public ntuple_content{
public:

    ntuple_pfCands():ntuple_content(){}

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup);

    void setSVToken(const edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> & t){
        svToken_=t;
    }


    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);

private:
	//config info
	unsigned int ncollinear_ = 0;

	//other

    edm::ESHandle<TransientTrackBuilder> builder;

    unsigned int n_Cpfcand_;
    float nCpfcand_;

    static constexpr size_t max_pfcand_=100;

    float  Cpfcan_pt_[max_pfcand_];
    float  Cpfcan_eta_[max_pfcand_];
    float  Cpfcan_phi_[max_pfcand_];
    float  Cpfcan_ptrel_[max_pfcand_];
    float  Cpfcan_erel_[max_pfcand_];
    float  Cpfcan_phirel_[max_pfcand_];
    float  Cpfcan_etarel_[max_pfcand_];
    float  Cpfcan_deltaR_[max_pfcand_];
    float  Cpfcan_puppiw_[max_pfcand_];
    float   Cpfcan_VTX_ass_[max_pfcand_];

    float   Cpfcan_fromPV_[max_pfcand_];

    float Cpfcan_vertexChi2_[max_pfcand_];
    float Cpfcan_vertexNdof_[max_pfcand_];
    float Cpfcan_vertexNormalizedChi2_[max_pfcand_];
    float Cpfcan_vertex_rho_[max_pfcand_];
    float Cpfcan_vertex_phirel_[max_pfcand_];
    float Cpfcan_vertex_etarel_[max_pfcand_];
    float Cpfcan_vertexRef_mass_[max_pfcand_];

    // covariance
    float  Cpfcan_dz_[max_pfcand_];
    float  Cpfcan_dxy_[max_pfcand_];

    float  Cpfcan_dxyerrinv_[max_pfcand_];
    float  Cpfcan_dxysig_[max_pfcand_];

    float  Cpfcan_dptdpt_[max_pfcand_];
    float  Cpfcan_detadeta_[max_pfcand_];
    float  Cpfcan_dphidphi_[max_pfcand_];
    float  Cpfcan_dxydxy_[max_pfcand_];
    float  Cpfcan_dzdz_[max_pfcand_];
    float  Cpfcan_dxydz_[max_pfcand_];
    float  Cpfcan_dphidxy_[max_pfcand_];
    float  Cpfcan_dlambdadz_[max_pfcand_];



    float Cpfcan_BtagPf_trackMomentum_[max_pfcand_];
    float Cpfcan_BtagPf_trackEta_[max_pfcand_];
    float Cpfcan_BtagPf_trackEtaRel_[max_pfcand_];
    float Cpfcan_BtagPf_trackPtRel_[max_pfcand_];
    float Cpfcan_BtagPf_trackPPar_[max_pfcand_];
    float Cpfcan_BtagPf_trackDeltaR_[max_pfcand_];
    float Cpfcan_BtagPf_trackPtRatio_[max_pfcand_];
    float Cpfcan_BtagPf_trackPParRatio_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip3dVal_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip3dSig_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip2dVal_[max_pfcand_];
    float Cpfcan_BtagPf_trackSip2dSig_[max_pfcand_];

    float Cpfcan_BtagPf_trackDecayLen_[max_pfcand_];

    float Cpfcan_BtagPf_trackJetDistVal_[max_pfcand_];
    float Cpfcan_BtagPf_trackJetDistSig_[max_pfcand_];

    // ID, skipped "charged hadron" as that is true if now the other
    // TODO (comment of Markus Stoye) add reco information
    float Cpfcan_isMu_[max_pfcand_]; // pitty that the quality is missing
    float Cpfcan_isEl_[max_pfcand_]; // pitty that the quality is missing
    float Cpfcan_charge_[max_pfcand_];

    // track quality
    float Cpfcan_lostInnerHits_[max_pfcand_];
    float Cpfcan_chi2_[max_pfcand_];
    float Cpfcan_quality_[max_pfcand_];

    float Cpfcan_drminsv_[max_pfcand_];

    //Neutral Pf candidates
    unsigned int n_Npfcand_;
    float nNpfcand_;
    float  Npfcan_pt_[max_pfcand_];
    float  Npfcan_eta_[max_pfcand_];
    float  Npfcan_phi_[max_pfcand_];
    float  Npfcan_ptrel_[max_pfcand_];
    float  Npfcan_erel_[max_pfcand_];
    float  Npfcan_puppiw_[max_pfcand_];
    float  Npfcan_phirel_[max_pfcand_];
    float  Npfcan_etarel_[max_pfcand_];
    float  Npfcan_deltaR_[max_pfcand_];
    float  Npfcan_isGamma_[max_pfcand_];
    float  Npfcan_HadFrac_[max_pfcand_];
    float  Npfcan_drminsv_[max_pfcand_];

    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
    edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;

    float mindrsvpfcand(const std::vector<reco::VertexCompositePtrCandidate> svs, const pat::PackedCandidate* pfcand);

};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PFCANDS_H_ */
