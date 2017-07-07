/*
 * ntuple_pfCands.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_pfCands.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "../interface/sorting_modules.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TVector3.h"
//should set seed?
#include "TRandom3.h"

class TrackInfoBuilder{
public:
    TrackInfoBuilder(edm::ESHandle<TransientTrackBuilder> & build):
        builder(build),
        trackMomentum_(0),
        trackEta_(0),
        trackEtaRel_(0),
        trackPtRel_(0),
        trackPPar_(0),
        trackDeltaR_(0),
        trackPtRatio_(0),
        trackPParRatio_(0),
        trackSip2dVal_(0),
        trackSip2dSig_(0),
        trackSip3dVal_(0),
        trackSip3dSig_(0),

        trackJetDistVal_(0),
        trackJetDistSig_(0)
{


}

    void buildTrackInfo(const pat::PackedCandidate* PackedCandidate_ ,const math::XYZVector&  jetDir, GlobalVector refjetdirection, const reco::Vertex & pv){
        const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();
        //                       const reco::Track * pf_track=PackedCandidate_->bestTrack();
        //edm::ESHandle<TransientTrackBuilder> builder;
        //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
        reco::TransientTrack transientTrack;
        transientTrack=builder->build(PseudoTrack);
        Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, refjetdirection, pv).second;
        Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, refjetdirection, pv).second;
        Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, refjetdirection, pv).second;
        math::XYZVector trackMom = PseudoTrack.momentum();
        double trackMag = std::sqrt(trackMom.Mag2());
        TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());
        TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());


        trackMomentum_=std::sqrt(trackMom.Mag2());
        trackEta_= trackMom.Eta();
        trackEtaRel_=reco::btau::etaRel(jetDir, trackMom);
        trackPtRel_=trackMom3.Perp(jetDir3);
        trackPPar_=jetDir.Dot(trackMom);
        trackDeltaR_=reco::deltaR(trackMom, jetDir);
        trackPtRatio_=trackMom3.Perp(jetDir3) / trackMag;
        trackPParRatio_=jetDir.Dot(trackMom) / trackMag;
        trackSip2dVal_=static_cast<float>(meas_ip2d.value());

        trackSip2dSig_=static_cast<float>(meas_ip2d.significance());
        trackSip3dVal_=static_cast<float>(meas_ip3d.value());


        trackSip3dSig_=static_cast<float>(meas_ip3d.significance());
        trackJetDistVal_=static_cast<float>(jetdist.value());
        trackJetDistSig_=static_cast<float>(jetdist.significance());

    }

    const float& getTrackDeltaR() const {return trackDeltaR_;}
    const float& getTrackEta() const {return trackEta_;}
    const float& getTrackEtaRel() const {return trackEtaRel_;}
    const float& getTrackJetDistSig() const {return trackJetDistSig_;}
    const float& getTrackJetDistVal() const {return trackJetDistVal_;}
    const float& getTrackMomentum() const {return trackMomentum_;}
    const float& getTrackPPar() const {return trackPPar_;}
    const float& getTrackPParRatio() const {return trackPParRatio_;}
    const float& getTrackPtRatio() const {return trackPtRatio_;}
    const float& getTrackPtRel() const {return trackPtRel_;}
    const float& getTrackSip2dSig() const {return trackSip2dSig_;}
    const float& getTrackSip2dVal() const {return trackSip2dVal_;}
    const float& getTrackSip3dSig() const {return trackSip3dSig_;}
    const float& getTrackSip3dVal() const {return trackSip3dVal_;}

private:

    edm::ESHandle<TransientTrackBuilder> builder;

    float trackMomentum_;
    float trackEta_;
    float trackEtaRel_;
    float trackPtRel_;
    float trackPPar_;
    float trackDeltaR_;
    float trackPtRatio_;
    float trackPParRatio_;
    float trackSip2dVal_;
    float trackSip2dSig_;
    float trackSip3dVal_;
    float trackSip3dSig_;

    float trackJetDistVal_;
    float trackJetDistSig_;

};




void ntuple_pfCands::readSetup(const edm::EventSetup& iSetup){

    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

}

void ntuple_pfCands::getInput(const edm::ParameterSet& iConfig){
	ncollinear_ = (iConfig.getParameter<unsigned int>("nCollinear"));
}

void ntuple_pfCands::initBranches(TTree* tree){

    addBranch(tree,"n_Cpfcand", &n_Cpfcand_,"n_Cpfcand_/i");

    addBranch(tree,"nCpfcand", &nCpfcand_,"nCpfcand_/f");

    addBranch(tree,"Cpfcan_pt", &Cpfcan_pt_,"Cpfcan_pt_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_eta", &Cpfcan_eta_,"Cpfcan_eta_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_phi", &Cpfcan_phi_,"Cpfcan_phi_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_ptrel", &Cpfcan_ptrel_,"Cpfcan_ptrel_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_erel", &Cpfcan_erel_,"Cpfcan_erel_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_phirel",&Cpfcan_phirel_,"Cpfcan_phirel_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_etarel",&Cpfcan_etarel_,"Cpfcan_etarel_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_deltaR",&Cpfcan_deltaR_,"Cpfcan_deltaR_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_puppiw",&Cpfcan_puppiw_,"Cpfcan_puppiw_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_dxy",&Cpfcan_dxy_,"Cpfcan_dxy_[n_Cpfcand_]/f");

    addBranch(tree,"Cpfcan_dxyerrinv",&Cpfcan_dxyerrinv_,"Cpfcan_dxyerrinv_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_dxysig",&Cpfcan_dxysig_,"Cpfcan_dxysig_[n_Cpfcand_]/f");

    addBranch(tree,"Cpfcan_dz",&Cpfcan_dz_,"Cpfcan_dz_[n_Cpfcand_]/f");

    addBranch(tree,"Cpfcan_VTX_ass",&Cpfcan_VTX_ass_,"Cpfcan_VTX_ass_[n_Cpfcand_]/f");

    addBranch(tree,"Cpfcan_fromPV",&Cpfcan_fromPV_,"Cpfcan_fromPV_[n_Cpfcand_]/f");

    addBranch(tree,"Cpfcan_drminsv",&Cpfcan_drminsv_,"Cpfcan_drminsv_[n_Cpfcand_]/f");

    //commented ones don't work
    /**///addBranch(tree,"Cpfcan_vertexChi2",&Cpfcan_vertexChi2_,"Cpfcan_vertexChi2_[n_Cpfcand_]/f");
    /**///addBranch(tree,"Cpfcan_vertexNdof",&Cpfcan_vertexNdof_,"Cpfcan_vertexNdof_[n_Cpfcand_]/f");
    /**///addBranch(tree,"Cpfcan_vertexNormalizedChi2",&Cpfcan_vertexNormalizedChi2_,"Cpfcan_vertexNormalizedChi2_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_vertex_rho",&Cpfcan_vertex_rho_,"Cpfcan_vertex_rho_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_vertex_phirel",&Cpfcan_vertex_phirel_,"Cpfcan_vertex_phirel_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_vertex_etarel",&Cpfcan_vertex_etarel_,"Cpfcan_vertex_etarel_[n_Cpfcand_]/f");
    /**///addBranch(tree,"Cpfcan_vertexRef_mass",&Cpfcan_vertexRef_mass_,"Cpfcan_vertexRef_mass_[n_Cpfcand_]/f");

    addBranch(tree,"Cpfcan_dptdpt",&Cpfcan_dptdpt_,"Cpfcan_dptdpt_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_detadeta",&Cpfcan_detadeta_,"Cpfcan_detadeta_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_dphidphi",&Cpfcan_dphidphi_,"Cpfcan_dphidphi_[n_Cpfcand_]/f");


    addBranch(tree,"Cpfcan_dxydxy",&Cpfcan_dxydxy_,"Cpfcan_dxydxy_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_dzdz",&Cpfcan_dzdz_,"Cpfcan_dzdz_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_dxydz",&Cpfcan_dxydz_,"Cpfcan_dxydz_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_dphidxy",&Cpfcan_dphidxy_,"Cpfcan_dphidxy_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_dlambdadz",&Cpfcan_dlambdadz_,"Cpfcan_dlambdadz_[n_Cpfcand_]/f");




    addBranch(tree,"Cpfcan_BtagPf_trackMomentum",&Cpfcan_BtagPf_trackMomentum_,"Cpfcan_BtagPf_trackMomentum_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackEta",&Cpfcan_BtagPf_trackEta_,"Cpfcan_BtagPf_trackEta_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackEtaRel",&Cpfcan_BtagPf_trackEtaRel_,"Cpfcan_BtagPf_trackEtaRel_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackPtRel",&Cpfcan_BtagPf_trackPtRel_,"Cpfcan_BtagPf_trackPtRel_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackPPar",&Cpfcan_BtagPf_trackPPar_,"Cpfcan_BtagPf_trackPPar_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackDeltaR",&Cpfcan_BtagPf_trackDeltaR_,"Cpfcan_BtagPf_trackDeltaR_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackPtRatio",&Cpfcan_BtagPf_trackPtRatio_,"Cpfcan_BtagPf_trackPtRatio_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackPParRatio",&Cpfcan_BtagPf_trackPParRatio_,"Cpfcan_BtagPf_trackPParRatio[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackSip3dVal",&Cpfcan_BtagPf_trackSip3dVal_,"Cpfcan_BtagPf_trackSip3dVal_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackSip3dSig",&Cpfcan_BtagPf_trackSip3dSig_,"Cpfcan_BtagPf_trackSip3dSig_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackSip2dVal",&Cpfcan_BtagPf_trackSip2dVal_,"Cpfcan_BtagPf_trackSip2dVal_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackSip2dSig",&Cpfcan_BtagPf_trackSip2dSig_,"Cpfcan_BtagPf_trackSip2dSig_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackDecayLen",&Cpfcan_BtagPf_trackDecayLen_,"Cpfcan_BtagPf_trackDecayLen_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackJetDistVal",&Cpfcan_BtagPf_trackJetDistVal_,"Cpfcan_BtagPf_trackJetDistVal_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_BtagPf_trackJetDistSig",&Cpfcan_BtagPf_trackJetDistSig_,"Cpfcan_BtagPf_trackJetDistSig_[n_Cpfcand_]/f");




    addBranch(tree,"Cpfcan_isMu",&Cpfcan_isMu_,"Cpfcan_isMu_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_isEl",&Cpfcan_isEl_,"Cpfcan_isEl_[n_Cpfcand_]/f");

    //in16 conversion broken
    // addBranch(tree,"Cpfcan_lostInnerHits",&Cpfcan_lostInnerHits_,"Cpfcan_lostInnerHits_[n_Cpfcand_]/i");
    addBranch(tree,"Cpfcan_chi2",&Cpfcan_chi2_,"Cpfcan_chi2_[n_Cpfcand_]/f");
    addBranch(tree,"Cpfcan_quality",&Cpfcan_quality_,"Cpfcan_quality_[n_Cpfcand_]/f");

    // did not give integers !!
    //  addBranch(tree,"Cpfcan_charge",&Cpfcan_charge_,"Cpfcan_charge_[n_Cpfcand_]/i");

    //Neutral Pf candidates
    addBranch(tree,"n_Npfcand", &n_Npfcand_,"n_Npfcand_/i");
    addBranch(tree,"nNpfcand", &nNpfcand_,"nNpfcand/f");

    addBranch(tree,"Npfcan_pt", &Npfcan_pt_,"Npfcan_pt_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_eta", &Npfcan_eta_,"Npfcan_eta_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_phi", &Npfcan_phi_,"Npfcan_phi_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_ptrel", &Npfcan_ptrel_,"Npfcan_ptrel_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_erel", &Npfcan_erel_,"Npfcan_erel_[n_Npfcand_]/f");

    addBranch(tree,"Npfcan_puppiw", &Npfcan_puppiw_,"Npfcan_puppiw_[n_Npfcand_]/f");


    addBranch(tree,"Npfcan_phirel",&Npfcan_phirel_,"Npfcan_phirel_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_etarel",&Npfcan_etarel_,"Npfcan_etarel_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_deltaR",&Npfcan_deltaR_,"Npfcan_deltaR_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_isGamma",&Npfcan_isGamma_,"Npfcan_isGamma_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_HadFrac",&Npfcan_HadFrac_,"Npfcan_HadFrac_[n_Npfcand_]/f");
    addBranch(tree,"Npfcan_drminsv",&Npfcan_drminsv_,"Npfcan_drminsv_[n_Npfcand_]/f");


}

void ntuple_pfCands::readEvent(const edm::Event& iEvent){

    iEvent.getByToken(svToken_, secVertices);

}



//use either of these functions

bool ntuple_pfCands::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

    reco::VertexCompositePtrCandidateCollection cpvtx=*secVertices;

    float etasign = 1.;
    if (jet.eta()<0) etasign =-1.;
    math::XYZVector jetDir = jet.momentum().Unit();
    GlobalVector jetRefTrackDir(jet.px(),jet.py(),jet.pz());
    const reco::Vertex & pv = vertices()->at(0);

    std::vector<sorting::sortingClass<pat::PackedCandidate> > sortedcharged, sortedneutrals, sortedall;

    TrackInfoBuilder trackinfo(builder);

		//create a buffer of collinear candidates to be used
		std::vector<pat::PackedCandidate*> collinears;
		std::vector<unsigned int> coll_idxs;
		collinears.reserve(2*ncollinear_);
		coll_idxs.reserve(ncollinear_);
		//draw the indexes to be split
		for(unsigned int i = 0; i <  ncollinear_; i++){		
			coll_idxs.push_back(gRandom->Integer(jet.chargedMultiplicity()));
		}
		std::sort(coll_idxs.begin(), coll_idxs.end());

    //create collection first, to be able to do some sorting
		unsigned int icharged=0;
		unsigned int icollinear=0;
    for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++){
        const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
        if(PackedCandidate){

            trackinfo.buildTrackInfo(PackedCandidate,jetDir,jetRefTrackDir,pv);

						if(PackedCandidate->charge() != 0 && icollinear != coll_idxs.size()) {
							//check if this charged must be split
							if(icharged == coll_idxs[icollinear]) {
								icollinear++; //increment collinear index

								//make two copies of the candidate and keep them in a tmp buffer
								pat::PackedCandidate* coll1 = PackedCandidate->clone();
								pat::PackedCandidate* coll2 = PackedCandidate->clone();
								collinears.push_back(coll1); 
								collinears.push_back(coll2);

								float split = gRandom->Rndm();

								//set candidate 1 momentum
								coll1->setP4(PackedCandidate->p4() * split);
								sortedall.push_back(
									sorting::sortingClass<pat::PackedCandidate>(
										coll1, trackinfo.getTrackSip2dSig(),
										-mindrsvpfcand(cpvtx,coll1), coll1->pt()/jet.pt()
										)
									);

								
								//set candidate 2 momentum and charge
								coll2->setP4(PackedCandidate->p4() * (1-split));
								coll2->setPdgId(111); //charge is imposed by the pdgid... use pi0
								sortedall.push_back(
									sorting::sortingClass<pat::PackedCandidate>(
										coll2, 0.,//neutrals do not have tracks
										-mindrsvpfcand(cpvtx,coll2), coll2->pt()/jet.pt()
										)
									);

								icharged++; //increment charged counter
								continue; //nothing to do here on
							}

							icharged++; //increment charged counter
						}

            sortedall.push_back(sorting::sortingClass<pat::PackedCandidate>
																(PackedCandidate, trackinfo.getTrackSip2dSig(),
																 -mindrsvpfcand(cpvtx,PackedCandidate), PackedCandidate->pt()/jet.pt()));

        }
    }

    std::sort(sortedall.begin(),sortedall.end(),sorting::sortingClass<pat::PackedCandidate>::compareByABCInv);

    if(false){//DEBUG
        for (const auto& s:sortedall){
            std::cout << s.sortValA <<' '<<s.sortValB<<' '<<s.sortValC<<' '<<s.get()->charge()<< std::endl;
        }
        std::cout << std::endl;
    }
    //sort by dxy significance - many infs and nans - check  - is this the reason for worse performance? FIXME


    // counts neutral and charged candicates
    n_Cpfcand_ = 0;
    n_Npfcand_ = 0;



    for (const auto& s:sortedall){

        const auto& PackedCandidate_=s.get();
        if(!PackedCandidate_)continue;

        // get the dr with the closest sv
        float drminpfcandsv_ = mindrsvpfcand(cpvtx,PackedCandidate_);

        /// This might include more than PF candidates, e.g. Reco muons and could
        /// be double counting. Needs to be checked.!!!!
        ///
        /// Split to charged and neutral candidates
        if(PackedCandidate_->charge()!=0 && max_pfcand_>n_Cpfcand_){



            Cpfcan_pt_[n_Cpfcand_] = PackedCandidate_->pt();
            Cpfcan_eta_[n_Cpfcand_] = PackedCandidate_->eta();
            Cpfcan_phi_[n_Cpfcand_] = PackedCandidate_->phi();
            Cpfcan_ptrel_[n_Cpfcand_] = catchInfsAndBound(PackedCandidate_->pt()/jet.pt()-1,0,-1,0);
            Cpfcan_erel_[n_Cpfcand_] = catchInfsAndBound(PackedCandidate_->energy()/jet.energy()-1,0,-1,0);
            Cpfcan_phirel_[n_Cpfcand_] = catchInfsAndBound(fabs(reco::deltaPhi(PackedCandidate_->phi(),jet.phi()))-0.5,0,-2,0);
            Cpfcan_etarel_[n_Cpfcand_] = catchInfsAndBound(fabs(PackedCandidate_->eta()-jet.eta())-0.5,0,-2,0);
            Cpfcan_deltaR_[n_Cpfcand_] =catchInfsAndBound(reco::deltaR(*PackedCandidate_,jet)-0.6,0,-0.6,0);
            Cpfcan_dxy_[n_Cpfcand_] = catchInfsAndBound(fabs(PackedCandidate_->dxy()),0,-50,50);


            Cpfcan_dxyerrinv_[n_Cpfcand_]=catchInfsAndBound(1/PackedCandidate_->dxyError(),0,-1, 10000.);

            Cpfcan_dxysig_[n_Cpfcand_]=catchInfsAndBound(fabs(PackedCandidate_->dxy()/PackedCandidate_->dxyError()),0.,-2000,2000);


            Cpfcan_dz_[n_Cpfcand_] = PackedCandidate_->dz();
            Cpfcan_VTX_ass_[n_Cpfcand_] = PackedCandidate_->pvAssociationQuality();

            Cpfcan_fromPV_[n_Cpfcand_] = PackedCandidate_->fromPV();

            float tempdontopt=PackedCandidate_->vx();
            tempdontopt++;

            Cpfcan_vertexChi2_[n_Cpfcand_]=PackedCandidate_->vertexChi2();
            Cpfcan_vertexNdof_[n_Cpfcand_]=PackedCandidate_->vertexNdof();
            //divided
            Cpfcan_vertexNormalizedChi2_[n_Cpfcand_]=PackedCandidate_->vertexNormalizedChi2();
            Cpfcan_vertex_rho_[n_Cpfcand_]=catchInfsAndBound(PackedCandidate_->vertex().rho(),0,-1,50);
            Cpfcan_vertex_phirel_[n_Cpfcand_]=reco::deltaPhi(PackedCandidate_->vertex().phi(),jet.phi());
            Cpfcan_vertex_etarel_[n_Cpfcand_]=etasign*(PackedCandidate_->vertex().eta()-jet.eta());
            Cpfcan_vertexRef_mass_[n_Cpfcand_]=PackedCandidate_->vertexRef()->p4().M();


            Cpfcan_puppiw_[n_Cpfcand_] = PackedCandidate_->puppiWeight();

            const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();
            reco::Track::CovarianceMatrix myCov = PseudoTrack.covariance ();
            //https://github.com/cms-sw/cmssw/blob/CMSSW_9_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h#L394

            Cpfcan_dptdpt_[n_Cpfcand_] =    catchInfsAndBound(myCov[0][0],0,-1,1);
            Cpfcan_detadeta_[n_Cpfcand_]=   catchInfsAndBound(myCov[1][1],0,-1,0.01);
            Cpfcan_dphidphi_[n_Cpfcand_]=   catchInfsAndBound(myCov[2][2],0,-1,0.1);

            /*
             * what makes the most sense here if a track is used in the fit... cerntainly no btag
             * for now leave it a t zero
             * infs and nans are set to poor quality
             */
            Cpfcan_dxydxy_[n_Cpfcand_] =    catchInfsAndBound(myCov[3][3],7.,-1,7); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dzdz_[n_Cpfcand_] =      catchInfsAndBound(myCov[4][4],6.5,-1,6.5); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dxydz_[n_Cpfcand_] =     catchInfsAndBound(myCov[3][4],6.,-6,6); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dphidxy_[n_Cpfcand_] =   catchInfs(myCov[2][3],-0.03); //zero if pvAssociationQuality ==7 ?
            Cpfcan_dlambdadz_[n_Cpfcand_]=  catchInfs(myCov[1][4],-0.03); //zero if pvAssociationQuality ==7 ?


            trackinfo.buildTrackInfo(PackedCandidate_,jetDir,jetRefTrackDir,pv);

            Cpfcan_BtagPf_trackMomentum_[n_Cpfcand_]   =catchInfsAndBound(trackinfo.getTrackMomentum(),0,0 ,1000);
            Cpfcan_BtagPf_trackEta_[n_Cpfcand_]        =catchInfsAndBound(trackinfo.getTrackEta()   ,  0,-5,5);
            Cpfcan_BtagPf_trackEtaRel_[n_Cpfcand_]     =catchInfsAndBound(trackinfo.getTrackEtaRel(),  0,-5,15);
            Cpfcan_BtagPf_trackPtRel_[n_Cpfcand_]      =catchInfsAndBound(trackinfo.getTrackPtRel(),   0,-1,4);
            Cpfcan_BtagPf_trackPPar_[n_Cpfcand_]       =catchInfsAndBound(trackinfo.getTrackPPar(),    0,-1e5,1e5 );
            Cpfcan_BtagPf_trackDeltaR_[n_Cpfcand_]     =catchInfsAndBound(trackinfo.getTrackDeltaR(),  0,-5,5 );
            Cpfcan_BtagPf_trackPtRatio_[n_Cpfcand_]    =catchInfsAndBound(trackinfo.getTrackPtRatio(), 0,-1,10 );
            Cpfcan_BtagPf_trackPParRatio_[n_Cpfcand_]  =catchInfsAndBound(trackinfo.getTrackPParRatio(),0,-10,100);
            Cpfcan_BtagPf_trackSip3dVal_[n_Cpfcand_]   =catchInfsAndBound(trackinfo.getTrackSip3dVal(), 0, -1,1e5 );
            Cpfcan_BtagPf_trackSip3dSig_[n_Cpfcand_]   =catchInfsAndBound(trackinfo.getTrackSip3dSig(), 0, -1,4e4 );
            Cpfcan_BtagPf_trackSip2dVal_[n_Cpfcand_]   =catchInfsAndBound(trackinfo.getTrackSip2dVal(), 0, -1,70 );
            Cpfcan_BtagPf_trackSip2dSig_[n_Cpfcand_]   =catchInfsAndBound(trackinfo.getTrackSip2dSig(), 0, -1,4e4 );
            Cpfcan_BtagPf_trackDecayLen_[n_Cpfcand_]   =0;
            Cpfcan_BtagPf_trackJetDistVal_[n_Cpfcand_] =catchInfsAndBound(trackinfo.getTrackJetDistVal(),0,-20,1 );
            Cpfcan_BtagPf_trackJetDistSig_[n_Cpfcand_] =catchInfsAndBound(trackinfo.getTrackJetDistSig(),0,-1,1e5 );

            // TO DO: we can do better than that by including reco::muon informations
            Cpfcan_isMu_[n_Cpfcand_] = 0;
            if(abs(PackedCandidate_->pdgId())==13) {
                Cpfcan_isMu_[n_Cpfcand_] = 1;
            }
            // TO DO: we can do better than that by including reco::electron informations
            Cpfcan_isEl_[n_Cpfcand_] = 0;
            if(abs(PackedCandidate_->pdgId())==11) {
                Cpfcan_isEl_[n_Cpfcand_] = 1;

            }

            Cpfcan_charge_[n_Cpfcand_] = PackedCandidate_->charge();
            Cpfcan_lostInnerHits_[n_Cpfcand_] = catchInfs(PackedCandidate_->lostInnerHits(),2);
            Cpfcan_chi2_[n_Cpfcand_] = catchInfsAndBound(PseudoTrack.normalizedChi2(),300,-1,300);
            //for some reason this returns the quality enum not a mask.
            Cpfcan_quality_[n_Cpfcand_] = PseudoTrack.qualityMask();

            Cpfcan_drminsv_[n_Cpfcand_] = catchInfsAndBound(drminpfcandsv_-5,0,-5,0);

            n_Cpfcand_++;
        }
        else if(max_pfcand_>n_Npfcand_){// neutral candidates
            Npfcan_pt_[n_Npfcand_] = PackedCandidate_->pt();
            Npfcan_eta_[n_Npfcand_] = PackedCandidate_->eta();
            Npfcan_phi_[n_Npfcand_] = PackedCandidate_->phi();
            Npfcan_ptrel_[n_Npfcand_] = catchInfsAndBound(PackedCandidate_->pt()/jet.pt()-1,0,-1,0);
            Npfcan_erel_[n_Npfcand_] = catchInfsAndBound(PackedCandidate_->energy()/jet.energy()-1,0,-1,0);
            Npfcan_puppiw_[n_Npfcand_] = PackedCandidate_->puppiWeight();
            Npfcan_phirel_[n_Npfcand_] = catchInfsAndBound(fabs(reco::deltaPhi(PackedCandidate_->phi(),jet.phi()))-0.5,0,-2,0);
            Npfcan_etarel_[n_Npfcand_] = catchInfsAndBound(fabs(PackedCandidate_->eta()-jet.eta())-0.5,0,-2,0);
            Npfcan_deltaR_[n_Npfcand_] = catchInfsAndBound(reco::deltaR(*PackedCandidate_,jet)-0.6,0,-0.6,0);
            Npfcan_isGamma_[n_Npfcand_] = 0;
            if(fabs(PackedCandidate_->pdgId())==22)  Npfcan_isGamma_[n_Npfcand_] = 1;
            Npfcan_HadFrac_[n_Npfcand_] = PackedCandidate_->hcalFraction();

            Npfcan_drminsv_[n_Npfcand_] = catchInfsAndBound(drminpfcandsv_-5,0,-5,0);

            n_Npfcand_++;
        }

    } // end loop over jet.numberOfDaughters()

    nCpfcand_=n_Cpfcand_;
    nNpfcand_=n_Npfcand_;

		//clear memory
		for(auto& ptr : collinears) {
			delete ptr;
		} 
		collinears.clear();

    return true; //for making cuts
}


float ntuple_pfCands::mindrsvpfcand(const std::vector<reco::VertexCompositePtrCandidate> svs, const pat::PackedCandidate* pfcand) {

    float mindr_ = 999.;
    for (unsigned int i0=0; i0<svs.size(); ++i0) {

        float tempdr_ = reco::deltaR(svs[i0],*pfcand);
        if (tempdr_<mindr_) { mindr_ = tempdr_; }

    }

    return mindr_;
}
