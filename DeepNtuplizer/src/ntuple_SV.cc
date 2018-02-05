/*
 * ntuple_SV.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../interface/ntuple_SV.h"
// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"


const reco::Vertex * ntuple_SV::spvp_;

ntuple_SV::ntuple_SV(std::string prefix, double jetR):ntuple_content(jetR),sv_num_(0){
    prefix_ = prefix;
}
ntuple_SV::~ntuple_SV(){}


void ntuple_SV::getInput(const edm::ParameterSet& iConfig){

}

void ntuple_SV::initBranches(TTree* tree){
    // SV candidates
    addBranch(tree,(prefix_+"n_sv").c_str()         ,&sv_num_         ,(prefix_+"sv_num_/i").c_str()     );
    addBranch(tree,(prefix_+"nsv").c_str()          ,&nsv_          ,(prefix_+"nsv_/f").c_str()         );
    addBranch(tree,(prefix_+"sv_pt").c_str()          ,&sv_pt_          ,(prefix_+"sv_pt_["+prefix_+"sv_num_]/f").c_str()        );
    addBranch(tree,(prefix_+"sv_eta").c_str()          ,&sv_eta_          ,(prefix_+"sv_eta_["+prefix_+"sv_num_]/f").c_str()        );
    addBranch(tree,(prefix_+"sv_phi").c_str()          ,&sv_phi_          ,(prefix_+"sv_phi_["+prefix_+"sv_num_]/f").c_str()        );
    addBranch(tree,(prefix_+"sv_etarel").c_str()         ,&sv_etarel_         ,(prefix_+"sv_etarel_["+prefix_+"sv_num_]/f").c_str()         );
    addBranch(tree,(prefix_+"sv_phirel").c_str()         ,&sv_phirel_         ,(prefix_+"sv_phirel_["+prefix_+"sv_num_]/f").c_str()         );
    addBranch(tree,(prefix_+"sv_deltaR").c_str()         ,&sv_deltaR_         ,(prefix_+"sv_deltaR_["+prefix_+"sv_num_]/f").c_str()         );
    addBranch(tree,(prefix_+"sv_mass").c_str()        ,&sv_mass_        ,(prefix_+"sv_mass_["+prefix_+"sv_num_]/f").c_str()        );
    addBranch(tree,(prefix_+"sv_ntracks").c_str()     ,&sv_ntracks_     ,(prefix_+"sv_ntracks_["+prefix_+"sv_num_]/f").c_str()     );
    addBranch(tree,(prefix_+"sv_chi2").c_str()        ,&sv_chi2_        ,(prefix_+"sv_chi2_["+prefix_+"sv_num_]/f").c_str()        );
    addBranch(tree,(prefix_+"sv_ndf").c_str()         ,&sv_ndf_         ,(prefix_+"sv_ndf_["+prefix_+"sv_num_]/f").c_str()         );
    addBranch(tree,(prefix_+"sv_normchi2").c_str()    ,&sv_normchi2_   ,(prefix_+"sv_normchi2_["+prefix_+"sv_num_]/f").c_str()     );
    addBranch(tree,(prefix_+"sv_dxy").c_str()         ,&sv_dxy_         ,(prefix_+"sv_dxy_["+prefix_+"sv_num_]/f").c_str()         );
    addBranch(tree,(prefix_+"sv_dxyerr").c_str()      ,&sv_dxyerr_      ,(prefix_+"sv_dxyerr_["+prefix_+"sv_num_]/f").c_str()      );
    addBranch(tree,(prefix_+"sv_dxysig").c_str()      ,&sv_dxysig_      ,(prefix_+"sv_dxysig_["+prefix_+"sv_num_]/f").c_str()      );
    addBranch(tree,(prefix_+"sv_d3d").c_str()         ,&sv_d3d_         ,(prefix_+"sv_d3d_["+prefix_+"sv_num_]/f").c_str()         );
    addBranch(tree,(prefix_+"sv_d3derr").c_str()      ,&sv_d3derr_      ,(prefix_+"sv_d3err_["+prefix_+"sv_num_]/f").c_str()       );
    addBranch(tree,(prefix_+"sv_d3dsig").c_str()      ,&sv_d3dsig_      ,(prefix_+"sv_d3dsig_["+prefix_+"sv_num_]/f").c_str()      );
    addBranch(tree,(prefix_+"sv_costhetasvpv").c_str(),&sv_costhetasvpv_,(prefix_+"sv_costhetasvpv_["+prefix_+"sv_num_]/f").c_str());
    addBranch(tree,(prefix_+"sv_enratio").c_str()     ,&sv_enratio_     ,(prefix_+"sv_enratio_["+prefix_+"sv_num_]/f").c_str());


}


void ntuple_SV::readEvent(const edm::Event& iEvent){


}


bool ntuple_SV::compareDxyDxyErr(const reco::VertexCompositePtrCandidate &sva,const reco::VertexCompositePtrCandidate &svb){
    reco::Vertex pv=*spvp_;
    float adxy= ntuple_SV::vertexDxy(sva,pv).value();
    float bdxy= ntuple_SV::vertexDxy(svb,pv).value();
    float aerr=ntuple_SV::vertexDxy(sva,pv).error();
    float berr=ntuple_SV::vertexDxy(svb,pv).error();

    float asig=ntuple_SV::catchInfs(adxy/aerr,0.);
    float bsig=ntuple_SV::catchInfs(bdxy/berr,0.);
    return bsig<asig;
}

bool ntuple_SV::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::Event& iEvent, const  edm::View<pat::Jet> * coll){


    const float jet_uncorr_e=jet.correctedJet("Uncorrected").energy();

    const reco::Vertex & pv =    vertices()->at(0);

    sv_num_ = 0;

    reco::VertexCompositePtrCandidateCollection cpvtx=*secVertices();

    spvp_ =   & vertices()->at(0);
    std::sort(cpvtx.begin(),cpvtx.end(),ntuple_SV::compareDxyDxyErr);

    float etasign=1;
    etasign++; //avoid unused warning
    if(jet.eta()<0)etasign=-1;

    double jet_radius = jetR();
    if (jet_radius<0){
      // subjets: use maxDR(subjet, pfcand)
      for (unsigned idau=0; idau<jet.numberOfDaughters(); ++idau){
        double dR = reco::deltaR(*jet.daughter(idau), jet);
        if (dR>jet_radius)
          jet_radius = dR;
      }
    }

    for (const reco::VertexCompositePtrCandidate &sv : cpvtx) {

        if (reco::deltaR(sv,jet)>jet_radius) { continue; }
        if((int)max_sv>sv_num_){

            sv_pt_[sv_num_]           = sv.pt();
            sv_eta_[sv_num_]          = sv.eta();
            sv_phi_[sv_num_]          = sv.phi();
            sv_etarel_[sv_num_]       = catchInfsAndBound(fabs(sv.eta()-jet.eta())-0.5,0,-2,0);
            sv_phirel_[sv_num_]       = catchInfsAndBound(fabs(reco::deltaPhi(sv.phi(),jet.phi()))-0.5,0,-2,0);
            sv_deltaR_[sv_num_]       = catchInfsAndBound(fabs(reco::deltaR(sv,jet))-0.5,0,-2,0);
            sv_mass_[sv_num_]         = sv.mass();
            sv_ntracks_[sv_num_]      = sv.numberOfDaughters();
            sv_chi2_[sv_num_]         = sv.vertexChi2();
            sv_ndf_[sv_num_]          = sv.vertexNdof();
            sv_normchi2_[sv_num_]     = catchInfsAndBound(sv_chi2_[sv_num_]/sv_ndf_[sv_num_],1000,-1000,1000);
            sv_dxy_[sv_num_]          = vertexDxy(sv,pv).value();
            sv_dxyerr_[sv_num_]       = catchInfsAndBound(vertexDxy(sv,pv).error()-2,0,-2,0);
            sv_dxysig_[sv_num_]       = catchInfsAndBound(sv_dxy_[sv_num_]/vertexDxy(sv,pv).error() ,0,-1,800);
            sv_d3d_[sv_num_]          = vertexD3d(sv,pv).value();
            sv_d3derr_[sv_num_]       = catchInfsAndBound(vertexD3d(sv,pv).error()-2,0,-2,0);
            sv_d3dsig_[sv_num_]       = catchInfsAndBound(vertexD3d(sv,pv).value()/vertexD3d(sv,pv).error(),0,-1,800);
            sv_costhetasvpv_[sv_num_] = vertexDdotP(sv,pv); // the pointing angle (i.e. the angle between the sum of the momentum
            // of the tracks in the SV and the flight direction betwen PV and SV)

            sv_enratio_[sv_num_]=sv.energy()/jet_uncorr_e;



            sv_num_++;
        }
    } // end of looping over the secondary vertices
    nsv_=sv_num_;

    return true;
}









///helpers seldomly touched



Measurement1D ntuple_SV::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

Measurement1D ntuple_SV::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

float ntuple_SV::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}

