// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//ROOT includes
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include <TRandom.h>

//CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#if defined( __GXX_EXPERIMENTAL_CXX0X__)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#endif

struct MagneticField;


class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit DeepNtuplizer(const edm::ParameterSet&);
  ~DeepNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const;
  Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const ;
  float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const ;


  // ----------member data --------------------------- 
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
  edm::EDGetTokenT<pat::JetCollection>     jetToken_;
  const double                    jetPtMin_;
  const double                    jetPtMax_;
  const double                    jetAbsEtaMin_;
  const double                    jetAbsEtaMax_;
  const double                    gluonReduction_;
  TRandom TRandom_;



  TFile *file_ = new TFile("output.root","recreate");
  TTree *tree_ = new TTree("tree","tree");

 // labels (MC truth)
  // regressions pt, Deta, Dphi
  float gen_pt_;
  //classification
  int isB_;
  int isC_;
  int isUDS_;
  int isG_;

  // global variables 
  unsigned int npv_;

  // jet variables
  float jet_pt_;
  float  jet_eta_;

  // variables with own  coordinates (eta phi)

  // CPF charged candidate variables
  unsigned int n_Cpfcand_;
  float  Cpfcan_pt_[100];
  float  Cpfcan_phirel_[100];
  float  Cpfcan_etarel_[100];
  float  Cpfcan_puppiw_[100];
  int   Cpfcan_VTX_ass_[100];

  // covariance
  float  Cpfcan_dz_[100];
  float  Cpfcan_dxy_[100];
  float  Cpfcan_dptdpt_[100];
  float  Cpfcan_detadeta_[100];
  float  Cpfcan_dphidphi_[100];
  float  Cpfcan_dxydxy_[100];
  float  Cpfcan_dzdz_[100];
  float  Cpfcan_dxydz_[100];
  float  Cpfcan_dphidxy_[100];
  float  Cpfcan_dlambdadz_[100];
 
  // ID, skipped "charged hadron" as that is true if now the other
  int Cpfcan_isMu_[100]; // pitty that the quality is missing
  int Cpfcan_isEl_[100]; // pitty that the quality is missing
  int Cpfcan_charge_[100];

  // track quality
  int Cpfcan_lostInnerHits_[100]; 
  float Cpfcan_chi2_[100]; 
  int Cpfcan_highPurity_[100]; 

  //Neutral Pf candidates
  int n_Npfcand_;
  float  Npfcan_pt_[100]; 
  float  Npfcan_phirel_[100]; 
  float  Npfcan_etarel_[100]; 
  int  Npfcan_isGamma_[100]; 
  float  Npfcan_HadFrac_[100]; 
  
  // SV candidates
  int   sv_num_;
  float sv_pt_[100];
  float sv_eta_[100];
  float sv_phi_[100];
  float sv_mass_[100];
  int   sv_ntracks_[100];
  float sv_chi2_[100];
  float sv_ndf_[100];
  float sv_dxy_[100];
  float sv_dxyerr_[100];
  float sv_d3d_[100];
  float sv_d3derr_[100];
  float sv_costhetasvpv_[100];
  
};


DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  svToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetPtMin_(iConfig.getParameter<double>("jetPtMin")),
  jetPtMax_(iConfig.getParameter<double>("jetPtMax")),
  jetAbsEtaMin_(iConfig.getParameter<double>("jetAbsEtaMin")),
  jetAbsEtaMax_(iConfig.getParameter<double>("jetAbsEtaMax")),
  gluonReduction_(iConfig.getParameter<double>("gluonReduction"))

{
  //now do what ever initialization is needed
  usesResource("TFileService");
  // truthe labels
  tree_->Branch("gen_pt"    ,&gen_pt_    ,"gen_pt_/f"    );
  tree_->Branch("isB",&isB_, "isB_/i");
  tree_->Branch("isC",&isC_, "isC_/i");
  tree_->Branch("isUDS",&isUDS_, "isUDS_/i");
  tree_->Branch("isG",&isG_, "isG_/i");

  tree_->Branch("npv"    ,&npv_    ,"npv/i"    );

  // jet variables
  tree_->Branch("jet_pt", &jet_pt_);
  tree_->Branch("jet_eta", &jet_eta_);

  // Cpfcanditates per jet
  tree_->Branch("n_Cpfcand", &n_Cpfcand_,"n_Cpfcand_/i");
  tree_->Branch("Cpfcan_pt", &Cpfcan_pt_,"Cpfcan_pt_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_phirel",&Cpfcan_phirel_,"Cpfcan_phirel_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_etarel",&Cpfcan_etarel_,"Cpfcan_etarel_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_puppiw",&Cpfcan_puppiw_,"Cpfcan_puppiw_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_dxy",&Cpfcan_dxy_,"Cpfcan_dxy_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_dz",&Cpfcan_dz_,"Cpfcan_dz_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_VTX_ass",&Cpfcan_VTX_ass_,"Cpfcan_VTX_ass_[n_Cpfcand_]/i");
  tree_->Branch("Cpfcan_dptdpt",&Cpfcan_dptdpt_,"Cpfcan_dptdpt_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_detadeta",&Cpfcan_detadeta_,"Cpfcan_detadeta_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_dphidphi",&Cpfcan_dphidphi_,"Cpfcan_dphidphi_[n_Cpfcand_]/f");

  // FIXME gave INFs?
  //  tree_->Branch("Cpfcan_dxydxy",&Cpfcan_dxydxy_,"Cpfcan_dxydxy_[n_Cpfcand_]/f");
  //  tree_->Branch("Cpfcan_dzdz",&Cpfcan_dzdz_,"Cpfcan_dzdz_[n_Cpfcand_]/f");
  // tree_->Branch("Cpfcan_dxydz",&Cpfcan_dxydz_,"Cpfcan_dxydz_[n_Cpfcand_]/f");
  // tree_->Branch("Cpfcan_dphidxy",&Cpfcan_dphidxy_,"Cpfcan_dphidxy_[n_Cpfcand_]/f");
  // tree_->Branch("Cpfcan_dlambdadz",&Cpfcan_dlambdadz_,"Cpfcan_dlambdadz_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_isMu",&Cpfcan_isMu_,"Cpfcan_isMu_[n_Cpfcand_]/i");
  tree_->Branch("Cpfcan_isEl",&Cpfcan_isEl_,"Cpfcan_isEl_[n_Cpfcand_]/i");
  // tree_->Branch("Cpfcan_lostInnerHits",&Cpfcan_lostInnerHits_,"Cpfcan_lostInnerHits_[n_Cpfcand_]/i");
  tree_->Branch("Cpfcan_chi2",&Cpfcan_chi2_,"Cpfcan_chi2_[n_Cpfcand_]/f");
  tree_->Branch("Cpfcan_highPurity",&Cpfcan_highPurity_,"Cpfcan_highPurity_[n_Cpfcand_]/i");

  // did not give integers !!
  //  tree_->Branch("Cpfcan_charge",&Cpfcan_charge_,"Cpfcan_charge_[n_Cpfcand_]/i");

  //Neutral Pf candidates 
  tree_->Branch("n_Npfcand", &n_Npfcand_,"n_Npfcand_/i");
  tree_->Branch("Npfcan_pt", &Npfcan_pt_,"Npfcan_pt_[n_Npfcand_]/f");
  tree_->Branch("Npfcan_phirel",&Npfcan_phirel_,"Npfcan_phirel_[n_Npfcand_]/f");
  tree_->Branch("Npfcan_etarel",&Npfcan_etarel_,"Npfcan_etarel_[n_Npfcand_]/f");
  tree_->Branch("Npfcan_isGamma",&Npfcan_isGamma_,"Npfcan_isGamma_[n_Npfcand_]/i");
  tree_->Branch("Npfcan_HadFrac",&Npfcan_HadFrac_,"Npfcan_HadFrac_[n_Npfcand_]/f");
  
  // SV candidates
  tree_->Branch("sv_num"         ,&sv_num_         ,"sv_num_/i"                  );
  tree_->Branch("sv_pt"          ,&sv_pt_          ,"sv_pt_[sv_num_]/f"          );
  tree_->Branch("sv_eta"         ,&sv_eta_         ,"sv_eta_[sv_num_]/f"         );
  tree_->Branch("sv_phi"         ,&sv_phi_         ,"sv_phi_[sv_num_]/f"         );
  tree_->Branch("sv_mass"        ,&sv_mass_        ,"sv_mass_[sv_num_]/f"        );
  tree_->Branch("sv_ntracks"     ,&sv_ntracks_     ,"sv_ntracks_[sv_num_]/i"     );
  tree_->Branch("sv_chi2"        ,&sv_chi2_        ,"sv_chi2_[sv_num_]/f"        );
  tree_->Branch("sv_ndf"         ,&sv_ndf_         ,"sv_ndf_[sv_num_]/f"         );
  tree_->Branch("sv_dxy"         ,&sv_dxy_         ,"sv_dxy_[sv_num_]/f"         );
  tree_->Branch("sv_dxyerr"      ,&sv_dxyerr_      ,"sv_dxyerr_[sv_num_]/f"      );
  tree_->Branch("sv_d3d"         ,&sv_d3d_         ,"sv_d3d_[sv_num_]/f"         );
  tree_->Branch("sv_d3derr"      ,&sv_d3derr_      ,"sv_d3err_[sv_num_]/f"       );
  tree_->Branch("sv_costhetasvpv",&sv_costhetasvpv_,"sv_costhetasvpv_[sv_num_]/f");
  
}


DeepNtuplizer::~DeepNtuplizer()
{
  file_->Close();
}


// ------------ method called for each event  ------------
void
DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
  iEvent.getByToken(svToken_, secVertices);


  npv_ = vertices->size();


  // loop over the jets
  for (const pat::Jet &jet : *jets) {

    if (npv_==0) { continue; } // we need a vertex

   // some cuts to contrin training region 
    if ( jet.pt() < jetPtMin_ ||  jet.pt() > jetPtMax_ ) continue;                  // apply jet pT cut
    if ( jet.eta() < fabs(jetAbsEtaMin_) ||jet.eta() > fabs(jetAbsEtaMax_) ) continue; // apply jet eta cut
    // often we have way to many gluons that we do not need. This randomply reduces the gluons
    if (gluonReduction_>0)
      {
	if(jet.partonFlavour()==21)
	  {
	    if(TRandom_.Uniform()>gluonReduction_) continue; 
	  }
      }

    // truth labels
    gen_pt_ = 0.;
    if(jet.genJet()!=NULL)   gen_pt_ =  jet.genJet()->pt();
    isB_= int(abs(jet.partonFlavour())==5);
    isC_= int(abs(jet.partonFlavour())==4);
    isUDS_= int( (abs(jet.partonFlavour())>0) && (abs(jet.partonFlavour())<4 ));
    isG_= int(jet.partonFlavour()==21);

    jet_pt_ = jet.pt();
    jet_eta_ = jet.eta();
      float etasign = 1.;
      if (jet.eta()<0) etasign =-1.;

      // counts neutral and charged candicates
      n_Cpfcand_ = 0;
      n_Npfcand_ = 0;
   
      for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++)
	{

	  /// This might include more than PF candidates, e.g. Reco muons and could
	  /// be double counting. Needs to be checked.!!!!
	  ///
	  /// Split to charged and neutral candidates

	  const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
	  if(PackedCandidate_->charge()!=0)
	    {
	      Cpfcan_pt_[n_Cpfcand_] = PackedCandidate_->pt();
	      Cpfcan_phirel_[n_Cpfcand_] = reco::deltaPhi(PackedCandidate_->phi(),jet.phi());
	      Cpfcan_etarel_[n_Cpfcand_] = etasign*(PackedCandidate_->eta()-jet.eta());
	      Cpfcan_dxy_[n_Cpfcand_] = PackedCandidate_->dxy();
	      Cpfcan_dz_[n_Cpfcand_] = PackedCandidate_->dz();
	      Cpfcan_VTX_ass_[n_Cpfcand_] = PackedCandidate_->pvAssociationQuality();
	      Cpfcan_puppiw_[n_Cpfcand_] = PackedCandidate_->puppiWeight();
	      
	      const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();
	      reco::Track::CovarianceMatrix myCov = PseudoTrack.covariance ();
	      //https://github.com/cms-sw/cmssw/blob/CMSSW_9_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h#L394
	      Cpfcan_dptdpt_[n_Cpfcand_] = myCov[0][0];
	      Cpfcan_detadeta_[n_Cpfcand_]= myCov[1][1];
	      Cpfcan_dphidphi_[n_Cpfcand_]= myCov[2][2];
	      Cpfcan_dxydxy_[n_Cpfcand_] =  myCov[3][3]; //zero if pvAssociationQuality ==7 ?
	      Cpfcan_dzdz_[n_Cpfcand_] =  myCov[4][4]; //zero if pvAssociationQuality ==7 ?
	      Cpfcan_dxydz_[n_Cpfcand_] =  myCov[3][4]; //zero if pvAssociationQuality ==7 ?
	      Cpfcan_dphidxy_[n_Cpfcand_] =  myCov[2][3]; //zero if pvAssociationQuality ==7 ?
	      Cpfcan_dlambdadz_[n_Cpfcand_] =  myCov[1][4]; //zero if pvAssociationQuality ==7 ?

	      // TO DO: we can do better than that by including reco::muon informations
	      Cpfcan_isMu_[n_Cpfcand_] = 0; 
	      if(abs(PackedCandidate_->pdgId())==13)    Cpfcan_isMu_[n_Cpfcand_] = 1;

	      // TO DO: we can do better than that by including reco::electron informations
	      Cpfcan_isEl_[n_Cpfcand_] = 0;
	      if(abs(PackedCandidate_->pdgId())==11)    Cpfcan_isEl_[n_Cpfcand_] = 1;

	      Cpfcan_charge_[n_Cpfcand_] = PackedCandidate_->charge();
	      Cpfcan_lostInnerHits_[n_Cpfcand_] = PackedCandidate_->lostInnerHits();
	      Cpfcan_chi2_[n_Cpfcand_] = PseudoTrack.normalizedChi2();
	      Cpfcan_highPurity_[n_Cpfcand_] = PseudoTrack.highPurity;
	      n_Cpfcand_++;
	}
	  else{// neutral candidates
	    Npfcan_pt_[n_Npfcand_] = PackedCandidate_->pt();
	    Npfcan_phirel_[n_Npfcand_] = reco::deltaPhi(PackedCandidate_->phi(),jet.phi());
	    Npfcan_etarel_[n_Npfcand_] = etasign*(PackedCandidate_->eta()-jet.eta());
	    Npfcan_isGamma_[n_Npfcand_] = 0;
	    if(fabs(PackedCandidate_->pdgId())==22)  Npfcan_isGamma_[n_Npfcand_] = 1;
	      Npfcan_HadFrac_[n_Npfcand_] = PackedCandidate_->hcalFraction();
	      n_Npfcand_++;
	  }
	  
	} // end loop over jet.numberOfDaughters()


      // get the SV that are in cone DR<0.4 of the jet
      const reco::Vertex & pv = (*vertices)[0];
      sv_num_ = 0; 
      for (const reco::VertexCompositePtrCandidate &sv : *secVertices) {
	
	if (reco::deltaR(sv,jet)>0.4) { continue; }

	sv_pt_[sv_num_]           = sv.pt();
	sv_eta_[sv_num_]          = sv.eta();
	sv_phi_[sv_num_]          = sv.phi();
	sv_mass_[sv_num_]         = sv.mass();
	sv_ntracks_[sv_num_]      = sv.numberOfDaughters();
	sv_chi2_[sv_num_]         = sv.vertexChi2();
	sv_ndf_[sv_num_]          = sv.vertexNdof();
	sv_dxy_[sv_num_]          = vertexDxy(sv,pv).value();
	sv_dxyerr_[sv_num_]       = vertexDxy(sv,pv).error();
	sv_d3d_[sv_num_]          = vertexD3d(sv,pv).value();
	sv_d3derr_[sv_num_]       = vertexD3d(sv,pv).error();
	sv_costhetasvpv_[sv_num_] = vertexDdotP(sv,pv); // the pointing angle (i.e. the angle between the sum of the momentum
	                                                // of the tracks in the SV and the flight direction betwen PV and SV)
	sv_num_++;

      } // end of looping over the secondary vertices





    tree_->Fill();
  }
}


Measurement1D DeepNtuplizer::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const {
  VertexDistanceXY dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

Measurement1D DeepNtuplizer::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) const {
  VertexDistance3D dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

float DeepNtuplizer::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) const {
  reco::Candidate::Vector p = sv.momentum();
  reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
  return p.Unit().Dot(d.Unit());
}


// ------------ method called once each job just before starting event loop  ------------
void
DeepNtuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DeepNtuplizer::endJob()
{
  file_->cd();
  tree_->Write();
  file_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
