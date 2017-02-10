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

class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit DeepNtuplizer(const edm::ParameterSet&);
  ~DeepNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data --------------------------- 
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> >     jetToken_; 
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchReclusterToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchWithNuToken_;
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
  float gen_pt_Recluster_;
  float gen_pt_WithNu_;
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

};


DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
  genJetMatchReclusterToken_(consumes<edm::Association<reco::GenJetCollection> >(iConfig.getParameter<edm::InputTag>( "genJetMatchRecluster" ))),
  genJetMatchWithNuToken_(consumes<edm::Association<reco::GenJetCollection> >(iConfig.getParameter<edm::InputTag>( "genJetMatchWithNu" ))),
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
  tree_->Branch("gen_pt_Recluster"    ,&gen_pt_Recluster_    ,"gen_pt_Recluster_/f"    );
  tree_->Branch("gen_pt_WithNu"    ,&gen_pt_WithNu_    ,"gen_pt_WithNu_/f"    );
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

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchRecluster;
  iEvent.getByToken(genJetMatchReclusterToken_, genJetMatchRecluster);

  edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchWithNu;
  iEvent.getByToken(genJetMatchWithNuToken_, genJetMatchWithNu);

  // clear vectors
  npv_ = vertices->size();


  // loop over the jets
  for (edm::View<pat::Jet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++) {
    unsigned int idx = itJet - jets->begin();       

   // some cuts to contrin training region 
    if ( itJet->pt() < jetPtMin_ ||  itJet->pt() > jetPtMax_ ) continue;                  // apply jet pT cut
    if ( itJet->eta() < fabs(jetAbsEtaMin_) ||itJet->eta() > fabs(jetAbsEtaMax_) ) continue; // apply jet eta cut
    // often we have way to many gluons that we do not need. This randomply reduces the gluons
    if (gluonReduction_>0)
      {
	if(itJet->partonFlavour()==21)
	  {
	    if(TRandom_.Uniform()>gluonReduction_) continue; 
	  }
      }

    // truth labels
    gen_pt_ = 0.;
    if(itJet->genJet()!=NULL)   gen_pt_ =  itJet->genJet()->pt();
    const edm::RefToBase<pat::Jet> patJetRef = jets->refAt(idx);      
    //    const edm::RefToBase<reco::Jet> jetRef(patJetRef.castTo<reco::JetRef>());     
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
    isB_= int(abs(itJet->partonFlavour())==5);
    isC_= int(abs(itJet->partonFlavour())==4);
    isUDS_= int( (abs(itJet->partonFlavour())>0) && (abs(itJet->partonFlavour())<4 ));
    isG_= int(itJet->partonFlavour()==21);

    jet_pt_ = itJet->pt();
    jet_eta_ = itJet->eta();
      float etasign = 1.;
      if (itJet->eta()<0) etasign =-1.;

      // counts neutral and charged candicates
      n_Cpfcand_ = 0;
      n_Npfcand_ = 0;
   
      for (unsigned int i = 0; i <  itJet->numberOfDaughters(); i++)
	{

	  /// This might include more than PF candidates, e.g. Reco muons and could
	  /// be double counting. Needs to be checked.!!!!
	  ///
	  /// Split to charged and neutral candidates

	  const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(itJet->daughter(i));
	  if(PackedCandidate_->charge()!=0)
	    {
	      Cpfcan_pt_[n_Cpfcand_] = PackedCandidate_->pt();
	      Cpfcan_phirel_[n_Cpfcand_] = reco::deltaPhi(PackedCandidate_->phi(),itJet->phi());
	      Cpfcan_etarel_[n_Cpfcand_] = etasign*(PackedCandidate_->eta()-itJet->eta());
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
	    Npfcan_phirel_[n_Npfcand_] = reco::deltaPhi(PackedCandidate_->phi(),itJet->phi());
	    Npfcan_etarel_[n_Npfcand_] = etasign*(PackedCandidate_->eta()-itJet->eta());
	    Npfcan_isGamma_[n_Npfcand_] = 0;
	    if(fabs(PackedCandidate_->pdgId())==22)  Npfcan_isGamma_[n_Npfcand_] = 1;
	      Npfcan_HadFrac_[n_Npfcand_] = PackedCandidate_->hcalFraction();
	      n_Npfcand_++;
	  }
	  
	} // end loop over itJet->numberOfDaughters()


    tree_->Fill();
  }
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
