/*
 * ntuple_DeepVertex.h
 *
 *  Created on: 23 June 2017
 *      Author: Seth Moortgat
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_DEEPVERTEX_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_DEEPVERTEX_H_

#include "ntuple_content.h"
#include "trackVars2.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class ntuple_DeepVertex: public ntuple_content{
public:

    ntuple_DeepVertex(double jetR = 0.4);
    ~ntuple_DeepVertex();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);


    void setCandidatesToken(const edm::EDGetTokenT<edm::View<pat::PackedCandidate> > & t){
        CandidateToken=t;
     }

private:

    // seed candidates
    static constexpr size_t max_seeds=10;
    
    int n_seeds=0;
    
    double seed_pt[max_seeds];
    double seed_eta[max_seeds];
    double seed_phi[max_seeds];
    double seed_mass[max_seeds];

    double seed_dz[max_seeds];
    double seed_dxy[max_seeds];
    double seed_3D_ip[max_seeds];
    double seed_3D_sip[max_seeds];
    double seed_2D_ip[max_seeds];
    double seed_2D_sip[max_seeds];
    double seed_3D_signedIp[max_seeds];
    double seed_3D_signedSip[max_seeds];
    double seed_2D_signedIp[max_seeds];
    double seed_2D_signedSip[max_seeds];
    //int seed_JetMatch[max_seeds];
 
    double seed_chi2reduced[max_seeds];
    double seed_nPixelHits[max_seeds];
    double seed_nHits[max_seeds];
    double seed_jetAxisDistance[max_seeds];
    double seed_jetAxisDlength[max_seeds];
    
    int seed_nNearTracks[max_seeds];
    
    
    //nearest track candidates
    static constexpr size_t max_nearestTrk=200; // 20 per seed
    
    int n_NearTracksTotal=0;
    
    double nearTracks_pt[max_nearestTrk];
    double nearTracks_eta[max_nearestTrk];
    double nearTracks_phi[max_nearestTrk];
    double nearTracks_mass[max_nearestTrk];
    double nearTracks_dz[max_nearestTrk];
    double nearTracks_dxy[max_nearestTrk];
    double nearTracks_3D_ip[max_nearestTrk];
    double nearTracks_3D_sip[max_nearestTrk];
    double nearTracks_2D_ip[max_nearestTrk];
    double nearTracks_2D_sip[max_nearestTrk];
    double nearTracks_PCAdist[max_nearestTrk];
    double nearTracks_PCAdsig[max_nearestTrk];      
    double nearTracks_PCAonSeed_x[max_nearestTrk];
    double nearTracks_PCAonSeed_y[max_nearestTrk];
    double nearTracks_PCAonSeed_z[max_nearestTrk];      
    double nearTracks_PCAonSeed_xerr[max_nearestTrk];
    double nearTracks_PCAonSeed_yerr[max_nearestTrk];
    double nearTracks_PCAonSeed_zerr[max_nearestTrk];      
    double nearTracks_PCAonTrack_x[max_nearestTrk];
    double nearTracks_PCAonTrack_y[max_nearestTrk];
    double nearTracks_PCAonTrack_z[max_nearestTrk];      
    double nearTracks_PCAonTrack_xerr[max_nearestTrk];
    double nearTracks_PCAonTrack_yerr[max_nearestTrk];
    double nearTracks_PCAonTrack_zerr[max_nearestTrk]; 
    double nearTracks_dotprodTrack[max_nearestTrk];
    double nearTracks_dotprodSeed[max_nearestTrk];
    double nearTracks_dotprodTrackSeed2D[max_nearestTrk];
    double nearTracks_dotprodTrackSeed3D[max_nearestTrk];
    double nearTracks_dotprodTrackSeedVectors2D[max_nearestTrk];
    double nearTracks_dotprodTrackSeedVectors3D[max_nearestTrk];      
    double nearTracks_PCAonSeed_pvd[max_nearestTrk];
    double nearTracks_PCAonTrack_pvd[max_nearestTrk];
    double nearTracks_PCAjetAxis_dist[max_nearestTrk];
    double nearTracks_PCAjetMomenta_dotprod[max_nearestTrk];
    double nearTracks_PCAjetDirs_DEta[max_nearestTrk];
    double nearTracks_PCAjetDirs_DPhi[max_nearestTrk];
    
    
    // IVF cut parameters (HARDCODED?? OR CONFIGURABLE IN PYTHON CONFIG)
    float min3DIPValue=0.005;
    float min3DIPSignificance=1.2;
    int max3DIPValue=9999.;
    int max3DIPSignificance=9999.;
    

    //tokens to be defined from main analyzer
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> > CandidateToken;

    //helper:
    edm::Handle<edm::View<pat::PackedCandidate> > tracks;
    
    // builder    
    edm::ESHandle<TransientTrackBuilder> builder;
    
    // temporary containers
    trackVars2 myTrack;
    std::vector<trackVars2> nearTracks; 
    std::multimap<double,std::pair<const reco::TransientTrack*,const std::vector<trackVars2> > > SortedSeedsMap;


};



#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_DEEPVERTEX_H_ */
