/*
 * treeReader.h
 *
 *  Created on: 7 Feb 2018
 *      Author: ebols
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_TREEREADER_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_TREEREADER_H_

#include "ntuple_content.h"

/*
 * For global jet info such as eta, pt, gen info
 */
class treeReader{
public:
    void initBranches(TTree* tree, bool MCtest);
    void initBranches(TTree* tree, std::string On, bool MCtest);

    static constexpr size_t max_jetN_=1000;
    //jet general
    float Jet_uncorrpt[max_jetN_];
    int nJet_;
    int Jet_nFirstSV_[max_jetN_];
    int Jet_nFirstTrkTagVarCSV_[max_jetN_];
    int Jet_nLastTrkTagVarCSV_[max_jetN_];
    int Jet_nFirstTrkEtaRelTagVarCSV_[max_jetN_];
    int Jet_nLastTrkEtaRelTagVarCSV_[max_jetN_];
    float Jet_pt_[max_jetN_];
    float Jet_eta_[max_jetN_];
    float trackJetPt_[max_jetN_];              // track-based jet transverse momentum
    int jetNTracks_[max_jetN_];              // tracks associated to jet
    float jetNSecondaryVertices_[max_jetN_];   // number of secondary vertices associated to the jet
    float trackSumJetEtRatio_[max_jetN_];      // ratio of track sum transverse energy over jet energy
    float trackSumJetDeltaR_[max_jetN_];       // pseudoangular distance between jet axis and track fourvector sum
    float trackSip2dValAboveCharm_[max_jetN_]; // track 2D signed impact parameter of first track lifting mass above charm
    float trackSip2dSigAboveCharm_[max_jetN_]; // track 2D signed impact parameter significance of first track lifting mass above charm
    float trackSip3dValAboveCharm_[max_jetN_]; // track 3D signed impact parameter of first track lifting mass above charm
    float trackSip3dSigAboveCharm_[max_jetN_]; // track 3D signed impact parameter significance of first track lifting mass above charm
    float vertexCategory_[max_jetN_];          // category of secondary vertex (Reco, Pseudo, No)
    //track info
    int   n_jetNTracksEtaRel_[max_jetN_];        // tracks associated to jet for which trackEtaRel is calculated
    int   n_jetNSelectedTracks_[max_jetN_];

    float jetNTracksEtaRel_[max_jetN_];
    float jetNSelectedTracks_[max_jetN_];

    static constexpr size_t max_jetNSelectedTracks_=10000;

    float trackMomentum_[max_jetNSelectedTracks_];    // track momentum
    float trackEta_[max_jetNSelectedTracks_];         // track pseudorapidity
    float trackPhi_[max_jetNSelectedTracks_];         // track polar angle
    float trackPtRel_[max_jetNSelectedTracks_];       // track transverse momentum, relative to the jet axis
    float trackPPar_[max_jetNSelectedTracks_];        // track parallel momentum, along the jet axis
    float trackDeltaR_[max_jetNSelectedTracks_];      // track pseudoangular distance from the jet axis
    float trackPtRatio_[max_jetNSelectedTracks_];     // track transverse momentum, relative to the jet axis, normalized to its energy
    float trackPParRatio_[max_jetNSelectedTracks_];   // track parallel momentum, along the jet axis, normalized to its energy
    float trackSip2dVal_[max_jetNSelectedTracks_];    // track 2D signed impact parameter
    float trackSip2dSig_[max_jetNSelectedTracks_];    // track 2D signed impact parameter significance
    float trackSip3dVal_[max_jetNSelectedTracks_];    // track 3D signed impact parameter
    float trackSip3dSig_[max_jetNSelectedTracks_];    // track 3D signed impact parameter significance
    float trackDecayLenVal_[max_jetNSelectedTracks_]; // track decay length
    float trackDecayLenSig_[max_jetNSelectedTracks_]; // track decay length significance
    float trackJetDistVal_[max_jetNSelectedTracks_];  // minimum track approach distance to jet axis
    float trackJetDistSig_[max_jetNSelectedTracks_];  // minimum track approach distance to jet axis significance
    float trackEtaRel_[max_jetNSelectedTracks_];      // track pseudorapidity, relative to the jet axis
    //SV info
    int   n_StoredVertices_[max_jetN_];
    float NStoredVertices_[max_jetN_];

    static constexpr size_t max_nStoredVertices_=50;

    float vertexMass_[max_nStoredVertices_];          // mass of track sum at secondary vertex
    float vertexNTracks_[max_nStoredVertices_];       // number of tracks at secondary vertex
    float vertexEnergyRatio_[max_nStoredVertices_];   // ratio of energy at secondary vertex over total energy
    float vertexJetDeltaR_[max_nStoredVertices_];     // pseudoangular distance between jet axis and secondary vertex direction
    float flightDistance2dVal_[max_nStoredVertices_]; // transverse distance between primary and secondary vertex
    float flightDistance2dSig_[max_nStoredVertices_]; // transverse distance significance between primary and secondary vertex
    float flightDistance3dVal_[max_nStoredVertices_]; // distance between primary and secondary vertex
    float flightDistance3dSig_[max_nStoredVertices_]; // distance significance between primary and secondary vertex

    Int_t Jet_hadronFlavour[max_jetN_];


    Float_t         Jet_DeepFlavourBDisc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourCvsLDisc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourCvsBDisc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourB[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourBB[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourLEPB[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourC[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourUDS[max_jetN_];   //[nJet]
    Float_t         Jet_DeepFlavourG[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVBDisc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVBDiscN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVBDiscP[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVCvsLDisc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVCvsLDiscN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVCvsLDiscP[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVCvsBDisc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVCvsBDiscN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVCvsBDiscP[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVb[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVl[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVbb[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVcc[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVbN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVcN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVlN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVbbN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVccN[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVbP[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVcP[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVlP[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVbbP[max_jetN_];   //[nJet]
    Float_t         Jet_DeepCSVccP[max_jetN_];   //[nJet]

    Int_t           nBitTrigger;
    Int_t           BitTrigger[3];
    Float_t         pthat;
    Int_t           nPV;


private:

    template <class T>
    void addBranch(TTree* t, const char* name,  T* address){
      t->SetBranchStatus(name, 1);
      t->SetBranchAddress(name,address);
    }

};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_TREEREADER_H_ */
