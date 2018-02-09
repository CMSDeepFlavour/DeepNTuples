#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include <TRandom.h>
#include <TH2F.h>
#include <TH1F.h>
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>

using namespace std;

TString gentype;
bool qcdtype;
int sqrtstev;
float sum_xs; 
float x_section[12]; 
float nmc_evt_vect[12];

float Jet_pt_[1000]; 
float Jet_eta_[1000]; 
int nJet_;

void Fill_nevent(double n15,double n20,double n30,double n50,double n80,double n120,double n170,double n300,double n470,double n600, double n800, double n1000){

 
  nmc_evt_vect[0]=n15;
  nmc_evt_vect[1]=n20;
  nmc_evt_vect[2]=n30;  
  nmc_evt_vect[3]=n50;  
  nmc_evt_vect[4]=n80;  
  nmc_evt_vect[5]=n120;
  nmc_evt_vect[6]=n170;  
  nmc_evt_vect[7]=n300;  
  nmc_evt_vect[8]=n470;   
  nmc_evt_vect[9]=n600;
  nmc_evt_vect[10]=n800;
  nmc_evt_vect[11]=n1000;
  
}


bool passTrigger(TString trigger, int pttrig, int BitTrigger[], int nJet, float Jet_pt[], float Jet_eta[]) {


  // FOR 2012 Trigger ! Not valid for 2011...

  bool passTrig=false;
  //bool Jet30  = false, Jet60  = false, Jet150 = false, Jet190 = false, Jet240 = false;
  bool Jet40  = false, Jet60=false,  Jet80  = false, Jet140 = false;
  bool Jet200 = false, Jet260 = false, Jet320 = false;
  bool Jet400 = false, Jet450 = false, Jet500 = false;
  bool Jet20  = false, Jet70  = false, Jet110 = false, Jet170 = false, Jet300 = false; 
  int triggerIdx = 0, bitIdx = 0;

  if ( trigger=="jet") {
    triggerIdx = 0;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;

    triggerIdx = 1;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet60  = true;
   
    triggerIdx = 2;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet80  = true;
   
    triggerIdx = 3;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet140 = true;
   
    triggerIdx = 4;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet200 = true;
   
    triggerIdx = 5;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet260 = true;
   
    triggerIdx = 6;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet320 = true;

    triggerIdx = 7;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet400 = true;

    triggerIdx = 8;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet450 = true;

    triggerIdx = 9;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet500 = true;


    if ( pttrig ==  40 && Jet40 )  passTrig=true;
    if ( pttrig ==  60 && Jet60 )  passTrig=true;
    if ( pttrig ==  80 && Jet80 )  passTrig=true;
    if ( pttrig == 140 && Jet140 ) passTrig=true;
    if ( pttrig == 200 && Jet200 ) passTrig=true;
    if ( pttrig == 260 && Jet260 ) passTrig=true;
    if ( pttrig == 320 && Jet320 ) passTrig=true;
    if ( pttrig == 400 && Jet400 ) passTrig=true;
    if ( pttrig == 450 && Jet450 ) passTrig=true;
    if ( pttrig == 500 && Jet500 ) passTrig=true;

    if (!passTrig) { return false; }

    //-----------------------------------
    //Determine if there is at least 
    //one jet which pass the trigger
    //in the event => away from the TO
    //-----------------------------------

    bool JetPtCut = false;
    for (int ijet=0; ijet<nJet ; ijet++) {
      float ptjet = Jet_pt[ijet];
      float etajet = fabs(Jet_eta[ijet]);
      //            if (      pttrig ==  40 && ptjet >  60. && etajet < 2.4 ) JetPtCut = true;
      if (      pttrig ==  40 && ptjet >  50. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig ==  60 && ptjet >  70. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig ==  80 && ptjet > 100. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig == 140 && ptjet > 160. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig == 200 && ptjet > 220. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig == 260 && ptjet > 300. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig == 320 && ptjet > 360. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig == 400 && ptjet > 425. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig == 450 && ptjet > 475. && etajet < 2.4 ) JetPtCut = true;
      else if ( pttrig == 500 && ptjet > 525. && etajet < 2.4 ) JetPtCut = true;
    }
    if (passTrig && JetPtCut) {return true;}
    else {return false;}
  }
   
  else if ( trigger=="btag" ) {
    triggerIdx = 32;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet20  = true;
   
    triggerIdx = 33;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;
   
    triggerIdx = 34;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet70  = true;
   
    triggerIdx = 35;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet110 = true;
   
    triggerIdx = 36;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet170 = true;

    triggerIdx = 37;
    bitIdx = int(triggerIdx/32);
    if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet300 = true;

    if ( pttrig ==  20 && Jet20 )  passTrig=true;
    if ( pttrig ==  40 && Jet40 )  passTrig=true;
    if ( pttrig ==  70 && Jet70 )  passTrig=true;
    if ( pttrig == 110 && Jet110 ) passTrig=true;
    if ( pttrig == 170 && Jet170 ) passTrig=true;
    if ( pttrig == 300 && Jet300 ) passTrig=true;
    
    
    if (!passTrig) { return false; }
    
    //-----------------------------------
    //Determine if there is at least 
    //two jets which pass the trigger
    //in the event => away from the TO
    //-----------------------------------

    int njtrig=0;
    if (pttrig ==300) njtrig+=1;
    for (int ijet = 0; ijet < nJet; ijet++) {
      float ptjet = Jet_pt[ijet];
      float etajet = fabs(Jet_eta[ijet]);
      //            if ( pttrig ==20  &&  ptjet > 60. && etajet < 2.4 )  njtrig++;
      if ( pttrig == 20  &&  ptjet > 30. && etajet < 2.4 )  njtrig++;
      if ( pttrig == 40  &&  ptjet > 50. && etajet < 2.4 )  njtrig++;
      if ( pttrig == 70  &&  ptjet > 80. && etajet < 2.4 )  njtrig++;
      if ( pttrig == 110  &&  ptjet > 120. && etajet < 2.4 ) njtrig++;
      if ( pttrig == 170  &&  ptjet > 190. && etajet < 2.4 ) njtrig++;
      if ( pttrig == 300  &&  ptjet > 330. && etajet < 2.4 ) njtrig++;
    }
    if (passTrig && njtrig>1) {return true;}
    else {return false;}

  }
  return false;
}

void SetXS(TString generator, bool MuEnriched, int TeV){
  /*
      pythia Inclusive  8TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000        --> 10 fichiers
      pythia MuEnriched 8TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 12 fichies
      herwig Inclusive  8TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 11 fichiers
      herwig MuEnriched 8TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000, n1000 --> 12 fichies
      pythia Inclusive  7TeV : 0     , n15-30, n30-50, n50-80,  n80-120,  n120-170, n170-300,  n300-470,  n470-600, n600-800, n800-1000;       --> 10 fichiers
      pythia MuEnriched 7TeV : n15-20, n20-30, n30-50, n50-80,  n80-120,  n120-150, n150-plus                                                  --> 7 fichiers
  */


  // http://cms.cern.ch/iCMS/jsp/mcprod/admin/requestmanagement.jsp?pwg=BTV&campid=Summer12_DR53X
  //                              0  &     15-30
  //                          15-20        20-30         30-50         50-80      80-120       120-170      170-300      300-470     470-600    600-800   800-1000   10000         
  double pythia_xs8  [12]={         0.,  9.8828742E8,  66285328.0,   8148778.0,  1033680.0,    156293.3,   34138.15,    1759.549,   113.8791,      26.9921, 3.550036,     0.};  // PTHAT 0, 15-30
  double pythia_xs8MU[12]={2.73858e+06,   1.8655e+06,   806298.0,    176187.6,      40448,     7463.94,    2299.752,    151.8048,    11.79648,    2.690196,   0.368781,     0.0849078};
  
  double herwig_xs8  [12]={         0.,        7.9E8,     5.32E7,   6545700.0,   833630.0,    126490.0,    27935.0,      1461.0,      95.25,     22.73,     2.997,     0.665}; // PTHAT 0, 15-30
  double herwig_xs8MU[12]={1.56419e+06,      783598.,     334139,     81821.2,   13754.9,     3162.25,    751.452,     53.4726,    3.22898,   1.00467,     0.128871,  0.025935};
  
  double pythia_xs7  [12]={         0.,   815900000.0, 53120000.0,    6359000.0,  784300.0,    115100.0,     24260.0,       1168.0,     70.22,      15.55,   1.844,    0 };
  double pythia_xs7MU[12]={  1471168.0,     1224034.0,   578463.0,    144421.74,   29048.7,   4440.2215,  2837.6712,          0,        0.,          0.,     0.,     0.};  
  //remark for pythia_xs7MU, it is for PTHAT 120-150 & 150-inf (and not 120-170 and 170-300)

  //    double pythia_xs13[12]   = { 0., 2.237E9, 1.615E8, 2.211E7, 3000114.3, 493200., 120300., 7475., 587.1, 167., 28.25, 0. };
  //    double pythia_xs13MU[12] = { 1.576E9*0.0039, 6.753E8*0.0065, 1.643E8*0.00816, 2.181E7*0.01522, 2.999E6*0.02424, 493200.*0.0473, 120300.*0.0676, 7475.*0.0864, 587.1*0.1024, 167.*0.0996, 28.25*0.1033, 8.975*0.1097 };
  //RunIISpring15
  double pythia_xs13[12]   = { 0., 1.83741E9, 1.40932E8, 1.92043E7, 2762530., 471100., 117276., 7823., 648.2, 186.9, 32.293, 9.4183 };
  double pythia_xs13MU[12] = { 1.27319E9*0.003, 5.58528E8*0.0053, 1.39803E8*0.01182, 1.92225E7*0.02276, 2.758420E6*0.03844, 469797.*0.05362, 117989.*0.07335, 7820.25*0.10196, 645.528*0.12242, 187.109*0.13412, 32.3486*0.14552, 10.4305*0.15544 };

  // cout << generator << " " << MuEnriched << " " << TeV << endl;
 
  if      (generator=="pythia" && !MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs8[i]; }        }
  else if (generator=="pythia" &&  MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs8MU[i]; }        }

  else if (generator=="herwig" && !MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=herwig_xs8[i]; }        }
  else if (generator=="herwig" &&  MuEnriched && TeV==8){       for (int i=0; i<12; i++){ x_section[i]=herwig_xs8MU[i]; }        }

  else if (generator=="pythia" && !MuEnriched && TeV==7){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs7[i]; }        }
  else if (generator=="pythia" &&  MuEnriched && TeV==7){       for (int i=0; i<12; i++){ x_section[i]=pythia_xs7MU[i]; }        }
  
  else if (generator=="pythia" && !MuEnriched && TeV==13){      for (int i=0; i<12; i++){ x_section[i]=pythia_xs13[i]; }      }
  else if (generator=="pythia" &&  MuEnriched && TeV==13){      for (int i=0; i<12; i++){ x_section[i]=pythia_xs13MU[i]; }      }

}

void SetSumXS(){
  sum_xs=0.0;
  
  for (int i=0; i<12; i++) {
    if (nmc_evt_vect[i]>0)  sum_xs+=x_section[i];
  }

}

float GetEvtWeight(float pthat){

  float WeightXS=0.0;
  float  nevt  =0.0;
  float  xs    =0.0;
  
  
  if (qcdtype==0) { // inclusive
    if ( pthat >=  15. && pthat <  30. ){
      nevt  =nmc_evt_vect[1];
      xs=x_section[1];
    } 
  }
  else {
    if ( pthat >=  15. && pthat <  20. ){
      nevt  =nmc_evt_vect[0];
      xs=x_section[0];
    } 
    else if ( pthat >=  20. && pthat <  30. ){
      nevt  =nmc_evt_vect[1];
      xs=x_section[1];
    } 
  }
     
  if ( pthat >=  30. && pthat <  50. ){
    nevt  =nmc_evt_vect[2];
    xs=x_section[2];
  }
  else if ( pthat >=  50. && pthat <  80. ){
    nevt  =nmc_evt_vect[3];
    xs=x_section[3];
  }     
  else if ( pthat >=  80. && pthat <  120. ){
    nevt  =nmc_evt_vect[4];
    xs=x_section[4];
  } 
  if (gentype=="pythia" && qcdtype==1 && sqrtstev==7) {
    if ( pthat >=  120. && pthat <  150. ){
      nevt  =nmc_evt_vect[5];
      xs=x_section[5];
    } 
    else if ( pthat >=  150. ){
      nevt  =nmc_evt_vect[6];
      xs=x_section[6];
    }   
    else {
      nevt=0;
      xs=0;
    }
  }
  else {
    if ( pthat >=  120. && pthat <  170. ){
      nevt  =nmc_evt_vect[5];
      xs=x_section[5];
    }
    else if ( pthat >=  170. && pthat <  300. ){
      nevt  =nmc_evt_vect[6];
      xs=x_section[6];
    }    
    else if ( pthat >=  300. && pthat <  470. ){
      nevt  =nmc_evt_vect[7];
      xs=x_section[7];
    } 
    else if ( pthat >=  470. && pthat <  600. ){
      nevt  =nmc_evt_vect[8];
      xs=x_section[8];
    }
    else if ( pthat >=  600. && pthat <  800. ){
      nevt  =nmc_evt_vect[9];
      xs=x_section[9];
    }
    else if ( pthat >=  800. && pthat <  1000. ){
      nevt  =nmc_evt_vect[10];
      xs=x_section[10];
    }
    else if ( pthat >=  1000. ){
      nevt  =nmc_evt_vect[11];
      xs=x_section[11];
    }
  }     
  
  if (nevt>0) { WeightXS =xs/(sum_xs*nevt);  }
  else { WeightXS=0.; }
  
  if ( pthat <1. ) WeightXS=1. ;  // data

  return  WeightXS;
   
}

int main(int argc, char *argv[]){

  if(argc != 4){
    cout << "missing input arguments, should be input file, output file, and then MC or DATA" << endl;
    return 0;
  }
  
  const char* inputfile = argv[1]; 
  const char* outputfile = argv[2];
  bool isMC = false;
  string mcarg = argv[3]; 
  if(mcarg == "MC"){
    isMC = true;
  }
  
  
  TFile *ntfile = new TFile(inputfile,"READ");
  TTree *btagana = (TTree*)ntfile->Get("btagana/ttree");
  TFile *f = new TFile(outputfile,"RECREATE");
  TH1D* nPV_data                = new TH1D("nPV_data",              "nPV_data",              180,0,180);
  TH1D* nPV_mc_unw              = new TH1D("nPV_mc_unw",            "nPV_mc_unw",                180,0,180);
  TH1D* pt_hat_data                  = new TH1D("pt_hat_data",                "pt_hat_data",                100,0,1000);
  TH1D* pt_hat_mc                  = new TH1D("pt_hat_mc",                "pt_hat_mc",                100,0,1000);
  int ecount = btagana->GetEntries();
  float pthat;
  int nBitTrigger, nPV;
  int BitTrigger[3];


  double   n15    = 0.; 
  double   n20    = 0.; 
  double   n30    = 19503604.; 
  double   n50    = 18828383.; 
  double   n80    = 25098641.; 
  double   n120  = 25098641.; 
  double   n170  = 29394249.; 
  double   n300  = 29115849.; 
  double   n470  = 25011264.; 
  double   n600  = 0; 
  double   n800  = 0; 
  double   n1000 = 0; 
  Fill_nevent(n15,n20,n30,n50,n80,n120,n170,n300,n470,n600,n800,n1000);

  gentype = "pythia";
  qcdtype = false;
  sqrtstev = 13;
  SetXS(gentype,qcdtype,sqrtstev); 
  SetSumXS();
  btagana->SetBranchStatus("*", 0);

  btagana->SetBranchStatus("Jet_pt", 1);
  btagana->SetBranchStatus("Jet_eta", 1);
  btagana->SetBranchStatus("nJet", 1);
  btagana->SetBranchStatus("nBitTrigger", 1);
  btagana->SetBranchStatus("BitTrigger", 1);
  btagana->SetBranchStatus("nPV", 1);
  btagana->SetBranchStatus("pthat", 1);

  btagana->SetBranchAddress("Jet_pt", &Jet_pt_);
  btagana->SetBranchAddress("Jet_eta", &Jet_eta_);
  btagana->SetBranchAddress("nJet", &nJet_);
  btagana->SetBranchAddress("nBitTrigger", &nBitTrigger);
  btagana->SetBranchAddress("BitTrigger", &BitTrigger);
  btagana->SetBranchAddress("nPV", &nPV);
  btagana->SetBranchAddress("pthat", &pthat);

  for (Int_t i=0;i<ecount;i++) {
    btagana->GetEntry(i);
    bool isTrigOK = passTrigger("jet", 40, BitTrigger,nJet_,Jet_pt_,Jet_eta_);
    if (!isTrigOK) continue;
    float ww_woPU = 1.0;
    if(isMC){
      ww_woPU = GetEvtWeight(pthat);
      nPV_mc_unw->Fill(nPV,ww_woPU);    
      pt_hat_mc->Fill(pthat,ww_woPU);
    }
    else{
      nPV_data->Fill(nPV,ww_woPU);    
      pt_hat_data->Fill(pthat,ww_woPU);      
    }
  }
  f->Write();
  return 1;


}
