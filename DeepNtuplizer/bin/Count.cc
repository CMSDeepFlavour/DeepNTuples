#include "TTree.h"
#include <TFile.h>
#include <TChain.h>
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

Int_t fCurrent; //!current Tree number in a TChain

Long64_t LoadTree(Long64_t entry, TChain * fChain)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  return centry;
}


int main(){

  TChain *fChain30to50 = new TChain("btagana/ttree");
  fChain30to50->Add("/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/incl_qcd_30/180117_082536/0000/J*root");
  TChain *fChain50to80 = new TChain("btagana/ttree");
  fChain50to80->Add("/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/incl_qcd_50/180117_082611/0000/J*root");
  TChain *fChain80to120 = new TChain("btagana/ttree");
  fChain80to120->Add("/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/incl_qcd_120/180117_082718/0000/J*root");
  TChain *fChain120to170 = new TChain("btagana/ttree");
  fChain120to170->Add("/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/incl_qcd_120/180117_082718/0000/J*root");
  TChain *fChain170to300 = new TChain("btagana/ttree");
  fChain170to300->Add("/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/incl_qcd_170/180117_082753/0000/J*root");
  TChain *fChain300to470 = new TChain("btagana/ttree");
  fChain300to470->Add("/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/incl_qcd_300/180117_082827/0000/J*root");
  TChain *fChain470to600 = new TChain("btagana/ttree");
  fChain470to600->Add("/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/incl_qcd_470/180117_082905/0000/J*root");


  Float_t         pthat;
  TBranch        *b_pthat;
  fCurrent = -1;
  Long64_t           n15_20, n20_30, n30_50,n50_80,n80_120,n120_170,n170_300,n300_470,n470_600,n600_800, n800_1000, n1000_inf;
  Long64_t           n15_30, n120_150, n150_inf;


  n15_20   =0;
  n20_30   =0;
  n30_50   = fChain30to50->GetEntries();
  cout << "30 to 50 done" << endl;
  n50_80   = fChain50to80->GetEntries();
  cout << "50 to 80 done" << endl;
  n80_120  = fChain80to120->GetEntries();
  cout << "80 to 120 done" << endl;
  n120_170 = fChain120to170->GetEntries();
  cout << "120 to 170 done" << endl;
  n170_300 = fChain170to300->GetEntries();
  cout << "170 to 300 done" << endl;
  n300_470 = fChain300to470->GetEntries();
  cout << "300 to 470 done" << endl;
  n470_600 = fChain470to600->GetEntries();
  n600_800 =0;
  n800_1000=0;
  n1000_inf=0;
  n15_30   =0;
  n120_150 =0;
  n150_inf =0;

  int   Nevent = 0;
  
  // Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  
  int qcdtype = 0;
  int sqrtstev = 13;
  cout << " To write in runCode.C " << endl;
  if (qcdtype==0) { // inclusive qcd 
    cout << " double   n15    = 0. ; " <<endl;
    cout << " double   n20    = "<< n15_30 << "; " << endl;
  }
  else { // MuEnriched qcd
    cout << " double   n15    = "<< n15_20 << "; " << endl;
    cout << " double   n20    = "<< n20_30 << "; " << endl;
  }
  cout << " double   n30    = "<< n30_50 << "; " << endl;
  cout << " double   n50    = "<< n50_80 << "; " << endl;
  cout << " double   n80    = "<< n80_120 << "; " << endl;
  if (sqrtstev!=7) { // 8 TeV
    cout << " double   n120  = "<< n120_170 << "; " << endl;
    cout << " double   n170  = "<< n170_300 << "; " << endl;
    cout << " double   n300  = "<< n300_470 << "; " << endl;
    cout << " double   n470  = "<< n470_600 << "; " << endl;
    cout << " double   n600  = "<< n600_800 << "; " << endl;
    cout << " double   n800  = "<< n800_1000 << "; " << endl;
    cout << " double   n1000 = "<< n1000_inf << "; " << endl;
  }
  else if (qcdtype==1) { // 7 TeV MuEnriched qcd
    cout << " double   n120  = "<< n120_150 << "; " << endl;
    cout << " double   n170  = "<< n150_inf << "; " << endl;
    cout << " double   n300  = 0.; " << endl;
    cout << " double   n470  = 0.; " << endl; 
    cout << " double   n600  = 0.; " << endl;
    cout << " double   n800  = 0.; " << endl;
    cout << " double   n1000 = 0.; " << endl;
  }
  cout << " m.Fill_nevent(n15,n20,n30,n50,n80,n120,n170,n300,n470,n600,n800,n1000);" << endl;

  return 1;
}
