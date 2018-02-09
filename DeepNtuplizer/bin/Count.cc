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

  TChain *fChain = new TChain("btagana/ttree");
  fChain->Add("/afs/cern.ch/work/e/ebols/public/TestSample/MC*root");
  Float_t         pthat;
  TBranch        *b_pthat;
  fCurrent = -1;
  fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
  double           n15_20, n20_30, n30_50,n50_80,n80_120,n120_170,n170_300,n300_470,n470_600,n600_800, n800_1000, n1000_inf;
  double           n15_30, n120_150, n150_inf;
  bool            use15_20,use20_30, use30_50,use50_80,use80_120,use120_170,use170_300,use300_470,use470_600,use600_800, use800_1000, use1000_inf;
  bool            use15_30, use150_inf; 


  n15_20   =0;
  n20_30   =0;
  n30_50   =0;
  n50_80   =0;
  n80_120  =0;
  n120_170 =0;
  n170_300 =0;
  n300_470 =0;
  n470_600 =0;
  n600_800 =0;
  n800_1000=0;
  n1000_inf=0;
  n15_30   =0;
  n120_150 =0;
  n150_inf =0;

  int   Nevent = 0;
  if (fChain == 0) return 0;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry,fChain);
    if (ientry < 0) break;
    
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    Nevent++;
    
    if(Nevent%50000 ==0 && Nevent!=0) cout << " number of processed events is " << Nevent << endl;    
    if      ( pthat >=  15. && pthat <  20. ){
      use15_20=true;
      n15_20++;
    } 
    else if ( pthat >=  20. && pthat <  30. ){
      use20_30=true;    
      n20_30++;
    } 
    else if ( pthat >=  30. && pthat <  50. ){
      use30_50=true;    
      n30_50++;
    } 
    else if ( pthat >=  50. && pthat <  80. ){
      use50_80=true;    
      n50_80++;
    }
    else if ( pthat >=  80. && pthat < 120. ){
      use80_120=true;    
      n80_120++;
    }
    else if ( pthat >= 120. && pthat < 170. ){
      use120_170=true;    
      n120_170++;
    }
    else if ( pthat >= 170. && pthat < 300. ){
      use170_300=true;    
      n170_300++;
    }
    else if ( pthat >= 300. && pthat < 470. ){
      use300_470=true;    
      n300_470++;
    }
    else if ( pthat >= 470. && pthat < 600. ){
      use470_600=true;    
      n470_600++;
    }
    else if ( pthat >= 600. && pthat <= 800. ){
      use600_800=true;    
      n600_800++;
    }
    else if ( pthat >= 800. && pthat <= 1000. ){
      use800_1000=true;    
      n800_1000++;
    }
    else if ( pthat >= 1000. ){
      use1000_inf=true;    
      n1000_inf++;
    }

    if      ( pthat >=  15. && pthat <  30. ){
      use15_30=true;
      n15_30++;
    } 
    else if ( pthat >= 120. && pthat < 150. ){
      n120_150++;
    }
    else if ( pthat >= 150. ){
      use150_inf=true;    
      n150_inf++;
    }

  }  
  
  if  (use15_20)  cout << "Run with QCD15_20   sample with " << n15_20  << " events" <<endl; 
  if  (use20_30)  cout << "Run with QCD20_30   sample with " << n20_30  << " events" <<endl; 
  if  (use30_50)  cout << "Run with QCD30_50   sample with " << n30_50  << " events" <<endl;
  if  (use50_80)  cout << "Run with QCD50_80   sample with " << n50_80  << " events" <<endl;
  if  (use80_120) cout << "Run with QCD80_120  sample with " << n80_120 << " events" <<endl;
  if  (use120_170)cout << "Run with QCD120_170 sample with " << n120_170<< " events" <<endl;
  if  (use170_300)cout << "Run with QCD170_300 sample with " << n170_300<< " events" <<endl;
  if  (use300_470)cout << "Run with QCD300_470 sample with " << n300_470<< " events" <<endl;
  if  (use470_600)cout << "Run with QCD470_600 sample with " << n470_600<< " events" <<endl;
  if  (use600_800)cout << "Run with QCD600_800 sample with " << n600_800<< " events" <<endl;
  if  (use800_1000)cout << "Run with QCD800_1000 sample with " << n800_1000<< " events" <<endl;
  if  (use1000_inf)cout << "Run with QCD1000 sample with " << n1000_inf<< " events" <<endl;
  cout << endl;
  if (use15_30) cout << "Run with QCD15_30   sample with " << n15_30  << " events" <<endl;
  if (use150_inf) cout << "Run with QCD150   sample with " << n150_inf  << " events" <<endl;

  cout << endl;
  cout << endl;

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
