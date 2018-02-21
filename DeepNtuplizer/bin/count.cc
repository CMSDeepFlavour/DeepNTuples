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
#include "Riostream.h"

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


int main(int argc, char *argv[]){

  if(argc != 2){
    cout << "missing input directories" << endl;
    return 0;
  }
  string dir = argv[1];
  string dir30to50,dir50to80,dir80to120,dir120to170,dir170to300,dir300to470,dir470to600;
  ifstream in;
  in.open(dir);
  std::getline(in,dir30to50);
  std::getline(in,dir50to80);
  std::getline(in,dir80to120);
  std::getline(in,dir120to170);
  std::getline(in,dir170to300);
  std::getline(in,dir300to470);
  std::getline(in,dir470to600);



  TChain *fChain30to50 = new TChain("btagana/ttree");
  fChain30to50->Add((dir30to50).c_str());
  TChain *fChain50to80 = new TChain("btagana/ttree");
  fChain50to80->Add((dir50to80).c_str());
  TChain *fChain80to120 = new TChain("btagana/ttree");
  fChain80to120->Add((dir80to120).c_str());
  TChain *fChain120to170 = new TChain("btagana/ttree");
  fChain120to170->Add((dir120to170).c_str());
  TChain *fChain170to300 = new TChain("btagana/ttree");
  fChain170to300->Add((dir170to300).c_str());
  TChain *fChain300to470 = new TChain("btagana/ttree");
  fChain300to470->Add((dir300to470).c_str());
  TChain *fChain470to600 = new TChain("btagana/ttree");
  fChain470to600->Add((dir470to600).c_str());


  fCurrent = -1;
  Long64_t           n15_20, n20_30, n30_50,n50_80,n80_120,n120_170,n170_300,n300_470,n470_600,n600_800, n800_1000, n1000_inf;
  

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
 
  
  
  cout << " To write in config file after the pvPath and MC configurations " << endl;
  
  cout << n15_20 << endl;
  cout << n20_30 << endl;
  cout << n30_50 << endl;
  cout << n50_80 << endl;
  cout << n80_120 << endl;
  cout << n120_170 << endl;
  cout << n170_300 << endl;
  cout << n300_470 << endl;
  cout << n470_600 << endl;
  cout << n600_800 << endl;
  cout << n800_1000 << endl;
  cout << n1000_inf << endl;
  
  return 1;
}
