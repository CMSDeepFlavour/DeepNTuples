


//getListOfBranches

#include "../interface/ntuple_JetInfo.h"
#include "../interface/ntuple_SV.h"
#include "../interface/ntuple_bTagVars.h"
#include "../interface/ntuple_pfCands.h"
#include "TFile.h"
#include <vector>
#include "TH1F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"

int main(int argc, char *argv[]){

	std::vector<ntuple_content*> branchinfos;
	branchinfos.push_back(new ntuple_JetInfo());
	branchinfos.push_back(new ntuple_SV());
	branchinfos.push_back(new ntuple_bTagVars());
	branchinfos.push_back(new ntuple_pfCands());

	if(argc < 3)
		return -1;
	TString samplefile=argv[1];
	TString outdir=argv[2];


	TFile* f = new TFile(samplefile,"READ");
	if(!f)
		return -2;
	TTree* t =(TTree*) f->Get("deepntuplizer/tree");

	if(!t)
		return -3;


	std::vector<TString> allbranches;
	for(const auto& b:branchinfos){
		std::vector<TString> tb=b->getListOfBranches();
		allbranches.insert(allbranches.end(),tb.begin(),tb.end());
	}

	TCanvas cv;
	for(const auto& b:allbranches){
		t->SetLineColor(kRed);
		t->Draw(b+">>"+b+"B","isB+isBB+isLeptonicB+isLeptonicB_C","normalized");
		TH1F *histo = (TH1F*)gROOT->FindObject(b+"B");
		histo->Draw("hist");
		t->SetLineColor(kGreen);
		t->Draw(b+">>"+b+"C","isC","same,normalized");
		histo = (TH1F*)gROOT->FindObject(b+"C");
		histo->Draw("hist,same");
		t->SetLineColor(kBlue);
		t->Draw(b+">>"+b+"L","isUD+isS+isG","same,normalized");
		histo = (TH1F*)gROOT->FindObject(b+"L");
		histo->Draw("hist,same");
		cv.Print(outdir+"/"+b+".pdf");
	}



	return 0;

}
