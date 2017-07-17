//getListOfBranches

#include "../interface/ntuple_JetInfo.h"
#include "../interface/ntuple_SV.h"
#include "../interface/ntuple_bTagVars.h"
#include "../interface/ntuple_pfCands.h"
#include "../interface/ntuple_FatJetInfo.h"
#include "TFile.h"
#include <vector>
#include "TH1F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"

int main(int argc, char *argv[]) {

    std::vector<ntuple_content*> branchinfos;
    branchinfos.push_back(new ntuple_JetInfo());
    branchinfos.push_back(new ntuple_SV());
    branchinfos.push_back(new ntuple_bTagVars());
    branchinfos.push_back(new ntuple_pfCands());
    branchinfos.push_back(new ntuple_FatJetInfo());

    if (argc < 3)
        return -1;
    TString samplefile = argv[1];
    TString outdir = argv[2];

    TFile* f = new TFile(samplefile, "READ");
    if (!f)
        return -2;
    TTree* t = (TTree*) f->Get("deepntuplizer/tree");

    if (!t)
        return -3;

    std::vector<TString> allbranches;
    for (const auto& b : branchinfos) {
        std::vector<TString> tb = b->getListOfBranches();
        allbranches.insert(allbranches.end(), tb.begin(), tb.end());
    }

    TCanvas cv;
    for (const auto& b : allbranches) {
        t->SetLineColor(kBlack);
        t->Draw(b + ">>" + b + "ALL", "", "normalized");
        TH1F *histoB = (TH1F*) gROOT->FindObject(b + "ALL");
        histoB->Draw("hist");

        t->SetLineColor(kRed);
        t->Draw(b + ">>" + b + "B", "isB+isBB+isLeptonicB+isLeptonicB_C",
                "same,normalized");
        TH1F* histo = (TH1F*) gROOT->FindObject(b + "B");
        histo->Draw("hist,same");
        t->SetLineColor(kGreen);
        t->Draw(b + ">>" + b + "C", "isC", "same,normalized");
        histo = (TH1F*) gROOT->FindObject(b + "C");
        histo->Draw("hist,same");
        t->SetLineColor(kBlue);
        t->Draw(b + ">>" + b + "L", "isUD+isS", "same,normalized");
        histo = (TH1F*) gROOT->FindObject(b + "L");
        histo->Draw("hist,same");
        t->SetLineColor(kMagenta);
        t->Draw(b + ">>" + b + "G", "isG", "same,normalized");
        histo = (TH1F*) gROOT->FindObject(b + "G");
        histo->Draw("hist,same");

        cv.BuildLegend(0.12,0.7,0.3,0.87);

        cv.Print(outdir + "/" + b + ".pdf");
    }

    return 0;

}
