//getListOfBranches

#include "../interface/ntuple_JetInfo.h"
#include "../interface/ntuple_SV.h"
#include "../interface/ntuple_bTagVars.h"
#include "../interface/ntuple_pfCands.h"
#include "../interface/ntuple_FatJetInfo.h"
#include "../interface/ntuple_DeepVertex.h"
#include "TFile.h"
#include <vector>
#include "TH1F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"


int main(int argc, char *argv[]) {

    std::vector<ntuple_content*> branchinfos;
    branchinfos.push_back(new ntuple_JetInfo());
    branchinfos.push_back(new ntuple_SV());
    branchinfos.push_back(new ntuple_bTagVars());
    branchinfos.push_back(new ntuple_pfCands());
  //  branchinfos.push_back(new ntuple_FatJetInfo());
  //  branchinfos.push_back(new ntuple_DeepVertex());

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



    for (const auto& b : allbranches) {
        std::vector<TH1F *> hists;

        TFile * outtfile= new TFile(outdir + "/" + b + ".root","RECREATE");

        TCanvas cv;
        t->SetLineColor(kBlack);
        t->Draw(b + ">>" + b + "ALL", "", "normalized");
        hists.push_back((TH1F*) gROOT->FindObject(b + "ALL"));


        t->SetLineColor(kRed);
        t->Draw(b + ">>" + b + "B", "isB+isBB+isGBB+isLeptonicB+isLeptonicB_C",
                "same,normalized");
        hists.push_back((TH1F*) gROOT->FindObject(b + "B"));


        t->SetLineColor(kGreen);
        t->Draw(b + ">>" + b + "C", "isC+isCC+isGCC", "same,normalized");
        hists.push_back( (TH1F*) gROOT->FindObject(b + "C"));

        t->SetLineColor(kBlue);
        t->Draw(b + ">>" + b + "L", "isUD+isS", "same,normalized");
        hists.push_back((TH1F*) gROOT->FindObject(b + "L"));

        t->SetLineColor(kMagenta);
        t->Draw(b + ">>" + b + "G", "isG", "same,normalized");
        hists.push_back( (TH1F*) gROOT->FindObject(b + "G"));

        cv.Clear();
        gStyle->SetOptStat(0);
        cv.SetLeftMargin(0.15);
        cv.SetBottomMargin(0.15);

        float max=-1e10;
        float min=1e10;
        for(size_t i=0;i<hists.size();i++){
            if(i==0)
                hists.at(i)->Draw("hist");
            else
                hists.at(i)->Draw("hist,same");
            float thismax=hists.at(i)->GetMaximum();
            float thismin=hists.at(i)->GetMinimum();
            if(thismin<min)min=thismin;
            if(thismax>max)max=thismax;
        }
        max*=1.2;

        hists.at(0)->GetYaxis()->SetRangeUser(min,max);

        cv.BuildLegend(0.72,0.6,0.9,0.87);


        cv.Write();
        cv.Print(outdir + "/" + b + ".pdf");

        for(auto h:hists)
            h->Write();

        outtfile->Close();

    }

    return 0;

}
