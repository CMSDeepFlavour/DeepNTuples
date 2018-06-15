######
# takes data and MC histograms of pileup distribution and computes the pileup weights
######

import ROOT
import pdb

file_pileup_MC = ROOT.TFile(
    "/afs/cern.ch/work/d/dwalter/DeepNTuples/CMSSW_8_0_29/src/DeepNTuples/DeepNtuplizer/data/MyMCPileupHist.root")
file_pileup_Data = ROOT.TFile(
    "/afs/cern.ch/work/d/dwalter/DeepNTuples/CMSSW_8_0_29/src/DeepNTuples//DeepNtuplizer/data/MyDataPileupHist.root")

hist_pileup_MC = file_pileup_MC.Get("pileup")
hist_pileup_Data = file_pileup_Data.Get("pileup")

def makePileupWeights(hist_data, hist_mc):
    pupWeights = [1.0]
    hist_mc.Scale(1. / hist_mc.Integral())
    hist_data.Scale(1. / hist_data.Integral())
    for bin in range(1,hist_data.GetNbinsX()+2):
        if hist_mc.GetBinContent(bin) != 0:
            w = hist_data.GetBinContent(bin)/hist_mc.GetBinContent(bin)
            pupWeights.append(w)
        else:
            pupWeights.append(1)
    return pupWeights

print makePileupWeights(hist_pileup_Data,hist_pileup_MC)