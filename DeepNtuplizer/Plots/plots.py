#########
#   Make Plots from MiniAOD root file computed with tt_dilep_selector
#########
import pdb
import ROOT
import os
import datetime

import FWCore.ParameterSet.Config as cms

from DataFormats.FWLite import Events, Handle
from argparse import ArgumentParser

print("starting time: "+ str(datetime.datetime.now()))

parser = ArgumentParser('program to make plots from tt_dilep_selector output data')
#parser.add_argument("-i", help="set input file for data", metavar="FILE")
#parser.add_argument("-tt", help="set input file for tt MC", metavar="FILE")
#parser.add_argument("-dy", help="set input file for DY MC", metavar="FILE")

#args=parser.parse_args()
#datafile=args.i
#ttfile = args.tt
#dyfile = args.dy

#Dont know if this is necessary
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.gSystem.Load("libDataFormatsPatCandidates.so")

#
datafile = "data.txt"
ttfilelist = "tt.txt"
dy10to50file ="dy10to50.txt"
dy50filelist = "dy50.txt"
wwfile = "ww.txt"
wzfile = "wz.txt"
zzfile = "zz.txt"
wantitfilelist = "wantit.txt"
wtfilelist = "wt.txt"
wjetsfilelist = "wjets.txt"

print("create handles")
### Muons
# create handle outside of loop
muonhandle = Handle('edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >')
# a label is just a tuple of strings that is initialized just
# like and edm::InputTag
muonlabel = ("GoodMuon")

### Electrons
electronhandle = Handle('edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >')
electronlabel = ("GoodElectron")

### Jets
jethandle = Handle('vector<pat::Jet>')
jetlabel = ("GoodJets")

### LHE weights
LHEweighthandle = Handle('LHEEventProduct')
LHEweightlabel = ("externalLHEProducer")

ROOT.gROOT.SetBatch()  # don't pop up canvases
ROOT.gROOT.SetStyle('Plain')  # white background

# Create histograms, etc.
print("create hists")

from array import array
hist_ll_pt_data = ROOT.TH1F("ll_pt_data", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_data = ROOT.TH1F("ll_eta_data", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_data = ROOT.TH1F("tl_pt_data", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_data = ROOT.TH1F("tl_eta_data", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_data = ROOT.TH1F("number_data", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_tt = ROOT.TH1F("ll_pt_tt", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_tt = ROOT.TH1F("ll_eta_tt", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_tt = ROOT.TH1F("tl_pt_tt", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_tt = ROOT.TH1F("tl_eta_tt", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_tt = ROOT.TH1F("number_tt", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_dy50 = ROOT.TH1F("ll_pt_dy50", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_dy50 = ROOT.TH1F("ll_eta_dy50", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_dy50 = ROOT.TH1F("tl_pt_dy50", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_dy50 = ROOT.TH1F("tl_eta_dy50", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_dy50 = ROOT.TH1F("number_dy50", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_dy10to50 = ROOT.TH1F("ll_pt_dy10to50", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_dy10to50 = ROOT.TH1F("ll_eta_dy10to50", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_dy10to50 = ROOT.TH1F("tl_pt_dy10to50", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_dy10to50 = ROOT.TH1F("tl_eta_dy10to50", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_dy10to50 = ROOT.TH1F("number_dy10to50", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_ww = ROOT.TH1F("ll_pt_ww", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_ww = ROOT.TH1F("ll_eta_ww", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_ww = ROOT.TH1F("tl_pt_ww", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_ww = ROOT.TH1F("tl_eta_ww", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_ww = ROOT.TH1F("number_ww", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_wz = ROOT.TH1F("ll_pt_wz", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_wz = ROOT.TH1F("ll_eta_wz", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_wz = ROOT.TH1F("tl_pt_wz", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_wz = ROOT.TH1F("tl_eta_wz", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_wz = ROOT.TH1F("number_wz", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_zz = ROOT.TH1F("ll_pt_zz", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_zz = ROOT.TH1F("ll_eta_zz", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_zz = ROOT.TH1F("tl_pt_zz", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_zz = ROOT.TH1F("tl_eta_zz", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_zz = ROOT.TH1F("number_zz", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_wantit = ROOT.TH1F("ll_pt_wantit", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_wantit = ROOT.TH1F("ll_eta_wantit", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_wantit = ROOT.TH1F("tl_pt_wantit", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_wantit = ROOT.TH1F("tl_eta_wantit", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_wantit = ROOT.TH1F("number_wantit", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_wt = ROOT.TH1F("ll_pt_wt", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_wt = ROOT.TH1F("ll_eta_wt", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_wt = ROOT.TH1F("tl_pt_wt", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_wt = ROOT.TH1F("tl_eta_wt", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_wt = ROOT.TH1F("number_wt", "number of  jets", 10, -0.5, 9.5)

hist_ll_pt_wjets = ROOT.TH1F("ll_pt_wjets", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta_wjets = ROOT.TH1F("ll_eta_wjets", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt_wjets = ROOT.TH1F("tl_pt_wjets", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta_wjets = ROOT.TH1F("tl_eta_wjets", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_n_wjets = ROOT.TH1F("number_wjets", "number of  jets", 10, -0.5, 9.5)

#hist_jet_pt = ROOT.TH1F("jet_pt", "pt of all jets", 9, array('d', [0, 30, 35, 40, 45, 50, 60, 70, 100, 200]))


hists_data = hist_ll_pt_data,hist_ll_eta_data,hist_tl_pt_data,hist_tl_eta_data,hist_j_n_data
hists_tt = hist_ll_pt_tt,hist_ll_eta_tt,hist_tl_pt_tt,hist_tl_eta_tt,hist_j_n_tt
hists_dy50 = hist_ll_pt_dy50,hist_ll_eta_dy50,hist_tl_pt_dy50,hist_tl_eta_dy50,hist_j_n_dy50
hists_dy10to50 = hist_ll_pt_dy10to50,hist_ll_eta_dy10to50,hist_tl_pt_dy10to50,hist_tl_eta_dy10to50, hist_j_n_dy10to50
hists_ww = hist_ll_pt_ww,hist_ll_eta_ww,hist_tl_pt_ww,hist_tl_eta_ww, hist_j_n_ww
hists_wz = hist_ll_pt_wz,hist_ll_eta_wz,hist_tl_pt_wz,hist_tl_eta_wz, hist_j_n_wz
hists_zz = hist_ll_pt_zz,hist_ll_eta_zz,hist_tl_pt_zz,hist_tl_eta_zz, hist_j_n_zz
hists_wantit = hist_ll_pt_wantit,hist_ll_eta_wantit,hist_tl_pt_wantit,hist_tl_eta_wantit, hist_j_n_wantit
hists_wt = hist_ll_pt_wt,hist_ll_eta_wt,hist_tl_pt_wt,hist_tl_eta_wt, hist_j_n_wt
hists_wjets = hist_ll_pt_wjets,hist_ll_eta_wjets,hist_tl_pt_wjets,hist_tl_eta_wjets, hist_j_n_wjets

print("create histogram stacks")
hs_ll_pt = ROOT.THStack("hs", "ll pt")
hs_ll_eta = ROOT.THStack("hs", "ll eta")
hs_tl_pt = ROOT.THStack("hs", "tl pt")
hs_tl_eta = ROOT.THStack("hs", "tl eta")
hs_j_n = ROOT.THStack("hs", "jet number")

hs = hs_ll_pt, hs_ll_eta, hs_tl_pt, hs_tl_eta, hs_j_n



def fillHists(hists, infile, islist = True, useWeights = False):
    print("start to fill hist: " + str(datetime.datetime.now()))

    nevents = 0;

    directory = ''


    if islist is True:
        print("loop over filelist "+ str(infile))
        filelist = open(infile,'r')
        ifile = filelist.readline()
        ifile = ifile[:-1]
        directory = infile[:-len(infile.split('/')[-1])]
    else:
        ifile = infile

    # loop over events
    while ifile != '':
        print("process file ",ifile)
        events = Events(directory + ifile)

        for event in events:
            nevents += 1
            leadingPt = 10
            leadingEta = 0
            trailingPt = 10
            trailingEta = 0
            weight = 1

            if useWeights == True:
                event.getByLabel( LHEweightlabel, LHEweighthandle )
                weights = LHEweighthandle.product()
                weight = weights.weights()[0].wgt/abs(weights.weights()[0].wgt)
                #pdb.set_trace()
                #print("lhe weights ...")
                #for w in weights.weights():
                #    print(w.id, w.wgt)


            # Fill muon hists
            event.getByLabel(muonlabel, muonhandle)
            muons = muonhandle.product()    # get the product

            for muon in muons:

                if muon.pt() > leadingPt:
                    leadingPt = muon.pt()
                    leadingEta = muon.eta()



            # Fill electron hists
            event.getByLabel(electronlabel, electronhandle)
            electrons = electronhandle.product()
            for electron in electrons:

                if electron.pt() > leadingPt:
                    trailingPt = leadingPt
                    leadingPt = electron.pt()
                    trailingEta = leadingEta
                    leadingEta = electron.eta()
                elif electron.pt() > trailingPt:
                    trailingPt = electron.pt()
                    trailingEta = electron.eta()



            hists[0].Fill(leadingPt, weight)
            hists[1].Fill(leadingEta, weight)
            hists[2].Fill(trailingPt, weight)
            hists[3].Fill(trailingEta, weight)

            # Fill jet hists
            event.getByLabel(jetlabel, jethandle)
            jets = jethandle.product()

            numJets = len(jets)
            hists[4].Fill(numJets, weight)

        if islist is True:
            ifile = filelist.readline()
            ifile = ifile[:-1]
        else:
            break

    if islist is True:
        print("close filelist")
        filelist.close()
    print("processed events ", nevents)



# Divide bin contents of histogram hist by it's bin widths
def divideByBinWidth(hist):
    for ibin in range(0,len(hist)):
        content = hist.GetBinContent(ibin)
        content /= hist.GetBinWidth(ibin)
        hist.SetBinContent(ibin,content)
#
#     #real data

def ScaleHists(hists, sigma=1, lumi = 1, eff=1):

    divideByBinWidth(hists[0])
    divideByBinWidth(hists[2])

    hists[1].Scale(4)
    hists[3].Scale(4)
    for hist in hists:
        hist.Scale(lumi)                    #scale to luminosity of experiment to compare
        hist.Scale(sigma)                   #Scale for linearity to cross section
        hist.Scale(eff)                     #Scale for linearity to efficiency
        #hist.Scale(1./(hist.Integral()+1))

def colorHists(hists, color):
    for hist in hists:
        hist.SetFillColor(color)

def removeStatusBoxes(hists):
    for hist in hists:
        hist.SetStats(ROOT.kFALSE)

def makeStacks(hs, hists_list):
    for i in range(0, len(hs)):
        for j in range(0,len(hists_list)):
            hs[i].Add(hists_list[j][i])


def makeDataPlots(canvas, hists, name="data"):
    for i in range(0,len(hists)):
        canvas.Clear()
        hists[i].SetTitleFont(42,"t")
        hists[i].SetTitleFont(42,"xyz")
        hists[i].SetLabelFont(42,"xyz")
        hists[i].Draw("PE")
        canvas.Print(name + "_" + str(i) + ".png")

def makeMCPlots(canvas, hists, name="mc"):
    for i in range(0, len(hists)):
        #canvas.Clear()

        canvas = ROOT.TCanvas("name"+str(i), "title"+str(i),800,600)
        hists[i].SetTitleFont(42,"t")
        hists[i].SetTitleFont(42,"xyz")
        hists[i].SetLabelFont(42,"xyz")
        hists[i].Draw("HIST")
        canvas.Print(name + "_" + str(i) + ".png")

def makeFullPlots(canvas, h_data,hs_mc, name="all"):
    h_data[0].SetMinimum(0)
    h_data[0].SetMaximum(8000)
    h_data[1].SetMinimum(0)
    h_data[1].SetMaximum(160000)
    h_data[2].SetMinimum(0)
    h_data[2].SetMaximum(18500)
    h_data[3].SetMinimum(0)
    h_data[3].SetMaximum(150000)
    h_data[4].SetMinimum(0)
    h_data[4].SetMaximum(140000)


    leg = ROOT.TLegend(0.59, 0.49, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.AddEntry(hists_data[0], "data", "lep")
    leg.AddEntry(hists_tt[0], "t#bar{t}", "f")
    leg.AddEntry(hists_dy50[0], "DY", "f")
    leg.AddEntry(hists_ww[0], "VV", "f")
    leg.AddEntry(hists_wt[0], "Wt/W#bar{t}", "f")
    leg.AddEntry(hists_wjets[0],"WJets", "f")

    leg21 = ROOT.TLegend(0.19, 0.69, 0.39, 0.89)
    leg21.SetBorderSize(0)
    leg21.SetTextFont(42)
    leg21.AddEntry(hists_data[0], "data", "lep")
    leg21.AddEntry(hists_tt[0], "t#bar{t}", "f")
    leg21.AddEntry(hists_dy50[0], "DY", "f")

    leg22 = ROOT.TLegend(0.69, 0.69, 0.89, 0.89)
    leg22.SetBorderSize(0)
    leg22.SetTextFont(42)
    leg22.AddEntry(hists_ww[0], "VV", "f")
    leg22.AddEntry(hists_wt[0], "Wt/W#bar{t}", "f")
    leg22.AddEntry(hists_wjets[0],"WJets", "f")

    leg.SetFillStyle(0)     #make transparent
    leg21.SetFillStyle(0)
    leg22.SetFillStyle(0)

    h_data[0].SetTitleOffset(1.18,"y")
    for hist in h_data:
        hist.SetTitle("")
        #hist.SetTitleSize(0.05, "xyz")
        hist.SetTitleFont(42,"t")
        hist.SetTitleFont(42,"xyz")
        hist.SetLabelFont(42,"xyz")

    h_data[0].SetXTitle("leading lept. p_{T} (GeV)")
    h_data[1].SetXTitle("leading lept. eta")
    h_data[2].SetXTitle("trailing lept. p_{T} (GeV)")
    h_data[3].SetXTitle("trailing lept. eta")
    h_data[4].SetXTitle("N_{Jet}")

    h_data[0].SetYTitle("Events/GeV")
    h_data[1].SetYTitle("Events/0.25")
    h_data[2].SetYTitle("Events/GeV")
    h_data[3].SetYTitle("Events/0.25")
    h_data[4].SetYTitle("Events")

    ROOT.TGaxis.SetMaxDigits(4);   #Force scientific notation for numbers with more than 4 digits


    canvas = ROOT.TCanvas("final", "title", 800, 600)

    for i in (0,2,4):
        canvas.Clear()
        h_data[i].Draw("PE")
        hs_mc[i].Draw("HIST same")
        h_data[i].Draw("PE same")
        #ROOT.gPad.RedrawAxis()      #draw axis in foreground
        leg.Draw("same")
        canvas.Print(name + "_" + str(i) + ".png")

    for i in (1,3):
        canvas.Clear()
        h_data[i].Draw("PE")
        hs_mc[i].Draw("HIST same")
        h_data[i].Draw("PE same")
        leg21.Draw("same")
        leg22.Draw("same")
        canvas.Print(name + "_" + str(i) + ".png")

    canvasSum = ROOT.TCanvas("finalSum", "title", 800, 600)
    canvasSum.Divide(2,3)

    for i in (0,2,4):
        canvasSum.cd(i+1)
        h_data[i].Draw("PE")
        hs_mc[i].Draw("HIST same")
        h_data[i].Draw("PE same")
        leg.Draw("same")

    for i in (1,3):
        canvasSum.cd(i+1)
        h_data[i].Draw("PE")
        hs_mc[i].Draw("HIST same")
        h_data[i].Draw("PE same")
        leg21.Draw("same")
        leg22.Draw("same")
    canvasSum.Print("alltogether.png")

def computeDiffList(h_data, hlist_mc):
    diff = []
    for ihist in range(0,len(h_data)):              #loop over different hists
        for ibin in range(0,len(h_data[ihist])):    #loop over bins in hist
            n_mc = 0
            n_data = h_data[ihist].GetBinContent(ibin)
            for ihlist in range(0,len(hlist_mc)):     #loop over different mc hists
                n_mc += hlist_mc[ihlist][ihist].GetBinContent(ibin)
            diff.append(n_data/(n_mc+1))
    print("diff: ", diff)



#MAIN part


lumi_analysis = 35.9
lumi_data =      6.615

sigma_tt =            831760
sigma_dy50 =         6025200
sigma_dy10to50 =    22635100
sigma_ww =            118700
sigma_wz =             44900
sigma_zz =             15400
sigma_wantit =         35600
sigma_wt =             35600
sigma_wjets =       61526700

eff_tt = 1./(77081156 + 77867738)
eff_dy50 = 1./122055388
eff_dy10to50 = 1./65888233
eff_ww = 1./994012
eff_wz = 1./1000000
eff_zz = 1./998034
eff_wantit = 1./6933094
eff_wt = 1./6952830
eff_wjets = 1./24120319


fillHists(hists_data, infile=datafile)
fillHists(hists_tt, infile=ttfilelist)
fillHists(hists_dy50, infile=dy50filelist, useWeights =  True)
fillHists(hists_dy10to50, infile=dy10to50file, useWeights = True)
fillHists(hists_ww, infile=wwfile)
fillHists(hists_wz, infile=wzfile)
fillHists(hists_zz, infile=zzfile)
fillHists(hists_wantit, infile=wantitfilelist)
fillHists(hists_wt, infile=wtfilelist)
fillHists(hists_wjets, infile=wjetsfilelist, useWeights = True)


print("postprocessing")
ScaleHists(hists_data,                              lumi=lumi_analysis/lumi_data)
ScaleHists(hists_tt,        sigma = sigma_tt,       lumi = lumi_analysis,   eff = eff_tt)
ScaleHists(hists_dy50,      sigma = sigma_dy50,     lumi = lumi_analysis,   eff = eff_dy50)
ScaleHists(hists_dy10to50,  sigma = sigma_dy10to50, lumi = lumi_analysis,   eff = eff_dy10to50)
ScaleHists(hists_ww,        sigma = sigma_ww,       lumi = lumi_analysis,   eff = eff_ww)
ScaleHists(hists_wz,        sigma = sigma_wz,       lumi = lumi_analysis,   eff = eff_wz)
ScaleHists(hists_zz,        sigma = sigma_zz,       lumi = lumi_analysis,   eff = eff_zz)
ScaleHists(hists_wantit,    sigma = sigma_wantit,   lumi = lumi_analysis,   eff = eff_wantit)
ScaleHists(hists_wt,        sigma = sigma_wt,       lumi = lumi_analysis,   eff = eff_wt)
ScaleHists(hists_wjets,     sigma = sigma_wjets,    lumi = lumi_analysis,   eff = eff_wjets)


colorHists(hists_dy10to50,  ROOT.kBlue)
colorHists(hists_dy50,      ROOT.kBlue)
colorHists(hists_wantit,    ROOT.kMagenta)
colorHists(hists_wt,        ROOT.kMagenta)
colorHists(hists_ww,        ROOT.kYellow)
colorHists(hists_wz,        ROOT.kYellow)
colorHists(hists_zz,        ROOT.kYellow)
colorHists(hists_tt,        ROOT.kRed)
colorHists(hists_wjets,     ROOT.kGreen)

directory = os.path.dirname('./plots/')
# make a canvas, draw, and save it
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

print("draw and save")
c1 = ROOT.TCanvas()

makeDataPlots(c1, hists_data, name="data")
makeMCPlots(c1, hists_tt, name="tt")
makeMCPlots(c1, hists_dy50, name="dy50")
makeMCPlots(c1, hists_dy10to50, name="dy10to50")
makeMCPlots(c1, hists_ww, name = "ww")
makeMCPlots(c1, hists_wz, name = "wz")
makeMCPlots(c1, hists_zz, name = "zz")
makeMCPlots(c1, hists_wantit, name = "wantit")
makeMCPlots(c1, hists_wt, name = "wt")
makeMCPlots(c1, hists_wjets, name = "wjets")

histslist_mc = hists_dy10to50, hists_dy50, hists_wantit, hists_wt, hists_ww, hists_wz, hists_zz, hists_wjets, hists_tt
#histslist_mc = hists_dy50, hists_dy10to50, hists_ww, hists_wz, hists_tt

removeStatusBoxes(hists_data)
for hists in histslist_mc:
    removeStatusBoxes(hists)

makeStacks(hs=hs, hists_list=histslist_mc)

makeFullPlots(c1, h_data=hists_data, hs_mc=hs, name="all")

#computeDiffList(hists_data, histslist_mc)

print("finish at: " + str(datetime.datetime.now()))


