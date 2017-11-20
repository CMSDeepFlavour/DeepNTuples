#########
#   Make Plots from MiniAOD root file computed with tt_dilep_selector
#########

import ROOT
import os
import matplotlib.pyplot as plt
import numpy as np

from DataFormats.FWLite import Events, Handle
from argparse import ArgumentParser

parser = ArgumentParser('program to make plots from tt_dilep_selector output data')
parser.add_argument("-i", help="set input file for data", metavar="FILE")
parser.add_argument("-tt", help="set input file for tt MC", metavar="FILE")
parser.add_argument("-dy", help="set input file for DY MC", metavar="FILE")

args=parser.parse_args()
datafile=args.i
ttfile = args.tt
dyfile = args.dy

#Dont know if this is necessary
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.gSystem.Load("libDataFormatsPatCandidates.so")

datafile = "data.root"
#ttfile =
#ttbackupfile =
dy10to50file ="dy10to50.root"
dy50filelist = "./DY50/filelist.txt"

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


def makeHistsContent(infile, islist = False):
    nevents = 0;
    ll_pt = np.array([])
    ll_eta = np.array([])
    tl_pt = np.array([])
    tl_eta = np.array([])
    j_n = np.array([])

    directory = ''

    ROOT.gROOT.SetBatch()  # don't pop up canvases

    # Create histograms, etc.
    print("create hists")
    ROOT.gROOT.SetStyle('Plain')  # white background

    if islist is True:
        filelist = open(infile,'r')
        ifile = filelist.readline()
        ifile = ifile[:-1]
        directory = infile[:-len(infile.split('/')[-1])]
    else:
        ifile = infile

    print("loop over files")
    while ifile != '':
        print("process file ",ifile)
        events = Events(directory + ifile)

        for event in events:
            nevents += 1
            leadingPt = 10
            leadingEta = 0
            trailingPt = 10
            trailingEta = 0

            # Fill muon hists
            event.getByLabel(muonlabel, muonhandle)
            muons = muonhandle.product()    # get the product

            for muon in muons:

                if muon.pt() > leadingPt:
                    trailingPt = leadingPt
                    leadingPt = muon.pt()
                    trailingEta = leadingEta
                    leadingEta = muon.eta()
                elif muon.pt() > trailingPt:
                    trailingPt = muon.pt()
                    trailingEta = muon.eta()


            # Fill electron hists
            event.getByLabel(electronlabel, electronhandle)
            electrons = electronhandle.product()
            for electron in electrons:
        #        hist_Electron_pt.Fill(electron.pt())
        #        hist_Electron_eta.Fill(electron.eta())
                if electron.pt() > leadingPt:
                    trailingPt = leadingPt
                    leadingPt = electron.pt()
                    trailingEta = leadingEta
                    leadingEta = electron.eta()
                elif electron.pt() > trailingPt:
                    trailingPt = electron.pt()
                    trailingEta = electron.eta()


            ll_pt = np.append(ll_pt,leadingPt)
            ll_eta = np.append(ll_eta,leadingEta)
            tl_pt = np.append(tl_pt,trailingPt)
            tl_eta = np.append(tl_eta,trailingEta)

            # Fill jet hists
            event.getByLabel(jetlabel, jethandle)
            jets = jethandle.product()
        #    for jet in jets:
        #        hist_Jets_pt.Fill(jet.pt())
        #        hist_Jets_eta.Fill(jet.eta())
            numJets = len(jets)
            j_n = np.append(j_n,numJets)

        if islist is True:
            ifile = filelist.readline()
            ifile = ifile[:-1]
        else:
            break
    print("processed events ", nevents)

    return ll_pt, ll_eta, tl_pt, tl_eta, j_n




#htt_ll_pt, htt_ll_eta, htt_tl_pt, htt_tl_eta, htt_j_n = makeHists(ttfile)

#lists_dy50 = makeHistsContent(dy50filelist, True)

lists_dy10to50 = makeHistsContent(dy10to50file, False)

import pdb

y, binEdges = np.histogram(lists_dy10to50[0], bins=[0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200])
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
meanStd     = np.sqrt(y)
width      = binEdges[1:]-binEdges[:-1]
y /= width
pdb.set_trace()

plt.bar(bincenters, y, width=width, align='center', color='r', yerr=meanStd)

#hist_ll_eta = plt.hist(lists_dy50[1], bins=20, range=(-2.4, 2.4))
#hist_tl_pt = plt.hist(lists_dy50[2], bins=[0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200])
#hist_tl_eta = plt.hist(lists_dy50[3], bins=20, range=(-2.4, 2.4))
#hist_j_n = plt.hist(lists_dy50[4], bins=10, range=(-0.5, 9.5))

plt.savefig('test.png')

#hists_data = makeHists(datafile, False)



# # #Postprocessing the histograms
# print("postprocessing")
# #
# # Divide bin contents of histogram hist by it's bin widths
# def divideByBinWidth(hist):
#     for ibin in range(0,len(hist)):
#         content = hist.GetBinContent(ibin)
#         content /= hist.GetBinWidth(ibin)
#         hist.SetBinContent(ibin,content)
# #
# #     #real data
# divideByBinWidth(hists_data[0])
# divideByBinWidth(hists_data[2])
# #
#     #Scale the histogram to a given luminosity
# lumi_analysis = 35.9
# lumi_data = 6.615
# scaleFactor_data = lumi_analysis/lumi_data
#
# hists_data[1].Scale(4)
# hists_data[3].Scale(4)
# for hist in hists_data:
#     hist.Scale(scaleFactor_data)    #scale to luminosity of experiment to compare
#
#
#     #MC
#     #corss sections in fb
# efficiency = 0.7
# sigma_tt = 831.76
# sigma_dy50 = 6025.2
# sigma_dy10to50 = 22635.1
#
# divideByBinWidth(hists_dy50[0])
# divideByBinWidth(hists_dy50[2])
# for hist in hists_dy50:
#     hist.Scale(1./hist.GetEntries())     #normalize
#     hist.Scale(sigma_dy50)              #Scale for linearity to cross section
#     hist.Scale(lumi_analysis)           #Scale for linearity to luminosity
#     hist.Scale(efficiency)              #Scale for linearity to efficiency
#
# #hdy_ll_pt, hdy_ll_eta, hdy_tl_pt, hdy_tl_eta, hdy_j_n = hists_dy50
#
# divideByBinWidth(hists_dy10to50[0])
# divideByBinWidth(hists_dy10to50[2])
# for hist in hists_dy10to50:
#     hist.Scale(1./hist.GetEntries())     #normalize
#     hist.Scale(sigma_dy10to50)          #Scale for linearity to cross section
#     hist.Scale(lumi_analysis)           #Scale for linearity to luminosity
#     hist.Scale(efficiency)              #Scale for linearity to efficiency
#
# #hdy1_ll_pt, hdy1_ll_eta, hdy1_tl_pt, hdy1_tl_eta, hdy1_j_n = hists_dy10to50
#
#
# #create histogram stacks
# print("create histogram stacks")
# hs_ll_pt = ROOT.THStack("hs", "ll pt")
# hs_ll_eta = ROOT.THStack("hs", "ll eta")
# hs_tl_pt = ROOT.THStack("hs", "tl pt")
# hs_tl_eta = ROOT.THStack("hs", "tl eta")
# hs_j_n = ROOT.THStack("hs","jet number")
#
# hs = hs_ll_pt, hs_ll_eta, hs_tl_pt, hs_tl_eta, hs_j_n
#
# for i in range(0,len(hs)):
#     hists_dy50[i].SetFillColor(ROOT.kRed)
#     hs[i].Add(hists_dy50[i])
#     hists_dy10to50[i].SetFillColor(ROOT.kBlue)
#     hs[i].Add(hists_dy10to50[i])
#
#
#
# directory = os.path.dirname('./plots/')
# if not os.path.exists(directory):
#     os.makedirs(directory)
# os.chdir(directory)
#
# # make a canvas, draw, and save it
# print("draw and save")
# c1 = ROOT.TCanvas()
#
#
#
# for i in range(0,len(hs)):
#     c1.Clear
#     hs[i].Draw()
#     hists_data[i].Draw("SAME PE")
#     c1.Print("dy_"+str(i)+".png")
#
# for i in range(0,len(hs)):
#     c1.Clear
#     hists_data[i].Draw()
#     c1.Print("data_"+str(i)+".png")
