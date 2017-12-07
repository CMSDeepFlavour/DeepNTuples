#########
# Compare Jets from MiniAOD with Jets from DeepNtuplizer
#########
import pdb
import ROOT
import os

from DataFormats.FWLite import Events, Handle
from argparse import ArgumentParser


#Dont know if this is necessary
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.gSystem.Load("libDataFormatsPatCandidates.so")


#

jets_MiniAOD = "TT/output_0_1.root"
jets_tuple = ROOT.TFile("output_0.root")


print("create handles")

### Jets
jethandle = Handle('vector<pat::Jet>')
jetlabel = ("GoodJets")
#jetcorrlabel= ("updatedPatJetsUpdatedJEC")

ROOT.gROOT.SetBatch()  # don't pop up canvases

# Create histograms, etc.
print("create hists")

ROOT.gROOT.SetStyle('Plain')  # white background

from array import array

hist_pt_MiniAOD = ROOT.TH1F ("jet_pt_MiniAOD", "pt of  jets", 100,0,700) #, array('d',[0, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_eta_MiniAOD = ROOT.TH1F ("jet_eta", "eta of  jets", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))

hist_pt_ntuple = ROOT.TH1F ("jet_pt_ntuple", "pt of  jets", 100,0,700) #, array('d',[0, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_eta_ntuple = ROOT.TH1F ("jet_eta_ntuple", "eta of  jets", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))

tree = jets_tuple.Get("deepntuplizer/tree")

tree.Draw("jet_pt>>jet_pt_ntuple")
tree.Draw("jet_eta>>jet_eta_ntuple")


hists_MiniAOD = hist_pt_MiniAOD, hist_eta_MiniAOD

def fillHists_MiniAOD(hists, infile, islist = False):
    nevents = 0
    # loop over events
    directory = ''



    if islist is True:
        print("loop over filelist "+ str(infile))
        filelist = open(infile,'r')
        ifile = filelist.readline()
        ifile = ifile[:-1]
        directory = infile[:-len(infile.split('/')[-1])]
    else:
        ifile = infile


    while ifile != '':
        print("process file ",ifile)
        events = Events(directory + ifile)

        for event in events:
            nevents += 1


            # Fill jet hists
            event.getByLabel(jetlabel, jethandle)
            jets = jethandle.product()

            for jet in jets:
                hists[0].Fill(jet.pt())
                hists[1].Fill(jet.eta())

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

    for hist in hists:
        hist.Scale(lumi)                    #scale to luminosity of experiment to compare
        hist.Scale(sigma)                   #Scale for linearity to cross section
        hist.Scale(eff)                     #Scale for linearity to efficiency
        #hist.Scale(1./(hist.Integral()+1))


def makeDataPlots(canvas, hists, name="data"):
    for i in range(0,len(hists)):
        canvas.Clear()
        hists[i].SetTitleFont(42,"t")
        hists[i].SetTitleFont(42,"xyz")
        hists[i].SetLabelFont(42,"xyz")
        hists[i].Draw()
        canvas.Print(name + "_" + str(i) + ".png")

def makeComparePlots(canvas, hist1, hist2, name="comp"):        #Draw two hists in one canvas
    canvas.Clear()
    hist1.SetLineColor(ROOT.kRed)
    hist2.SetLineColor(ROOT.kBlue)
    hist1.Draw()
    hist2.Draw("same")
    canvas.Print(name + ".png")

#MAIN part


fillHists_MiniAOD(hists_MiniAOD, infile=jets_MiniAOD, islist=False)

print("postprocessing")
#ScaleHists(hists, lumi=lumi_analysis/lumi_data)


directory = os.path.dirname('./plots_jets/')
# make a canvas, draw, and save it
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

print("draw and save")

hist_pt_ntuple.SetLineColor(ROOT.kRed)
hist_eta_ntuple.SetLineColor(ROOT.kRed)

c1 = ROOT.TCanvas()
hist_pt_ntuple.Draw()
hist_pt_MiniAOD.Draw("SAME")
c1.Print("Jet_pt.png")

c1.Clear()
hist_pt_MiniAOD.Draw()
c1.Print("Jet_pt_MiniAOD.png")

c1.Clear()
hist_pt_ntuple.Draw()
c1.Print("Jet_pt_nTuple.png")

c1 = ROOT.TCanvas()
hist_eta_ntuple.Draw()
hist_eta_MiniAOD.Draw("SAME")
c1.Print("Jet_eta.png")

c1.Clear()
hist_eta_MiniAOD.Draw()
c1.Print("Jet_eta_MiniAOD.png")

c1.Clear()
hist_eta_ntuple.Draw()
c1.Print("Jet_eta_nTuple.png")
