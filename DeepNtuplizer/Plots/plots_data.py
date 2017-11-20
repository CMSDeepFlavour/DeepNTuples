#########
#   Make variety of Plots from MiniAOD root file computed with tt_dilep_selector
#   only from one kind of process
#########
import pdb
import ROOT
import os

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


#

#datafile = "data.txt"
datafile = "output_0.root"


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
jetcorrlabel= ("updatedPatJetsUpdatedJEC")

ROOT.gROOT.SetBatch()  # don't pop up canvases

# Create histograms, etc.
print("create hists")

ROOT.gROOT.SetStyle('Plain')  # white background

from array import array
hist_el_pt = ROOT.TH1F("el_pt_data", "pt ot electrons", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_el_eta = ROOT.TH1F("el_eta_data", "eta of electrons", 20, -2.4, 2.4)
hist_muon_pt = ROOT.TH1F("muon_pt_data", "pt of muons", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_muon_eta = ROOT.TH1F("muon_eta_data", "eta of muons", 20, -2.4, 2.4)
hist_ll_pt = ROOT.TH1F("ll_pt_data", "pt ot the leading lepton", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_eta = ROOT.TH1F("ll_eta_data", "eta of the leading leption", 20, -2.4, 2.4)
hist_tl_pt = ROOT.TH1F("tl_pt_data", "pt ot the trailing lepton", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_eta = ROOT.TH1F("tl_eta_data", "eta of the trailing leption", 20, -2.4, 2.4)
hist_j_pt = ROOT.TH1F ("jet_pt", "pt of  jets", 9, array('d', [0, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_j_eta = ROOT.TH1F ("jet_eta", "eta of  jets", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))
hist_j_n = ROOT.TH1F("jet_n", "number of  jets", 10, -0.5, 9.5)
hist_jcor_pt = ROOT.TH1F ("jetcor_pt", "pt of  jets", 9, array('d', [0, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_jcor_eta = ROOT.TH1F ("jetcor_eta", "eta of  jets", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))
hist_jcor_n = ROOT.TH1F("jetcor_n", "number of  jets", 10, -0.5, 9.5)

hist_ll_muon_pt = ROOT.TH1F("ll_pt_muon_data", "pt ot the leading muon", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_muon_eta = ROOT.TH1F("ll_eta_muon_data", "eta of the leading muon", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))
hist_tl_muon_pt = ROOT.TH1F("tl_pt_muon_data", "pt ot the trailing muon", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_muon_eta = ROOT.TH1F("tl_eta_muon_data", "eta of the trailing muon", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))

hist_ll_el_pt = ROOT.TH1F("ll_pt_el_data", "pt ot the leading electron", 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_ll_el_eta = ROOT.TH1F("ll_eta_el_data", "eta of the leading electron", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))
hist_tl_el_pt = ROOT.TH1F("tl_pt_el_data", "pt ot the trailing electron", 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
hist_tl_el_eta = ROOT.TH1F("tl_eta_el_data", "eta of the trailing electron", 14, array('d', [-2.8,-2.4,-2, -1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,2.8]))


hists = hist_el_pt,hist_el_eta,hist_muon_pt,hist_muon_eta,hist_ll_pt,hist_ll_eta,hist_tl_pt,hist_tl_eta,hist_j_pt,hist_j_eta,hist_j_n, hist_ll_muon_pt, hist_ll_muon_eta, hist_tl_muon_pt, hist_tl_muon_eta, hist_ll_el_pt, hist_ll_el_eta, hist_tl_el_pt, hist_tl_el_eta, hist_jcor_pt, hist_jcor_eta, hist_jcor_n

def fillHists(hists, infile, islist = False):
    nevents = 0;
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
            leadingPt = 10
            leadingEta = 0
            trailingPt = 10
            trailingEta = 0
            leadingIsMuon = True

            # Fill muon hists
            event.getByLabel(muonlabel, muonhandle)
            muons = muonhandle.product()    # get the product

            for muon in muons:
                hists[2].Fill(muon.pt())
                hists[3].Fill(muon.eta())
                if muon.pt() > leadingPt:
                    leadingPt = muon.pt()
                    leadingEta = muon.eta()
                    leadingIsMuon = True


        #    numMuons = len(muons)
        #    hist_Muon_number.Fill(numMuons)

            # Fill electron hists
            event.getByLabel(electronlabel, electronhandle)
            electrons = electronhandle.product()
            for electron in electrons:
                hists[0].Fill(electron.pt())
                hists[1].Fill(electron.eta())
                if electron.pt() > leadingPt:
                    trailingPt = leadingPt
                    leadingPt = electron.pt()
                    trailingEta = leadingEta
                    leadingEta = electron.eta()
                    leadingIsMuon = False
                elif electron.pt() > trailingPt:
                    trailingPt = electron.pt()
                    trailingEta = electron.eta()

        #    numElectrons = len(electrons)
        #    hist_Electron_number.Fill(numElectrons)

            hists[4].Fill(leadingPt)
            hists[5].Fill(leadingEta)
            hists[6].Fill(trailingPt)
            hists[7].Fill(trailingEta)

            if leadingIsMuon == True:
                hists[11].Fill(leadingPt)
                hists[12].Fill(leadingEta)
                hists[17].Fill(trailingPt)
                hists[18].Fill(trailingEta)
            else:
                hists[15].Fill(leadingPt)
                hists[16].Fill(leadingEta)
                hists[13].Fill(trailingPt)
                hists[14].Fill(trailingEta)

            # Fill jet hists
            event.getByLabel(jetlabel, jethandle)
            jets = jethandle.product()

            for jet in jets:
                hists[8].Fill(jet.pt())
                hists[9].Fill(jet.eta())
            numJets = len(jets)
            hists[10].Fill(numJets)

            event.getByLabel(jetcorrlabel, jethandle)
            corrjets = jethandle.product()
            for jet in corrjets:
                hists[19].Fill(jet.pt())
                hists[20].Fill(jet.eta())
            numCorrJets = len(corrjets)
            hists[21].Fill(numCorrJets)

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
    divideByBinWidth(hists[4])
    divideByBinWidth(hists[6])
    divideByBinWidth(hists[8])
    divideByBinWidth(hists[11])
    divideByBinWidth(hists[13])
    divideByBinWidth(hists[15])
    divideByBinWidth(hists[17])
    divideByBinWidth(hists[19])


    hists[1].Scale(4)
    hists[3].Scale(4)
    hists[5].Scale(4)
    hists[7].Scale(4)
    hists[12].Scale(4)
    hists[14].Scale(4)
    hists[16].Scale(4)
    hists[18].Scale(4)

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


lumi_analysis = 35.9
lumi_data =      6.615


fillHists(hists, infile=datafile, islist=False)


print("postprocessing")
ScaleHists(hists, lumi=lumi_analysis/lumi_data)


directory = os.path.dirname('./plots_data_upd/')
# make a canvas, draw, and save it
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

print("draw and save")
c1 = ROOT.TCanvas()

makeComparePlots(c1, hists[8], hists[19], "jet_pt")
makeComparePlots(c1, hists[9], hists[20], "jet_eta")
makeComparePlots(c1, hists[10], hists[21], "jet_n")

makeDataPlots(c1, hists, name="data")





