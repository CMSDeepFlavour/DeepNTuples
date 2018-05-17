
import pdb
from ROOT import TH1F, TFile, TCanvas, gROOT, gStyle, THStack, TLegend, TGaxis, gPad, TPad, TGraph, TLine
import numpy as np
import os
from array import array
import datetime
from DataFormats.FWLite import Events, Handle

gROOT.SetBatch()  # don't pop up canvases
gROOT.SetStyle('Plain')  # white background
gStyle.SetFillStyle(0)
TGaxis.SetMaxDigits(4)  # Force scientific notation for numbers with more than 4 digits
#gStyle.SetBorderSize(0)
#gStyle.SetTextFont(43)
#gStyle.SetTextSize(14)


class pileupweight:
    """ DOC """

    def __init__(self, file_pileup_MC, file_pileup_Data):
        hist_mc = file_pileup_MC.Get("pileup")
        hist_data = file_pileup_Data.Get("pileup")
        self.pupWeights = [0]
        hist_mc.Scale(1. / hist_mc.Integral())
        hist_data.Scale(1. / hist_data.Integral())
        for bin in range(1, hist_data.GetNbinsX() + 2):
            if hist_mc.GetBinContent(bin) != 0:
                w = hist_data.GetBinContent(bin) / hist_mc.GetBinContent(bin)
                self.pupWeights.append(w)
            else:
                self.pupWeights.append(1)

    def getPileupWeight(self,nInteractions):
        if nInteractions >= len(self.pupWeights):
            return 0
        return self.pupWeights[nInteractions]

class scaleFactor():
    """
     class for isolation and id scale factors
     hist is the 2d histogram with the scale factor with abs(eta) on x and pt on y
    """
    def __init__(self,hist):
        self.hist = hist
        self.xaxis = self.hist.GetXaxis()
        self.yaxis = self.hist.GetYaxis()
        self.binx = 0
        self.biny = 0
        self.scalefactor = 1.
        self.maxBinX = self.xaxis.GetLast()
        self.maxBinY = self.yaxis.GetLast()
        self.minX = self.xaxis.GetXmin()
        self.maxX = self.xaxis.GetXmax()
        self.minY = self.yaxis.GetXmin()
        self.maxY = self.yaxis.GetXmax()

    def getScalefactor(self,x,y):

        self.binx = self.xaxis.FindBin(x)
        self.biny = self.yaxis.FindBin(y)
        if x >= self.maxX: self.binx = self.maxBinX
        if y >= self.maxY: self.biny = self.maxBinY
        self.scalefactor = self.hist.GetBinContent(self.binx, self.biny)
        return self.scalefactor


class scaleFactor_tracking():
    """
     class for tracking scale factors
     the tracking scale factors are saved in a TGraphAsymmErrors object which has to be written in a histogram first
    """
    def __init__(self,graph, name):
        self.graph = graph
        self.npoints = self.graph.GetN()
        self.x_centers = self.graph.GetX()
        self.y_centers = self.graph.GetY()
        self.x_lows = np.zeros(self.npoints)
        self.x_highs = np.zeros(self.npoints)
        self.y_lows = np.zeros(self.npoints)
        self.y_highs = np.zeros(self.npoints)
        self.x_edges = np.zeros(self.npoints +1)

        for i in range(0,self.npoints):
            self.x_lows[i] = self.graph.GetErrorXlow(i)
            self.x_highs[i] = self.graph.GetErrorXhigh(i)
            self.y_lows[i] = self.graph.GetErrorYlow(i)
            self.y_lows[i] = self.graph.GetErrorYhigh(i)

            self.x_edges[i] = self.x_centers[i] - self.x_lows[i]
        self.x_edges[self.npoints] = self.x_centers[self.npoints - 1] + self.x_highs[self.npoints-1]
        self.hist = TH1F(name, name,self.npoints,self.x_edges)

        for i in range(0,self.npoints):
            self.hist.SetBinContent(i+1, self.y_centers[i])
            self.hist.SetBinError(i+1, max(self.y_lows[i], self.y_highs[i]))

        self.bin = 0
        self.axis = self.hist.GetXaxis()
        self.scalefactor = 1.

    def getScalefactor(self,eta):
        self.bin = self.axis.FindBin(abs(eta))
        self.scalefactor = self.hist.GetBinContent(self.bin)
        return self.scalefactor


class process:
    """ DOC """
    datadir = "/afs/desy.de/user/d/dwalter/CMSSW_8_0_29/src/DeepNTuples/DeepNtuplizer/data/"
    pupw = pileupweight(
        TFile(datadir + "MyMCPileupHist.root"),
        TFile(datadir + "MyDataPileupHist.root")
    )
    # muon scalefactor hists from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonWorkInProgressAndPagResults
    file_sf_emu_trigger = TFile(datadir + "triggerSummary_emu_ReReco2016_ttH.root")
    file_sf_muon_id_GH = TFile(datadir + "EfficienciesAndSF_ID_GH.root")
    file_sf_muon_iso_GH = TFile(datadir + "EfficienciesAndSF_ISO_GH.root")
    file_sf_muon_tracking = TFile(datadir + "Tracking_EfficienciesAndSF_BCDEFGH.root")

    # electron scalefactor hists from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2
    #   -> Electron cut-based 80XID WPs. Scale factors for 80X ->  Tight cut-based ID WP scale factor (root file)
    #   id and iso is combined in one file, there is no tracking for the electron
    file_sf_el_id = TFile(datadir + "egammaEffi.txt_EGM2D.root"
    )

    h_sf_emu_trigger = file_sf_emu_trigger.Get("scalefactor_eta2d_with_syst")
    h_sf_muon_id_GH = file_sf_muon_id_GH.Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")
    h_sf_muon_iso_GH = file_sf_muon_iso_GH.Get("TightISO_TightID_pt_eta/abseta_pt_ratio")
    g_sf_muon_tracking = file_sf_muon_tracking.Get("ratio_eff_aeta_dr030e030_corr")
    h_sf_el_id = file_sf_el_id.Get("EGamma_SF2D")

    sf_emu_trigger = scaleFactor(h_sf_emu_trigger)
    sf_muon_id = scaleFactor(h_sf_muon_id_GH)
    sf_muon_iso = scaleFactor(h_sf_muon_iso_GH)
    sf_muon_tracking = scaleFactor_tracking(g_sf_muon_tracking,"sf_muon_tracking")
    sf_el_id = scaleFactor(h_sf_el_id)

    muonhandle = Handle("vector<pat::Muon>")
#muonhandle = Handle('edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >')
    muonlabel = ("goodMuons")
#muonlabel = ("GoodMuon")
    electronhandle = Handle("vector<pat::Electron>")
#electronhandle = Handle('edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >')
    electronlabel = ("goodElectrons")
#electronlabel = ("GoodElectron")
    jethandle = Handle('vector<pat::Jet>')
    jetlabel = ('slimmedJets')
#jetlabel = ("GoodJets")
    ### LHE weights
    LHEweighthandle = Handle('LHEEventProduct')
    LHEweightlabel = ("externalLHEProducer")

    ## Pileup reweighting
    puphandle = Handle('vector<PileupSummaryInfo>')
    puplabel = ("slimmedAddPileupInfo")

    hs_ll_pt = THStack("hs_ll_pt", "ll pt")
    hs_ll_eta = THStack("hs_ll_eta", "ll eta")
    hs_tl_pt = THStack("hs_tl_pt", "tl pt")
    hs_tl_eta = THStack("hs_tl_eta", "tl eta")
    hs_j_n = THStack("hs_njets", "jet number")

    hs_list = hs_ll_pt, hs_ll_eta, hs_tl_pt, hs_tl_eta, hs_j_n

    collection = []
    hists_data = None

    tsize = 20
    tfont = 43

    def __init__(self, name, isData=False, color=1):
        self.name = name
        self.isData = isData
        self.h_ll_pt = TH1F(self.name + "_ll_pt", "pt of leading lepton in " + self.name, 10, array('d', [0, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
        self.h_ll_eta = TH1F(self.name + "_ll_eta", "eta of the leading leption in " + self.name, 20, -2.4, 2.4)
        self.h_tl_pt = TH1F(self.name + "_tl_pt", "pt of leading lepton in " + self.name, 11, array('d', [0, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 200]))
        self.h_tl_eta = TH1F(self.name + "_tl_eta", "eta of the leading leption in " + self.name, 20, -2.4, 2.4)
        self.h_njets = TH1F(self.name + "_njets", "number of  jets in " +self.name, 10, -0.5, 9.5)

        self.hists = self.h_ll_pt, self.h_ll_eta, self.h_tl_pt, self.h_tl_eta, self.h_njets

        self.integral = 0.

        #remove status boxes
        for hist in self.hists:
            hist.SetStats(False)

        if self.isData==False:
            process.collection.append(self)
            for hist in self.hists:
                hist.SetFillColor(color)
        else:
            process.hists_data = self

            self.h_ll_pt.SetXTitle("leading lept. p_{T} (GeV)")
            self.h_ll_eta.SetXTitle("leading lept. eta")
            self.h_tl_pt.SetXTitle("trailing lept. p_{T} (GeV)")
            self.h_tl_eta.SetXTitle("trailing lept. eta")
            self.h_njets.SetXTitle("N_{Jet}")

            self.h_ll_pt.SetYTitle("Events/GeV")
            self.h_ll_pt.GetYaxis().SetTitleOffset(1.2)
            self.h_ll_eta.SetYTitle("Events/0.25")
            self.h_tl_pt.SetYTitle("Events/GeV")
            self.h_tl_eta.SetYTitle("Events/0.25")
            self.h_njets.SetYTitle("Events")

            self.h_njets.GetXaxis().SetNdivisions(10,0,0)

    def fillHists(self, infile, islist=True, useLHEWeights=False, isData=False, eventIDs=[]):
        print("start to fill hist: " + str(datetime.datetime.now()))
        nevents = 0
        directory = ''

        if islist is True:
            print("loop over filelist " + str(infile))
            filelist = open(infile, 'r')
            ifile = filelist.readline()
            ifile = ifile[:-1]
            directory = infile[:-len(infile.split('/')[-1])]
        else:
            ifile = infile

        while ifile != '':
            print("process file ", ifile)
            events = Events(directory + ifile)

            ### EVENT LOOP
            for event in events:
                eventID = event.eventAuxiliary().event()

                if isData:  # avoid double counting
                    if eventID in eventIDs:
                        # print("this event was already processed")
                        continue
                    eventIDs = np.append(eventIDs, eventID)

                nevents += 1
                leadingPt = 10
                leadingEta = 0
                trailingPt = 10
                trailingEta = 0
                leadingMuonPt = 10
                leadingMuonEta = 0
                leadingElPt = 10
                leadingElEta = 0
                leadingElSuEta = 0
                weight = 1.

                # Fill muon hists
                event.getByLabel(process.muonlabel, process.muonhandle)
                muons = process.muonhandle.product()  # get the product

                for muon in muons:
                    if muon.pt() > leadingMuonPt:
                        leadingMuonPt = muon.pt()
                        leadingMuonEta = muon.eta()

                # Fill electron hists
                event.getByLabel(process.electronlabel, process.electronhandle)
                electrons = process.electronhandle.product()
                for electron in electrons:

                    if electron.pt() > leadingElPt:
                        leadingElPt = electron.pt()
                        leadingElEta = electron.eta()
                        leadingElSuEta = electron.superCluster().eta()

                if leadingMuonPt > leadingElPt:
                    leadingPt = leadingMuonPt
                    leadingEta = leadingMuonEta
                    trailingPt = leadingElPt
                    trailingEta = leadingElEta
                else:
                    leadingPt = leadingElPt
                    leadingEta = leadingElEta
                    trailingPt = leadingMuonPt
                    trailingEta = leadingMuonEta



                # Fill jet hists
                event.getByLabel(process.jetlabel, process.jethandle)
                jets = process.jethandle.product()

                numJets = len(jets)

                # Compute Event weights
                if useLHEWeights == True:
                    event.getByLabel(process.LHEweightlabel, process.LHEweighthandle)
                    weights = process.LHEweighthandle.product()
                    weight = weights.weights()[0].wgt / abs(weights.weights()[0].wgt)

                if isData == False:
                    event.getByLabel(process.puplabel, process.puphandle)
                    pupInfos = process.puphandle.product()
                    for pupInfo in pupInfos:
                        if (pupInfo.getBunchCrossing() == 0):
                            weight *= process.pupw.getPileupWeight(int(pupInfo.getTrueNumInteractions()))

                    weight *= process.sf_emu_trigger.getScalefactor(abs(leadingElEta),abs(leadingMuonEta))
                    weight *= process.sf_muon_id.getScalefactor(abs(leadingMuonEta), leadingMuonPt)
                    weight *= process.sf_muon_iso.getScalefactor(abs(leadingMuonEta), leadingMuonPt)
                    weight *= process.sf_muon_tracking.getScalefactor(abs(leadingMuonEta))
                    weight *= process.sf_el_id.getScalefactor(abs(leadingElSuEta), leadingElPt)


                self.h_ll_pt.Fill(leadingPt, weight)
                self.h_ll_eta.Fill(leadingEta, weight)
                self.h_tl_pt.Fill(trailingPt, weight)
                self.h_tl_eta.Fill(trailingEta, weight)
                self.h_njets.Fill(numJets, weight)

            if islist is True:
                ifile = filelist.readline()
                ifile = ifile[:-1]
            else:
                break

        if islist is True:
            print("close filelist")
            filelist.close()

        print("processed events ", nevents)
        self.save_hists()

    def scale_hists(self,sigma=1.,lumi=1.,eff=1.):
        for hist in self.hists:
            hist.Scale(sigma*lumi*eff)
        self.integral = self.hists[0].Integral()

        self.divideByBinWidth(self.hists[0])
        self.divideByBinWidth(self.hists[2])
        self.hists[1].Scale(4)
        self.hists[3].Scale(4)

    def divideByBinWidth(self,hist):
        for ibin in range(0, len(hist)):
            content = hist.GetBinContent(ibin)
            content /= hist.GetBinWidth(ibin)
            hist.SetBinContent(ibin, content)

    def save_hists(self):
        f1 = TFile(self.name + ".root", "RECREATE")
        for hist in self.hists:
            hist.Write()
        f1.Close()

    def load_hists(self):
        print "load " + self.name

        if(os.path.isfile(self.name + ".root")):
            f1 = TFile(self.name + ".root", "READ")
            self.h_ll_pt.Add(f1.Get(self.name + "_ll_pt"))
            self.h_ll_eta.Add(f1.Get(self.name + "_ll_eta"))
            self.h_tl_pt.Add(f1.Get(self.name + "_tl_pt"))
            self.h_tl_eta.Add(f1.Get(self.name + "_tl_eta"))
            self.h_njets.Add(f1.Get(self.name + "_njets"))

            f1.Close()
        else:
            print("no file for this process")

    def draw_hists(self):
        canvas = TCanvas(self.name, "", 800, 600)
        for hist in self.hists:
            canvas.Clear()
            hist.Draw("HIST")
            canvas.Print(hist.GetName()+".png")

    def setMinMax(self, y_min_ll_pt, y_max_ll_pt, y_min_ll_eta, y_max_ll_eta, y_min_tl_pt, y_max_tl_pt, y_min_tl_eta, y_max_tl_eta, y_min_njets, y_max_njets):
        self.hists[0].SetMinimum(y_min_ll_pt)
        self.hists[0].SetMaximum(y_max_ll_pt)
        self.hists[1].SetMinimum(y_min_ll_eta)
        self.hists[1].SetMaximum(y_max_ll_eta)
        self.hists[2].SetMinimum(y_min_tl_pt)
        self.hists[2].SetMaximum(y_max_tl_pt)
        self.hists[3].SetMinimum(y_min_tl_eta)
        self.hists[3].SetMaximum(y_max_tl_eta)
        self.hists[4].SetMinimum(y_min_njets)
        self.hists[4].SetMaximum(y_max_njets)

    @classmethod
    def make_stacks(cls):
        for proc in cls.collection:
            for i,hist in enumerate(proc.hists):
                cls.hs_list[i].Add(hist)


    @classmethod
    def print_stacks(cls):
        cls.make_stacks()

        leg = TLegend(0.59, 0.54, 0.89, 0.84)
        leg.SetBorderSize(0)
        leg.SetTextFont(cls.tfont)
        leg.AddEntry(cls.hists_data.hists[0], "data", "ep")
        leg.AddEntry(cls.collection[8].hists[0], "t#bar{t}", "f")
        leg.AddEntry(cls.collection[6].hists[0], "Wt/W#bar{t}", "f")
        leg.AddEntry(cls.collection[5].hists[0], "W+Jets", "f")
        leg.AddEntry(cls.collection[2].hists[0], "VV", "f")
        leg.AddEntry(cls.collection[0].hists[0],  "DY", "f")

        leg21 = TLegend(0.19, 0.69, 0.39, 0.84)
        leg21.SetBorderSize(0)
        leg21.SetTextFont(cls.tfont)
        leg21.AddEntry(cls.hists_data.hists[0], "data", "ep")
        leg21.AddEntry(cls.collection[8].hists[0], "t#bar{t}", "f")
        leg21.AddEntry(cls.collection[6].hists[0], "Wt/W#bar{t}", "f")

        leg22 = TLegend(0.69, 0.69, 0.89, 0.84)
        leg22.SetBorderSize(0)
        leg22.SetTextFont(cls.tfont)
        leg22.AddEntry(cls.collection[5].hists[0], "W+Jets", "f")
        leg22.AddEntry(cls.collection[2].hists[0], "VV", "f")
        leg22.AddEntry(cls.collection[0].hists[0],  "DY", "f")

        leg.SetFillStyle(0)  # make transparent
        leg21.SetFillStyle(0)
        leg22.SetFillStyle(0)

        canvas = TCanvas("c1","Example 2 pads (20,80)",200,10,800,600)
        #canvas.SetBottomMargin(0.25)
        pad1 = TPad("pad1", "The pad 80% of the height", 0, 0.2, 1, 1.0)
        pad2 = TPad("pad2", "The pad 20% of the height", 0, 0.05, 1, 0.2)

        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.1)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.1)

        pad1.SetBottomMargin(0.02)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.25)

        pad1.Draw()
        pad2.Draw()
        for i, h_data in enumerate(cls.hists_data.hists):

            h_data.SetTitle("")
            h_data.GetXaxis().SetLabelFont(cls.tfont)
            h_data.GetXaxis().SetLabelSize(0)
            h_data.GetXaxis().SetTitleSize(0)
            h_data.GetYaxis().SetLabelFont(cls.tfont)
            h_data.GetYaxis().SetLabelSize(cls.tsize)
            h_data.GetYaxis().SetTitleFont(cls.tfont)
            h_data.GetYaxis().SetTitleSize(cls.tsize)
            h_data.GetYaxis().SetTitleOffset(1.5)

            h_data.SetMarkerStyle(20)
            #pad1.Clear()
            #pad2.Clear()
            pad1.cd()

            h_data.Draw("E1")
            cls.hs_list[i].Draw("HIST same")
            h_data.Draw("E1 same")
            gPad.RedrawAxis()  # draw axis in foreground

            if i%2 == 0:
                leg.Draw("same")
            else:
                leg21.Draw("same")
                leg22.Draw("same")

            pad2.cd()
            h_ratio = TH1F(cls.hs_list[i].GetStack().Last())

            h_ratio.SetTitle("")
            h_ratio.GetXaxis().SetLabelFont(cls.tfont)
            h_ratio.GetXaxis().SetLabelSize(cls.tsize)
            h_ratio.GetYaxis().SetLabelFont(cls.tfont)
            h_ratio.GetYaxis().SetLabelSize(cls.tsize)
            h_ratio.GetXaxis().SetTitleFont(cls.tfont)
            h_ratio.GetXaxis().SetTitleSize(cls.tsize)
            h_ratio.GetYaxis().SetTitleFont(cls.tfont)
            h_ratio.GetYaxis().SetTitleSize(cls.tsize)
            h_ratio.GetXaxis().SetTitleOffset(7.)
            h_ratio.GetYaxis().SetTitleOffset(1.5)
            h_ratio.GetYaxis().SetNdivisions(403)
            h_ratio.GetXaxis().SetTickLength(0.15)

            if i ==4:
                h_ratio.GetXaxis().SetNdivisions(10, 0, 0)

            h_ratio.SetMarkerStyle(20)
            h_ratio.SetXTitle(h_data.GetXaxis().GetTitle())

            h_ratio.Divide(h_data)
            h_ratio.GetYaxis().SetTitle("mc/data")

            h_ratio.SetMinimum(0.75)
            h_ratio.SetMaximum(1.25)


            middleLine = TLine(h_ratio.GetXaxis().GetXmin(),1. ,h_ratio.GetXaxis().GetXmax(), 1.)

            h_ratio.Draw("e1")
            middleLine.Draw("same")
            gStyle.SetErrorX(0.)
            canvas.Print("stacks_" + str(i) + ".png")



dy10to50 = process("dy10to50", color=593)
dy50 = process("dy50", color=593)
ww = process("ww", color=393)
wz = process("wz", color=393)
zz = process("zz", color=393)
wjets = process("wjets", color=409)
wt = process("wt", color=609)
wantit = process("wantit", color=609)
tt = process("tt", color=625)
data = process("data",isData=True)


dy10to50.fillHists("dy10to50.txt",useLHEWeights=True)
dy50.fillHists("dy50.txt",useLHEWeights=True)
ww.fillHists("ww.txt")
wz.fillHists("wz.txt")
zz.fillHists("zz.txt")
wjets.fillHists("wjets.txt",useLHEWeights=True)
#wt.fillHists("wt.txt")
wantit.fillHists("wantit.txt")
tt.fillHists("tt.txt")
data.fillHists("data.txt",isData=True)

#for ip in process.collection:
#    ip.load_hists()
#process.hists_data.load_hists()

lumi_total = 35.9
lumi_data = 8.746

sigma_tt =            831760
sigma_dy50 =         6025200
sigma_dy10to50 =    22635100
sigma_ww =            118700
sigma_wz =             44900
sigma_zz =             15400
sigma_wantit =         35600
sigma_wt =             35600
sigma_wjets =       61526700

# Efficiency = 1/(number of events before cuts), in case of lheweights its the
# effective number of events (in this case: number of events with positive weights - number of events with negative weights)
eff_tt =        1./(77081156 + 77867738)
#eff_tt =    1./(77013872)
eff_dy50 =       1./81781052        # effective number of events
#eff_dy50 =      1./122055388        # true number of events
eff_dy10to50 =   1./47946519        # effective number of events
#eff_dy10to50 =   1./65888233        # true number of events
eff_ww =           1./994012
eff_wz =          1./1000000
eff_zz =           1./998034
eff_wantit =      1./6933094
eff_wt =          1./6952830
eff_wjets =      1./16497031        # effective number of events
#eff_wjets =      1./24120319        # true number of events

data.scale_hists(                           lumi=lumi_total/lumi_data)
tt.scale_hists(sigma = sigma_tt,            lumi = lumi_total,   eff = eff_tt)
dy50.scale_hists(sigma = sigma_dy50,        lumi = lumi_total,   eff = eff_dy50)
dy10to50.scale_hists(sigma = sigma_dy10to50,lumi = lumi_total,   eff = eff_dy10to50)
ww.scale_hists(sigma = sigma_ww,            lumi = lumi_total,   eff = eff_ww)
wz.scale_hists(sigma = sigma_wz,            lumi = lumi_total,   eff = eff_wz)
zz.scale_hists(sigma = sigma_zz,            lumi = lumi_total,   eff = eff_zz)
wantit.scale_hists(sigma = sigma_wantit,    lumi = lumi_total,   eff = eff_wantit)
wt.scale_hists(sigma = sigma_wt,            lumi = lumi_total,   eff = eff_wt)
wjets.scale_hists(sigma = sigma_wjets,      lumi = lumi_total,   eff = eff_wjets)

data.setMinMax(0,20000,0,300000,0,18500,0,150000,0,140000)
title='180327'

directory = os.path.dirname('./plots_'+title+'/')
# make a canvas, draw, and save it
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

integral_mc = 0.
integral_data = process.hists_data.integral
for ip in process.collection:
    integral_mc += ip.integral
print "integral of mc events is " + str(integral_mc)
print "integral of data events is " + str(integral_data)
print "data/mc = " +str(integral_data/integral_mc)

#for ip in process.collection:
#    ip.draw_hists()
#process.hists_data.draw_hists()

process.print_stacks()