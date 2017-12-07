#########
# Compare Jets from MiniAOD with Jets from DeepNtuplizer
#########
import pdb
import ROOT
import os
from array import array

class HistStacks:
    """ DOC """
    def __init__(self,name):
        self.hists = []
        self.dataHist = []
        self.minY = 0
        self.maxY = 0
        self.legX = 0.59
        self.legY = 0.49
        self.setDefault = True
        self.name = name
        self.stack = ROOT.THStack("stack_"+self.name, "stack of "+self.name+" hists")

    def addHist(self, hist):
        self.hists.append(hist)

    def addDataHist(self, hist):
        self.dataHist = hist

    def makeStack(self):
        for hist in self.hists:
            self.stack.Add(hist)

    def setMinMax(self, min, max):
        self.minY = min
        self.maxY = max
        self.setDefault = False

    def setLegPos(self, x,y):
        self.legX = x
        self.legY = y

    def removeStatusBoxes(self):
        self.dataHist.SetStats(ROOT.kFALSE)
        for hist in self.hists:
            hist.SetStats(ROOT.kFALSE)

    def drawStack(self):

        self.removeStatusBoxes()
        self.makeStack()

        leg = ROOT.TLegend(self.legX, self.legY, self.legX + 0.2, self.legY + 0.4)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.AddEntry(self.dataHist, "data", "lep")
        leg.AddEntry(self.hists[0], "t#bar{t}", "f")
        leg.AddEntry(self.hists[1], "DY", "f")
        leg.AddEntry(self.hists[5], "VV", "f")
        leg.AddEntry(self.hists[3], "Wt/W#bar{t}", "f")
        leg.AddEntry(self.hists[8], "WJets", "f")

        if not self.setDefault:
            self.dataHist.SetMinimum(self.minY)
            self.dataHist.SetMaximum(self.maxY)

        canvas = ROOT.TCanvas()
        self.dataHist.Draw("AXIS")
        self.stack.Draw("HIST same")
        self.dataHist.Draw("PE same")
        leg.Draw("same")
        ROOT.gPad.RedrawAxis()      #draw axis in foreground
        canvas.Print("stack_"+self.name+".png")


class Source:
    """ DOC """
    stack_pt = HistStacks("pt")
    stack_eta = HistStacks("eta")
    stack_nCpfcand = HistStacks("nCpfcand")
    stack_nNpfcand = HistStacks("nNpfcand")
    stack_Cpfcan_ptrel = HistStacks("Cpfcan_ptrel")
    stack_nsv = HistStacks("nsv")
    stack_npv = HistStacks("npv")
    stack_sv_ntracks = HistStacks("sv_ntracks")
    stack_sv_pt = HistStacks("sv_pt")

    nMC = 0
    nData = 0


    def __init__(self, name, color = ROOT.kBlack, min_weight=0., max_weight=10. , data=False):
        self.name = name
        self.data = data
        self.rootfile = ROOT.TFile(self.name+"_0.root")
        self.tree = self.rootfile.Get("deepntuplizer/tree")

        self.hist_jet_pt = ROOT.TH1F("jet_pt_"+self.name, "pt of "+self.name+" jets", 50, 0, 700)
        self.hist_jet_eta = ROOT.TH1F("jet_eta_"+self.name, "eta of "+self.name+" jets", 50,-2.5,2.5)
        self.hist_nCpfcand = ROOT.TH1F("nCpfcand"+self.name, "number of charged particles in "+self.name, 50,-0.5,49.5)
        self.hist_nNpfcand = ROOT.TH1F("nNpfcand"+self.name, "number of neutral particles in "+self.name, 50,-0.5,49.5)
        self.hist_Cpfcan_ptrel = ROOT.TH1F("Cpfcan_ptrel"+self.name, "pt of charged particles in "+self.name, 30,-1.5,0.5)
        self.hist_nsv = ROOT.TH1F("nsv"+self.name, "number of secondary vertices "+self.name, 5,-0.5,4.5)
        self.hist_npv = ROOT.TH1F("npv"+self.name, "number of primary vertices "+self.name, 50,-.5,49.5)
        self.hist_sv_ntracks = ROOT.TH1F("sv_ntracks"+self.name, "number of primary vertices "+self.name, 20,-.5,19.5)
        self.hist_sv_pt = ROOT.TH1F("sv_pt"+self.name, "number of primary vertices "+self.name, 50,0,200)
        self.hist_weights = ROOT.TH1F("jet_weight_"+self.name, "weight of "+self.name+" jets", 50, min_weight, max_weight)


        self.tree.Draw("jet_pt>>jet_pt_"+self.name, "jet_weight")
        self.tree.Draw("jet_eta>>jet_eta_"+self.name, "jet_weight")
        self.tree.Draw("nCpfcand>>nCpfcand"+self.name, "jet_weight")
        self.tree.Draw("nNpfcand>>nNpfcand"+self.name, "jet_weight")
        self.tree.Draw("Cpfcan_ptrel>>Cpfcan_ptrel"+self.name, "jet_weight")
        self.tree.Draw("nsv>>nsv"+self.name, "jet_weight")
        self.tree.Draw("npv>>npv"+self.name, "jet_weight")
        self.tree.Draw("sv_ntracks>>sv_ntracks"+self.name, "jet_weight")
        self.tree.Draw("sv_pt>>sv_pt"+self.name, "jet_weight")
        self.tree.Draw("jet_weight>>jet_weight_"+self.name)


        self.hists = self.hist_jet_pt, self.hist_jet_eta, self.hist_nCpfcand, self.hist_nNpfcand, self.hist_Cpfcan_ptrel, self.hist_nsv, self.hist_npv, self.hist_sv_ntracks, self.hist_sv_pt

        self.setHistColor(color)

        if data:
            Source.stack_pt.addDataHist(self.hist_jet_pt)
            Source.stack_eta.addDataHist(self.hist_jet_eta)
            Source.stack_nCpfcand.addDataHist(self.hist_nCpfcand)
            Source.stack_nNpfcand.addDataHist(self.hist_nNpfcand)
            Source.stack_Cpfcan_ptrel.addDataHist(self.hist_Cpfcan_ptrel)
            Source.stack_nsv.addDataHist(self.hist_nsv)
            Source.stack_npv.addDataHist(self.hist_npv)
            Source.stack_sv_ntracks.addDataHist(self.hist_sv_ntracks)
            Source.stack_sv_pt.addDataHist(self.hist_sv_pt)

        else:
            self.scaleHists(0.8180639286773081)

            Source.stack_pt.addHist(self.hist_jet_pt)
            Source.stack_eta.addHist(self.hist_jet_eta)
            Source.stack_nCpfcand.addHist(self.hist_nCpfcand)
            Source.stack_nNpfcand.addHist(self.hist_nNpfcand)
            Source.stack_Cpfcan_ptrel.addHist(self.hist_Cpfcan_ptrel)
            Source.stack_nsv.addHist(self.hist_nsv)
            Source.stack_npv.addHist(self.hist_npv)
            Source.stack_sv_ntracks.addHist(self.hist_sv_ntracks)
            Source.stack_sv_pt.addHist(self.hist_sv_pt)


        self.integral = self.hist_jet_pt.Integral()
        if data:
            Source.nData += self.integral
        else:
            Source.nMC += self.integral


        print("Integral of " + self.name + " is " + str(self.integral))

    def setRange_weight(self,min,max):
        self.hist_weights.SetAxisRange(min, max, "X")

    def setHistColor(self,color=ROOT.kBlack):
        for hist in self.hists:
            hist.SetFillColor(color)
        self.hist_weights.SetLineColor(color)

    def scaleHists(self,factor):
        for hist in self.hists:
            hist.Scale(factor)
        self.hist_weights.Scale(factor)

    def draw_all(self):
        i = 0
        for hist in self.hists:
            c1 = ROOT.TCanvas()
            if self.data:
                hist.Draw("PE")
            else:
                hist.Draw("HIST")
            c1.Print(self.name + "_Jet_"+str(i)+".png")
            i += 1

        c1 = ROOT.TCanvas()
        self.hist_weights.Draw()
        c1.Print(self.name + "_Jet_weight.png")


    @classmethod
    def getN_MC(cls):
        return cls.nMC

    @classmethod
    def getN_Data(cls):
        return cls.nData


    @classmethod
    def draw_stacks(cls):

        cls.stack_nNpfcand.setMinMax(0,50000)
        cls.stack_sv_pt.setMinMax(0,70000)
        cls.stack_sv_ntracks.setMinMax(0,200000)

        cls.stack_eta.setLegPos(0.4,0.2)

        cls.stack_pt.drawStack()
        cls.stack_eta.drawStack()
        cls.stack_nCpfcand.drawStack()
        cls.stack_nNpfcand.drawStack()
        cls.stack_Cpfcan_ptrel.drawStack()
        cls.stack_nsv.drawStack()
        cls.stack_npv.drawStack()
        cls.stack_sv_ntracks.drawStack()
        cls.stack_sv_pt.drawStack()




ROOT.gROOT.SetBatch()           # don't pop up canvases
ROOT.gROOT.SetStyle('Plain')    # white background
#ROOT.gStyle.SetFillStyle(0)     # TPave objects (e.g. legend) are transparent

ROOT.gStyle.SetTextFont(42)
ROOT.gStyle.SetTitleFont(42, "t")
ROOT.gStyle.SetTitleFont(42, "xyz")
ROOT.gStyle.SetLabelFont(42, "xyz")

globalScaleFactor = 0.8180639286773081

data = Source("data", data=True)
tt = Source("tt",ROOT.kRed, 0., 1.)
dy50 = Source("dy50",ROOT.kBlue, -5.,5.)
dy10to50 = Source("dy10to50",ROOT.kBlue, -25.,25.)
wantit = Source("wantit",ROOT.kMagenta,0.,1.)
wt = Source("wt",ROOT.kMagenta,0.,1.)
ww = Source("ww",ROOT.kYellow)
wz = Source("wz",ROOT.kYellow)
zz = Source("zz",ROOT.kYellow)
wjets = Source("wjets",ROOT.kGreen,-1000,1000)

print("total number on Jets in MC ", Source.getN_MC(), " and in Data ", Source.getN_Data(), " (Data-MC)/Data) ", (Source.getN_Data() - Source.getN_MC())/Source.getN_Data())



print("make directory and save plots")
directory = os.path.dirname('./plots_jets/')
# make a canvas, draw, and save it
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

data.draw_all()
tt.draw_all()
dy50.draw_all()
dy10to50.draw_all()
wantit.draw_all()
wt.draw_all()
ww.draw_all()
wz.draw_all()
zz.draw_all()
wjets.draw_all()

Source.draw_stacks()


