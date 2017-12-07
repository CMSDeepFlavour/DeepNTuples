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
        self.stack = ROOT.THStack("stack_"+self.name, "jets with label "+self.name)

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
        if self.dataHist != []:
            self.dataHist.SetStats(ROOT.kFALSE)
        for hist in self.hists:
            hist.SetStats(ROOT.kFALSE)

    def drawStack(self):

        self.removeStatusBoxes()
        self.makeStack()

        leg = ROOT.TLegend(self.legX, self.legY, self.legX + 0.2, self.legY + 0.4)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        #leg.AddEntry(self.dataHist, "data", "lep")
        leg.AddEntry(self.hists[0], "t#bar{t}", "f")
        leg.AddEntry(self.hists[1], "DY", "f")
        leg.AddEntry(self.hists[5], "VV", "f")
        leg.AddEntry(self.hists[3], "Wt/W#bar{t}", "f")
        leg.AddEntry(self.hists[8], "WJets", "f")

        if not self.setDefault:
            self.dataHist.SetMinimum(self.minY)
            self.dataHist.SetMaximum(self.maxY)

        canvas = ROOT.TCanvas()
        self.stack.Draw("HIST")
        leg.Draw("same")
        ROOT.gPad.RedrawAxis()      #draw axis in foreground
        canvas.Print("stack_"+self.name+".png")


class Source:
    """ DOC """
    stack_isB = HistStacks("isB")
    stack_isWithB = HistStacks("isWithB")

    def __init__(self, name, color = ROOT.kBlack, data=False):
        self.name = name
        self.data = data
        self.isWithB = 0
        self.isWithoutB = 0
        self.hists = []
        self.rootfile = ROOT.TFile(self.name+"_0.root")
        self.tree = self.rootfile.Get("deepntuplizer/tree")

        self.hist_isB = ROOT.TH1F("isB_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isBB = ROOT.TH1F("isBB_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isGBB = ROOT.TH1F("isGBB_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isLeptonicB = ROOT.TH1F("isLeptonicB_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isLeptonicB_C = ROOT.TH1F("isLeptonicB_C_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isC = ROOT.TH1F("isC_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isGCC = ROOT.TH1F("isGCC_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isCC = ROOT.TH1F("isCC_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isUD = ROOT.TH1F("isUD_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isS = ROOT.TH1F("isS_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isG = ROOT.TH1F("isG_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)
        self.hist_isUndefined = ROOT.TH1F("isUndefined_"+self.name, "labels of "+self.name+" jets", 2, -0.5, 1.5)

        self.hist_isWithB = ROOT.TH1F("isWithB_"+self.name, "labels of "+self.name+" jets is sth with b", 2, -0.5, 1.5)

        self.tree.Draw("isB>>isB_"+self.name, "jet_weight")
        self.tree.Draw("isBB>>isBB_"+self.name, "jet_weight")
        self.tree.Draw("isGBB>>isGBB_"+self.name, "jet_weight")
        self.tree.Draw("isLeptonicB>>isLeptonicB_"+self.name, "jet_weight")
        self.tree.Draw("isLeptonicB_C>>isLeptonicB_C_"+self.name, "jet_weight")
        self.tree.Draw("isC>>isC_"+self.name, "jet_weight")
        self.tree.Draw("isGCC>>isGCC_"+self.name, "jet_weight")
        self.tree.Draw("isCC>>isCC_"+self.name, "jet_weight")
        self.tree.Draw("isUD>>isUD_"+self.name, "jet_weight")
        self.tree.Draw("isS>>isS_"+self.name, "jet_weight")
        self.tree.Draw("isG>>isG_"+self.name, "jet_weight")
        self.tree.Draw("isUndefined>>isUndefined_"+self.name, "jet_weight")


        self.isWithB += self.hist_isB.GetBinContent(2)
        self.isWithB += self.hist_isBB.GetBinContent(2)
        self.isWithB += self.hist_isGBB.GetBinContent(2)
        self.isWithB += self.hist_isLeptonicB.GetBinContent(2)
        self.isWithB += self.hist_isLeptonicB_C.GetBinContent(2)

        self.isWithoutB += self.hist_isC.GetBinContent(2)
        self.isWithoutB += self.hist_isGCC.GetBinContent(2)
        self.isWithoutB += self.hist_isCC.GetBinContent(2)
        self.isWithoutB += self.hist_isUD.GetBinContent(2)
        self.isWithoutB += self.hist_isS.GetBinContent(2)
        self.isWithoutB += self.hist_isG.GetBinContent(2)
        self.isWithoutB += self.hist_isUndefined.GetBinContent(2)

        self.hist_isWithB.Fill(1,self.isWithB)
        self.hist_isWithB.Fill(0,self.isWithoutB)

        self.hists = self.hist_isB, self.hist_isBB, self.hist_isGBB,self.hist_isLeptonicB,self.hist_isLeptonicB_C,\
                     self.hist_isC,self.hist_isGCC,self.hist_isCC,self.hist_isUD,self.hist_isS,self.hist_isG,\
                     self.hist_isUndefined, self.hist_isWithB


        self.setHistColor(color)

        self.scaleHists(0.8180639286773081)


        Source.stack_isB.addHist(self.hist_isB)
        Source.stack_isWithB.addHist(self.hist_isWithB)


    def setHistColor(self,color=ROOT.kBlack):
        for hist in self.hists:
            hist.SetFillColor(color)

    def scaleHists(self,factor):
        for hist in self.hists:
            hist.Scale(factor)

    def draw_all(self):
        i = 0
        for hist in self.hists:
            c1 = ROOT.TCanvas()
            hist.SetMinimum(0)
            if self.data:
                hist.Draw("PE")
            else:
                hist.Draw("HIST")
            c1.Print(self.name + "_Jet_"+str(i)+".png")
            i += 1


    @classmethod
    def draw_stacks(cls):

        cls.stack_isB.setLegPos(0.2,0.2)
        cls.stack_isWithB.setLegPos(0.6,0.2)

        cls.stack_isB.drawStack()
        cls.stack_isWithB.drawStack()






ROOT.gROOT.SetBatch()           # don't pop up canvases
ROOT.gROOT.SetStyle('Plain')    # white background
#ROOT.gStyle.SetFillStyle(0)     # TPave objects (e.g. legend) are transparent

ROOT.gStyle.SetTextFont(42)
ROOT.gStyle.SetTitleFont(42, "t")
ROOT.gStyle.SetTitleFont(42, "xyz")
ROOT.gStyle.SetLabelFont(42, "xyz")


#data = Source("data", data=True)
tt = Source("tt",ROOT.kRed)
dy50 = Source("dy50",ROOT.kBlue)
dy10to50 = Source("dy10to50",ROOT.kBlue)
wantit = Source("wantit",ROOT.kMagenta)
wt = Source("wt",ROOT.kMagenta)
ww = Source("ww",ROOT.kYellow)
wz = Source("wz",ROOT.kYellow)
zz = Source("zz",ROOT.kYellow)
wjets = Source("wjets",ROOT.kGreen)


print("make directory and save plots")
directory = os.path.dirname('./plots_labels/')
# make a canvas, draw, and save it
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

#data.draw_all()
#tt.draw_all()
#dy50.draw_all()
#dy10to50.draw_all()
#wantit.draw_all()
#wt.draw_all()
#ww.draw_all()
#wz.draw_all()
#zz.draw_all()
#wjets.draw_all()

Source.draw_stacks()


