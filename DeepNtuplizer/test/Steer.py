import os
from subprocess import call
from fnmatch import fnmatch

outputdir = "/afs/cern.ch/work/e/ebols/public/hists2/"

datadir = ["JetHT/runB/", "JetHT/runC/", "JetHT/runD/", "JetHT/runE/", "JetHT/runF/"]

mcdir = ["QCD_Pt_30to50_TuneCP5_13TeV_pythia8/incl_qcd_30/180117_082536/0000/", "QCD_Pt_50to80_TuneCP5_13TeV_pythia8/incl_qcd_50/180117_082611/0000/", "QCD_Pt_80to120_TuneCP5_13TeV_pythia8/incl_qcd_80/180117_082645/0000/", "QCD_Pt_120to170_TuneCP5_13TeV_pythia8/incl_qcd_120/180117_082718/0000/", "QCD_Pt_170to300_TuneCP5_13TeV_pythia8/incl_qcd_170/180117_082753/0000/", "QCD_Pt_300to470_TuneCP5_13TeV_pythia8/incl_qcd_300/180117_082827/0000/","QCD_Pt_470to600_TuneCP5_13TeV_pythia8/incl_qcd_470/180117_082905/0000/"]

root = '/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/'
config = '/afs/cern.ch/work/e/ebols/public/CMSSW_9_2_10/src/DeepNTuples/DeepNtuplizer/data/config.txt'


outfile = "file"
outend = ".root"

nr = 0


#replace weightDistribution with ntupleConverter when making ntuples instead of nPV distributions 

for x in datadir:
    for root1, dirs, files in os.walk(root+x):
        for file in files:
            if file.endswith(".root"):
                index = str(nr)
                call(["weightDistribution", os.path.join(root1,file), outputdir+x+outfile+index+outend, "DATA"], config) 
                nr += 1
                if nr%20 == 0:
                    print(nr)


for x in mcdir:
    for file in os.listdir(root+x):
        if file.endswith(".root"):
            index = str(nr)
            call(["weightDistribution", os.path.join(root+x,file), outputdir+x+outfile+index+outend, "MC"], config)
            nr += 1
            if nr%20 == 0:
                print(nr)
