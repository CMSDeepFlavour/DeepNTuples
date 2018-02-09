import os
from subprocess import call
from fnmatch import fnmatch

outputdir = "/afs/cern.ch/work/e/ebols/public/tuples/"

datadir = ["JetHT/runB/171220_125115/0000", "JetHT/runC/171220_125350/0000", "JetHT/runC/171220_125350/0001", "JetHT/runD/171220_125424/0000", "JetHT/runE/171220_125459/0000", "JetHT/runE/171220_125459/0001", "JetHT/runF/171220_125533/0000", "JetHT/runF/171220_125533/0001"]

mcdir = ["QCD_Pt_30to50_TuneCP5_13TeV_pythia8/incl_qcd_30/180117_082536/0000/", "QCD_Pt_50to80_TuneCP5_13TeV_pythia8/incl_qcd_50/180117_082611/0000/", "QCD_Pt_80to120_TuneCP5_13TeV_pythia8/incl_qcd_80/180117_082645/0000/", "QCD_Pt_120to170_TuneCP5_13TeV_pythia8/incl_qcd_120/180117_082718/0000/", "QCD_Pt_170to300_TuneCP5_13TeV_pythia8/incl_qcd_170/180117_082753/0000/", "QCD_Pt_300to470_TuneCP5_13TeV_pythia8/incl_qcd_300/180117_082827/0000/"]

root = '/eos/cms/store/group/phys_btag/gpaspala/Commissioning_2018/'

outfile = "data"
outend = ".root"

blah = 0

for root, dirs, files in os.walk(root+"JetHT/runB"):
    for file in files:
        if file.endswith(".root"):
            index = str(blah)
            call(["WeightDistrubtions", os.path.join(root,file), outputdir+outfile+index+outend, "DATA"])
            blah += 1
            if blah%100 == 0:
                print(blah)

for x in mcdir:
    for file in os.listdir(root+x):
        if file.endswith(".root"):
            index = str(blah)
            call(["WeightDistrubtions", os.path.join(root+x,file), outputdir+outfile+index+outend, "MC"])
            blah += 1
            if blah%100 == 0:
                print(blah)
