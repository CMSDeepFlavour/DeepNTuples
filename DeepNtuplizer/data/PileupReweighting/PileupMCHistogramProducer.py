

import pdb
import ROOT
#from DataFormats.FWLite import Events, Handle
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()
options.register('inputScript','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"input Script")
options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('maxEvents',-1,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "total jobs")
options.register(
	'inputFiles','',
	VarParsing.VarParsing.multiplicity.list,
	VarParsing.VarParsing.varType.string,
	"input files (default is the tt RelVal)"
	)
import os
release=os.environ['CMSSW_VERSION'][6:11]
print("Using release "+release)

# Set up a process, named RECO in this case
process = cms.Process("PileupMCHistProd")

#default input
sampleListFile = 'DeepNTuples.DeepNtuplizer.samples.singleMuon_2016_cfg'
process.load(sampleListFile)
process.source.fileNames=['file:/afs/cern.ch/work/d/dwalter/data/ttbar/TT/output_0_1.root']

if options.inputFiles:
    process.source.fileNames = options.inputFiles

if options.inputScript != '' and options.inputScript != sampleListFile:
    process.load(options.inputScript)

numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
if options.nJobs > 1:
    print ("running over these files:")
    print (process.source.fileNames)

process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents  = cms.untracked.PSet(
    input = cms.untracked.int32 (options.maxEvents)
)


process.MINIAODSIMEventContent.outputCommands.extend([
    'keep *_NumInteractions_*_*',

])


process.NumInteractions = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("slimmedAddPileupInfo"),
    cut = cms.string("getBunchCrossing() == 0")
    )


# #Dont know if this is necessary
# ROOT.gSystem.Load("libFWCoreFWLite.so")
# ROOT.AutoLibraryLoader.enable()
# ROOT.gSystem.Load("libDataFormatsFWLite.so")
# ROOT.gSystem.Load("libDataFormatsPatCandidates.so")
#
# infile = "ttbar/TT/output_0_1.root"
# directory = ""
# islist = ("root" != infile[-len(infile.split('.')[-1]):])
#
# outfile = ROOT.TFile("MyMCPileupHistogram.root", "RECREATE")
#
# hist_nInt = ROOT.TH1F("numInteractions", "Number of Interactions per bunch corssing", 50, 0, 50)
#
# print("make handle")
# handle = Handle('vector<PileupSummaryInfo>')
# label = ("slimmedAddPileupInfo")
#
# if islist:
#     print("loop over filelist " + str(infile))
#     filelist = open(infile, 'r')
#     ifile = filelist.readline()
#     ifile = ifile[:-1]
#     directory = infile[:-len(infile.split('/')[-1])]
# else:
#     ifile = infile
#
# while ifile != '':
#     print("process file ",ifile)
#     events = Events(directory + ifile)
#     for event in events:
#         event.getByLabel(label, handle)
#         pupInfos = handle.product()
#
#         print(pupInfos[12].getBunchCrossing(),"   ", pupInfos[12].getPU_NumInteractions())
#         hist_nInt.Fill(pupInfos[12].getPU_NumInteractions())
#
#     if islist:
#         ifile = filelist.readline()
#         ifile = ifile[:-1]
#     else:
#         break
#
# hist_nInt.Write()


outFileName = options.outputFile + '_' + str(options.job) +  '.root'

# Configure the object that writes an output file
process.out = cms.OutputModule("PoolOutputModule",
    process.MINIAODSIMEventContent,
    SelectEvents=cms.untracked.PSet(
        SelectEvents=cms.vstring('mcHistProd')
    ),
    fileName=cms.untracked.string(outFileName),
    )
process.mcHistProd = cms.Path(process.NumInteractions)
# Configure a path and endpath to run the producer and output modules
process.ep = cms.EndPath(process.out)





