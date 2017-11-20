
import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
### parsing job options 
import sys


options = VarParsing.VarParsing()

options.register('inputScript','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"input Script")
options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('maxEvents',-1,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "total jobs")
options.register('isData', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch off generator jets")

import os
release=os.environ['CMSSW_VERSION'][6:11]
print("Using release "+release)

options.register(
	'inputFiles','',
	VarParsing.VarParsing.multiplicity.list,
	VarParsing.VarParsing.varType.string,
	"input files (default is the tt RelVal)"
	)

if hasattr(sys, "argv"):
    options.parseArguments()


process = cms.Process("ttbarSelector")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if options.isData == True:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'     # For Data Jet Energy correction
if options.isData == False:
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'               # For MC Jet Energy correction

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
if options.inputScript == '': #this is probably for testing
    process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),
   wantSummary=cms.untracked.bool(True)
)

#process.source.fileNames=['file:/afs/cern.ch/work/d/dwalter/data/ttbar/TT/output_0_1.root']   #store/data/Run2016H/SingleMuon/MINIAOD/18Apr2017-v1/00000/00E02A09-853C-E711-93FF-3417EBE644A7.root
#process.source.fileNames=['file:./000C6E52-8BEC-E611-B3FF-0025905C42FE.root']  #

sampleListFile = 'DeepNTuples.DeepNtuplizer.samples.singleMuon_2016_cfg'
process.load(sampleListFile) #default input

if options.inputFiles:
    process.source.fileNames = options.inputFiles

if options.inputScript != '' and options.inputScript != sampleListFile:
    process.load(options.inputScript)

process.source.fileNames=['file:/afs/cern.ch/work/d/dwalter/data/ttbar/TT/output_0_1.root']   #store/data/Run2016H/SingleMuon/MINIAOD/18Apr2017-v1/00000/00E02A09-853C-E711-93FF-3417EBE644A7.root


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

outFileName = options.outputFile + '_' + str(options.job) +  '.root'
print ('Using output file ' + outFileName)

process.MINIAODSIMEventContent.outputCommands.extend([
    'drop *',
    'keep *_NumInteractions_*_*'
])

process.outmod = cms.OutputModule("PoolOutputModule",
                                process.MINIAODSIMEventContent,
                                SelectEvents = cms.untracked.PSet(
                                    SelectEvents = cms.vstring('PUPpath')
                                ),
                                dropMetaData = cms.untracked.string("DROPPED"),
                                fileName = cms.untracked.string(outFileName)
                                )





process.NumInteractions = cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("slimmedElectrons"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("z"),
    eventInfo=cms.untracked.bool(False),
    variables=cms.VPSet(
        cms.PSet(
            tag=cms.untracked.string("Pt"),
            quantity=cms.untracked.string("pt")
        )

    )
)




if options.isData:
    process.PUPpath = cms.Path(process.NumInteractions)
else:
    process.PUPpath = cms.Path(process.NumInteractions)

process.endp = cms.EndPath(process.outmod)
