
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



sampleListFile = 'DeepNTuples.DeepNtuplizer.samples.singleMuon_2016_cfg'
process.load(sampleListFile) #default input

if options.inputFiles:
    process.source.fileNames = options.inputFiles

if options.inputScript != '' and options.inputScript != sampleListFile:
    process.load(options.inputScript)

#process.source.fileNames=['file:./00E02A09-853C-E711-93FF-3417EBE644A7.root']   #store/data/Run2016H/SingleMuon/MINIAOD/18Apr2017-v1/00000/00E02A09-853C-E711-93FF-3417EBE644A7.root
#process.source.fileNames=['file:./000C6E52-8BEC-E611-B3FF-0025905C42FE.root']  #

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
    'keep *_GoodElectron_*_*',
    'keep *_GoodMuon_*_*',
    'keep *_GoodJets_*_*',
    'keep *_GoodOFLeptonPair_*_*',
    'keep *_GoodFinalSel_*_*'
])

process.outmod = cms.OutputModule("PoolOutputModule",
                                process.MINIAODSIMEventContent,
                                SelectEvents = cms.untracked.PSet(
                                    SelectEvents = cms.vstring('ttbaremupath')
                                ),
                                fileName = cms.untracked.string(outFileName)
                                )

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = [
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff'
]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#HighLevelTrigger
HLTlistSM = cms.vstring("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                        "HLT_Ele27_WPTight_Gsf_v*",
                        "HLT_IsoTkMu24_v*",
                        "HLT_IsoMu24_v*"
                        )
process.hltHighLevelSM = cms.EDFilter("HLTHighLevel",
                                       TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                       HLTPaths = cms.vstring(HLTlistSM),
                                       eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                       andOr = cms.bool(True),   # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                       throw = cms.bool(True)    # throw exception on unknown path names
    )


# Electron Selection
process.EleIdEmbed = cms.EDProducer("ElectronIdAdder",
                                    src=cms.InputTag("slimmedElectrons"),
                                    vSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    idMap=cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")
                                   )

process.GoodElectron = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("EleIdEmbed"),
    cut = cms.string("pt > 20.0 && abs(eta)<2.4 "
                     "&& userFloat('tightcutbased')"
                     "&& userFloat('notInEtaVetoRegion')"
                     "&& userFloat('inAbsD0')"
                     "&& userFloat('inAbsDz')"
                     )
    )

### Muon Selection
process.MuonIdEmbed = cms.EDProducer("MuonIdAdder",
                                     src=cms.InputTag("slimmedMuons"),
                                     vSrc=cms.InputTag("offlineSlimmedPrimaryVertices")
                                     )

process.GoodMuon = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("MuonIdEmbed"),
    cut = cms.string("pt > 20.0 && abs(eta)<2.4 "
                     "&& userFloat('tightcutbased')"
                     "&& (pfIsolationR04().sumChargedHadronPt + max(0.0, (pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt -0.5*pfIsolationR04().sumPUPt)))/pt < 0.15"
    )
)
### Jets

# Jet Energy Corrections
if options.isData == True:
    corrections = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
if options.isData == False:
    corrections = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', corrections, 'None')
)

process.GoodJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("updatedPatJetsUpdatedJEC"),
    preselection = cms.string("pt>30 && abs(eta) < 2.4 && neutralHadronEnergyFraction < 0.99 "
                     "&& neutralEmEnergyFraction < 0.99 && (chargedMultiplicity+neutralMultiplicity) > 1 "
                     "&& chargedHadronEnergyFraction > 0.0 "
                     "&& chargedMultiplicity > 0.0 && chargedEmEnergyFraction < 0.99 "),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("GoodMuon"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.4),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps cause the jet to be discarded
        ),
        electrons = cms.PSet(
           src       = cms.InputTag("GoodElectron"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.4),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps cause the jet to be discarded
        )
    ),
    finalCut = cms.string('')
  )

process.GoodOFLeptonPair = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("GoodMuon@+ GoodElectron@-"),
    cut = cms.string("20.0 < mass"
                     "&& (daughter(0).pt > 25 || daughter(1).pt > 25)")
  )


process.FinalSel = cms.EDFilter("CandViewCountFilter",
     src = cms.InputTag("GoodOFLeptonPair"),
     minNumber = cms.uint32(1),
  )



if options.isData:
    process.ttbaremupath = cms.Path(process.hltHighLevelSM + process.FinalSel)
else:
    process.ttbaremupath = cms.Path(process.hltHighLevelSM + process.FinalSel)

process.endp = cms.EndPath(process.outmod)
