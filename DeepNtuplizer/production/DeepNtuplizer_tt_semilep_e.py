# this can be run with CMSSW 8_0_29; in CMSSW 8_0_25 the module 'cutBasedElectronID_Summer16_80X_V1_cff' is missing
#basically deepntuplizer with tt semileptonic single electron selection

import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
### parsing job options
import sys

options = VarParsing.VarParsing()

options.register('inputScript', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"input Script")
options.register('outputFile', 'output', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "output File (w/o .root)")
options.register('maxEvents', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "total jobs")
options.register('gluonReduction', 0.0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "gluon reduction")
options.register('selectJets', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool,"select jets with good gen match")
options.register('globalTag', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"global tag for jet energy correction")
options.register('isData', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "switch off generator jets")
options.register('deepNtuplizer',True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "run deepNtuplizer or just the ttbar selection")
options.register('lheWeights',False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use LHE weights")
options.register('crossSection', 1., VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float,"cross section")
options.register('nEvents', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"number of events to process")


import os

release = os.environ['CMSSW_VERSION'][6:12]
print("Using release " + release)

options.register(
    'inputFiles', '',
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.string,
    "input files (default is the tt RelVal)"
)

if hasattr(sys, "argv"):
    options.parseArguments()

process = cms.Process("semilepElectronSelectedDNNFiller")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if options.globalTag == '':
    if options.isData == False:
        process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'     # For MC Jet Energy correction
    if options.isData == True:
        process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' # For Data Jet Energy correction
else:
    process.GlobalTag.globaltag = options.globalTag

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))
#process.PoolSource = cms.untracked.PSet(
#    firstEvent = cms.untracked.uint32(33285817),
#)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
if options.inputScript == '':  # this is probably for testing
    process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.load('DeepNTuples.DeepNtuplizer.samples.TTJetsPhase1_cfg')  # default input

if options.inputFiles:
    process.source.fileNames = options.inputFiles

if options.inputScript != '' and options.inputScript != 'DeepNTuples.DeepNtuplizer.samples.TTJetsPhase1_cfg':
    process.load(options.inputScript)


#process.source.fileNames=['file:./00CC509E-0C3B-E711-98F2-0242AC130004.root']   #isData=True
#process.source.fileNames=['file:./0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root']    #isData=False
#process.source.fileNames=['file:./EE95DEDC-96BE-E611-B45D-A0000420FE80.root']    #isData=False
#process.source.fileNames=['file:./0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root']    #store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/
#process.source.fileNames=['file:./00E02A09-853C-E711-93FF-3417EBE644A7.root']    #store/data/Run2016H/SingleMuon/MINIAOD/18Apr2017-v1/00000/
#process.source.fileNames=['file:./F21AE451-7EBA-E611-9399-0025905B858E.root']    #store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/

#process.source.fileNames=['file:./00CC509E-0C3B-E711-98F2-0242AC130004.root']   #Run2016B/SingleMuon/MINIAOD/18Apr2017_ver2-v1/
numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
#if options.nJobs > 1:
print ("running over these files:")
print (process.source.fileNames)

process.globalInfo = cms.EDAnalyzer('globalInfo',
                                lheInfo = cms.InputTag("externalLHEProducer"),
                                useLHEWeights=cms.bool(options.lheWeights)
                               )

process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents)
)

if int(release.replace("_", "")) >= 8400 or int(release.replace("_", "")) == 8029:
    bTagInfos = [
        'pfImpactParameterTagInfos',
        'pfInclusiveSecondaryVertexFinderTagInfos',
        'pfDeepCSVTagInfos']
else:
    bTagInfos = [
        'pfImpactParameterTagInfos',
        'pfInclusiveSecondaryVertexFinderTagInfos',
        'deepNNTagInfos',
    ]

if int(release.replace("_", "")) >= 8400 or int(release.replace("_", "")) == 8029:
    bTagDiscriminators = [
        'softPFMuonBJetTags',
        'softPFElectronBJetTags',
        'pfJetBProbabilityBJetTags',
        'pfJetProbabilityBJetTags',
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probudsg',  # to be fixed with new names
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probbb',
        'pfDeepCSVJetTags:probcc',
    ]
else:
    bTagDiscriminators = [
        'softPFMuonBJetTags',
        'softPFElectronBJetTags',
        'pfJetBProbabilityBJetTags',
        'pfJetProbabilityBJetTags',
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'deepFlavourJetTags:probudsg',  # to be fixed with new names
        'deepFlavourJetTags:probb',
        'deepFlavourJetTags:probc',
        'deepFlavourJetTags:probbb',
        'deepFlavourJetTags:probcc',
    ]

###### semilep selection
outFileName = options.outputFile + '_' + str(options.job) + '.root'
print ('Using output file ' + outFileName)

if options.deepNtuplizer:
    process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(outFileName))
else:
    process.MINIAODSIMEventContent.outputCommands.extend([
        'keep *_goodElectrons_*_*',
        'keep *_goodMuons_*_*',
        'keep *_goodJets_*_*',
    ])

    process.outmod = cms.OutputModule("PoolOutputModule",
                                    process.MINIAODSIMEventContent,
                                    SelectEvents = cms.untracked.PSet(
                                        SelectEvents = cms.vstring('p')
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
HLTlistSM = cms.vstring("HLT_Ele27_WPTight_Gsf_v*"
                        )
process.TriggerSel = cms.EDFilter("HLTHighLevel",
                                       TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                       HLTPaths = cms.vstring(HLTlistSM),
                                       eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                       andOr = cms.bool(True),   # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                       throw = cms.bool(True)    # throw exception on unknown path names
)
### Electron Selection
process.goodElectrons = cms.EDProducer("ElectronIdAdder",
                                    src=cms.InputTag("slimmedElectrons"),
                                    vSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    idMap=cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
                                    minPt=cms.double(20.0),
                                    maxAbsEta=cms.double(2.4),
                                    )

### Muon Selection
process.goodMuons = cms.EDProducer("MuonIdAdder",
                                     src=cms.InputTag("slimmedMuons"),
                                     vSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                     minPt=cms.double(30.0),
                                     maxAbsEta=cms.double(2.4),
                                     maxRMI=cms.double(0.15)    #RMI = relative muon isolation
                                     )

### Jet Selection
process.goodJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("slimmedJets"),
    preselection = cms.string("pt>30 && abs(eta) < 2.4 && neutralHadronEnergyFraction < 0.99 "
                     "&& neutralEmEnergyFraction < 0.99 && (chargedMultiplicity+neutralMultiplicity) > 1 "
                     "&& chargedHadronEnergyFraction > 0.0 "
                     "&& chargedMultiplicity > 0.0 && chargedEmEnergyFraction < 0.99 "),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("goodMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.4),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps cause the jet to be discarded
        ),
        electrons = cms.PSet(
           src       = cms.InputTag("goodElectrons"),
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


### Event Filter
process.ttsemilepFilter = cms.EDFilter("ttsemilepFilter",
                                       src_muons = cms.InputTag("goodMuons"),
                                       src_electrons = cms.InputTag("goodElectrons"),
                                       src_jets = cms.InputTag("goodJets"),
                                       src_mets = cms.InputTag("slimmedMETs"),
                                       cut_minMT_leptonMETpair = cms.double(50.0),
                                       nElectrons=cms.uint32(1),
                                       nMuons=cms.uint32(0)
                                       )



### end selection

# Jet Energy Corrections
if options.isData:
    jetCorrections = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
else:
    jetCorrections = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
        process,
        labelName="DeepFlavour",
        #         jetSource=cms.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked', 'SubJets'),  # 'subjets from AK8'
        jetSource=cms.InputTag('goodJets'),  # 'ak4Jets'
        jetCorrections=('AK4PFchs', jetCorrections, 'None'),
        pfCandidates=cms.InputTag('packedPFCandidates'),
        pvSource=cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource=cms.InputTag('slimmedSecondaryVertices'),
        muSource=cms.InputTag('slimmedMuons'),
        elSource=cms.InputTag('slimmedElectrons'),
        btagInfos=bTagInfos,
        btagDiscriminators=bTagDiscriminators,
        explicitJTA=False
        )

if hasattr(process, 'updatedPatJetsTransientCorrectedDeepFlavour'):
    process.updatedPatJetsTransientCorrectedDeepFlavour.addTagInfos = cms.bool(True)
    process.updatedPatJetsTransientCorrectedDeepFlavour.addBTagInfo = cms.bool(True)
else:
    raise ValueError(
        'I could not find updatedPatJetsTransientCorrectedDeepFlavour to embed the tagInfos, please check the cfg')

# QGLikelihood
process.load("DeepNTuples.DeepNtuplizer.QGLikelihood_cfi")
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource", "QGPoolDBESSource")
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag("selectedUpdatedPatJetsDeepFlavour")
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets

process.ak4GenJetsWithNu = ak4GenJets.clone(src='packedGenParticles')

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src=cms.InputTag("packedGenParticles"),
                                                     cut=cms.string(
                                                         "abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
## Define GenJets
process.ak4GenJetsRecluster = ak4GenJets.clone(src='packedGenParticlesForJetsNoNu')

process.patGenJetMatchWithNu = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
                                              src=cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),
                                              # RECO jets (any View<Jet> is ok)
                                              matched=cms.InputTag("ak4GenJetsWithNu"),
                                              # GEN jets  (must be GenJetCollection)
                                              mcPdgId=cms.vint32(),  # n/a
                                              mcStatus=cms.vint32(),  # n/a
                                              checkCharge=cms.bool(False),  # n/a
                                              maxDeltaR=cms.double(0.4),  # Minimum deltaR for the match
                                              # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
                                              resolveAmbiguities=cms.bool(True),
                                              # Forbid two RECO objects to match to the same GEN object
                                              resolveByMatchQuality=cms.bool(False),
                                              # False = just match input in order; True = pick lowest deltaR pair first
                                              )

process.patGenJetMatchRecluster = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
                                                 src=cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),
                                                 # RECO jets (any View<Jet> is ok)
                                                 matched=cms.InputTag("ak4GenJetsRecluster"),
                                                 # GEN jets  (must be GenJetCollection)
                                                 mcPdgId=cms.vint32(),  # n/a
                                                 mcStatus=cms.vint32(),  # n/a
                                                 checkCharge=cms.bool(False),  # n/a
                                                 maxDeltaR=cms.double(0.4),  # Minimum deltaR for the match
                                                 # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
                                                 resolveAmbiguities=cms.bool(True),
                                                 # Forbid two RECO objects to match to the same GEN object
                                                 resolveByMatchQuality=cms.bool(False),
                                                 # False = just match input in order; True = pick lowest deltaR pair first
                                                 )

process.genJetSequence = cms.Sequence(process.packedGenParticlesForJetsNoNu
                                      * process.ak4GenJetsWithNu
                                      * process.ak4GenJetsRecluster
                                      * process.patGenJetMatchWithNu
                                      * process.patGenJetMatchRecluster)

# Very Loose IVF SV collection
from PhysicsTools.PatAlgos.tools.helpers import loadWithPrefix

loadWithPrefix(process, 'RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff', "looseIVF")
process.looseIVFinclusiveCandidateVertexFinder.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFinclusiveCandidateVertexFinder.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLen2DSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLenSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.fitterSigmacut = 20

process.looseIVFcandidateVertexArbitrator.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFcandidateVertexArbitrator.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFcandidateVertexArbitrator.secondaryVertices = cms.InputTag("looseIVFcandidateVertexMerger")
process.looseIVFcandidateVertexArbitrator.fitterSigmacut = 20

#datapath=os.environ['CMSSW_BASE']+'/src/DeepNTuples/DeepNtuplizer/data/'
datapath=''

# DeepNtuplizer configurations
process.load("DeepNTuples.DeepNtuplizer.DeepNtuplizer_cfi")
process.deepntuplizer.jets = cms.InputTag('selectedUpdatedPatJetsDeepFlavour')
process.deepntuplizer.bDiscriminators = bTagDiscriminators
process.deepntuplizer.bDiscriminators.append('pfCombinedMVAV2BJetTags')
process.deepntuplizer.LooseSVs = cms.InputTag("looseIVFinclusiveCandidateSecondaryVertices")
process.deepntuplizer.applySelection = cms.bool(options.selectJets)

process.deepntuplizer.isData = cms.bool(options.isData)

process.deepntuplizer.removeUndefined = cms.bool(False)

if not options.isData:

    process.deepntuplizer.useLHEWeights = cms.bool(options.lheWeights)

    process.deepntuplizer.pileupData=cms.string(datapath+"pileup_data_2016.root")
    process.deepntuplizer.pileupMC=cms.string(datapath+"pileup_MC_2016.root")

    process.deepntuplizer.sfElectrons = cms.InputTag("goodElectrons")

    process.deepntuplizer.periods=cms.vstring("2016All",)
    process.deepntuplizer.lumis=cms.vdouble(1.,)
    process.deepntuplizer.crossSection=cms.double(options.crossSection)
    process.deepntuplizer.nEvents=cms.uint32(options.nEvents)

    process.deepntuplizer.sfTrigger_e=cms.vstring(datapath + "TriggerSF_Run2016All_v1.root",)
    process.deepntuplizer.sfTrigger_e_Hist=cms.vstring("Ele27_WPTight_Gsf",)

    process.deepntuplizer.sfElIdAndIso=cms.vstring(datapath+"egammaEffi.txt_EGM2D.root",)
    process.deepntuplizer.sfElIdAndIso_Hist=cms.vstring("EGamma_SF2D",)

    process.deepntuplizer.gluonReduction = cms.double(options.gluonReduction)

if int(release.replace("_", "")) >= 840:
    process.deepntuplizer.tagInfoName = cms.string('pfDeepCSV')


# 1631
#process.ProfilerService = cms.Service(
#    "ProfilerService",
#    firstEvent=cms.untracked.int32(1631),
#    lastEvent=cms.untracked.int32(1641),
#    paths=cms.untracked.vstring('p')
#)


if options.deepNtuplizer:
    if options.isData:
        process.p = cms.Path(process.globalInfo + process.TriggerSel + process.ttsemilepFilter + process.QGTagger + process.deepntuplizer)
    else:
        process.p = cms.Path(process.globalInfo + process.TriggerSel + process.ttsemilepFilter + process.QGTagger + process.genJetSequence * process.deepntuplizer)
else:
        process.p = cms.Path(process.TriggerSel + process.ttsemilepFilter)

        process.endp = cms.EndPath(process.outmod)

