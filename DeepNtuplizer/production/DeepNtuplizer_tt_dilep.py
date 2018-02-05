##########
#   tested with CMSSW 8.0.25
#   extract jets from tt_dilep selected root file
#   the input MiniAOD files are specified in an extra config file,
#   for example "DeepNtuplizer_tt_dilep.py" which can be initialized with the option config=tt
##########

import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
### parsing job options
import sys

options = VarParsing.VarParsing()

options.register('inputScript', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"input Script")
options.register('config','tt',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"config for the kind of process under investigation with parameters for deepntuplizer and the input file list")
options.register('outputFile', 'output', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "output File (w/o .root)")
options.register('maxEvents', -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"total jobs")
options.register('gluonReduction', 0.0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.float, "gluon reduction")
options.register('selectJets', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool,"select jets with good gen match")
options.register('isData', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool,"switch off generator jets")



import os

release = os.environ['CMSSW_VERSION'][6:11]
print("Using release " + release)

options.register(
    'inputFiles', '',
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.string,
    "input files (default is the tt RelVal)"
)

if hasattr(sys, "argv"):
    options.parseArguments()

process = cms.Process("ntuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if options.isData:
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'  # For Data Jet Energy correction
else:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'  # For MC Jet Energy correction

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
if options.inputScript == '':  # this is probably for testing
    process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False)
)

if options.config == '':  # default input
    sampleListFile = "DeepNTuples.DeepNtuplizer.samples.samples_template"
    process.load(sampleListFile)
    process.load("DeepNTuples.DeepNtuplizer.DeepNtuplizer_cfi")

else:
    configFile = "DeepNTuples.DeepNtuplizer.DeepNtuplizer_"+options.config+"_cfi"
    print(configFile)
    process.load(configFile)


if options.inputFiles:
    process.source.fileNames = options.inputFiles

if options.inputScript != '' and options.inputScript != sampleListFile:
    process.load(options.inputScript)

#process.source.fileNames = ['file:/afs/cern.ch/work/d/dwalter/data/ttbar/data.root']

numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
if options.nJobs > 1:
    print ("running over these files:")
    print (process.source.fileNames)

process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents)
)

if int(release.replace("_", "")) >= 840:
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

if int(release.replace("_", "")) >= 840:
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

if options.isData:
    jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
else:
    jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], 'None')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
    process,
    labelName="DeepFlavour",
    #         jetSource=cms.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked', 'SubJets'),  # 'subjets from AK8'
    jetSource=cms.InputTag('GoodJets'),  # 'ak4Jets'
    jetCorrections=jetCorrectionsAK4,
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

# from RecoJets.JetProducers.ak4GenJets_cfi import ak4TrackJets


process.ak4GenJetsWithNu = ak4GenJets.clone(src='packedGenParticles')
# Process.ak4GenJetsWithNu = ak4TrackJets.clone(src = 'packedGenParticles')

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
                                                     src=cms.InputTag("packedGenParticles"),
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

outFileName = options.config + '_2_' + str(options.job) + '.root'
print ('Using output file ' + outFileName)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(outFileName))

# DeepNtuplizer
process.deepntuplizer.jets = cms.InputTag('selectedUpdatedPatJetsDeepFlavour')
process.deepntuplizer.bDiscriminators = bTagDiscriminators
process.deepntuplizer.bDiscriminators.append('pfCombinedMVAV2BJetTags')
process.deepntuplizer.LooseSVs = cms.InputTag("looseIVFinclusiveCandidateSecondaryVertices")
process.deepntuplizer.applySelection = cms.bool(options.selectJets)


if int(release.replace("_", "")) >= 840:
    process.deepntuplizer.tagInfoName = cms.string('pfDeepCSV')

process.deepntuplizer.gluonReduction = cms.double(options.gluonReduction)

if options.isData:
    process.p = cms.Path(process.QGTagger + process.deepntuplizer)
else:
    process.p = cms.Path(process.QGTagger + process.genJetSequence * process.deepntuplizer)