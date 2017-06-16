import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
### parsing job options 
import sys, copy

options = VarParsing.VarParsing()

options.register('inputScript','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"input Script")
options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('maxEvents',-1,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "total jobs")
options.register('gluonReduction', 0.0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "gluon reduction")

options.register(
  'inputFiles','',
  VarParsing.VarParsing.multiplicity.list,
  VarParsing.VarParsing.varType.string,
  "input files (default is the tt RelVal)"
  )

if hasattr(sys, "argv"):
    options.parseArguments()




process = cms.Process("DNNFiller")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
if options.inputScript == '': #this is probably for testing
  process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),  
   wantSummary=cms.untracked.bool(False)
)

from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValTTbarPileUpMINIAODSIM

#process.load('DeepNTuples.DeepNtuplizer.samples.TTJetsPhase1_cfg') #default input
process.load('DeepNTuples.DeepNtuplizer.samples.QCD_Pt_600to800Phase1_cfg') #default input
#process.load('DeepNTuples.DeepNtuplizer.samples.BGHH4b_narrow_M3000_Phase1_cfg') #default input


if options.inputFiles:
  process.source.fileNames = options.inputFiles

if options.inputScript != '' and options.inputScript != 'DeepNTuples.DeepNtuplizer.samples.TTJetsPhase1_cfg':
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


bTagInfos = [
  'pfImpactParameterTagInfos',
  'pfInclusiveSecondaryVertexFinderTagInfos',
  'deepNNTagInfos',
]
bTagDiscriminators = [
  'softPFMuonBJetTags',
  'softPFElectronBJetTags',
  'pfJetBProbabilityBJetTags',
  'pfJetProbabilityBJetTags',
  'pfCombinedInclusiveSecondaryVertexV2BJetTags',
  'deepFlavourJetTags:probudsg', #to be fixed with new names
  'deepFlavourJetTags:probb', 
  'deepFlavourJetTags:probc', 
  'deepFlavourJetTags:probbb', 
  'deepFlavourJetTags:probcc',
]
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
        process,
        labelName = "DeepFlavour",
        jetSource = cms.InputTag('slimmedJets'),#'ak4Jets'
        jetCorrections = jetCorrectionsAK4,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
)

updateJetCollection(
        process,
        labelName = "AK8DeepFlavour",
        jetSource = cms.InputTag('slimmedJetsAK8'),#'ak4Jets'
        jetCorrections = jetCorrectionsAK8,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
)

updateJetCollection(
    process,
    labelName='SoftDropSubjetsDeepFlavour',
    jetSource=cms.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked:SubJets'),
    jetCorrections = jetCorrectionsAK4,
    pfCandidates = cms.InputTag('packedPFCandidates'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    muSource = cms.InputTag('slimmedMuons'),
    elSource = cms.InputTag('slimmedElectrons'),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    explicitJTA = True,   # needed for subjet b tagging
    svClustering = False, # needed for subjet b tagging (IMPORTANT: Needs to be set to False to disable ghost-association which does not work with slimmed jets)
    fatJets = cms.InputTag('slimmedJetsAK8'), # needed for subjet b tagging
    rParam=0.8, # needed for subjet b tagging
    algo='AK',              # has to be defined but is not used since svClustering=False
)

if hasattr(process,'updatedPatJetsTransientCorrectedDeepFlavour'):
  process.updatedPatJetsTransientCorrectedDeepFlavour.addTagInfos = cms.bool(True) 
  process.updatedPatJetsTransientCorrectedDeepFlavour.addBTagInfo = cms.bool(True)
else:
  raise ValueError('I could not find updatedPatJetsTransientCorrectedDeepFlavour to embed the tagInfos, please check the cfg')

if hasattr(process,'updatedPatJetsTransientCorrectedSoftDropSubjetsDeepFlavour'):
  process.updatedPatJetsTransientCorrectedSoftDropSubjetsDeepFlavour.addTagInfos = cms.bool(True) 
  process.updatedPatJetsTransientCorrectedSoftDropSubjetsDeepFlavour.addBTagInfo = cms.bool(True)
  print ">>>Adding ValueMaps for updatedPatJetsTransientCorrectedSoftDropSubjetsDeepFlavour"
else:
  raise ValueError('I could not find updatedPatJetsTransientCorrectedSoftDropSubjetsDeepFlavour to embed the tagInfos, please check the cfg')

if hasattr(process,'updatedPatJetsTransientCorrectedAK8DeepFlavour'):
  process.updatedPatJetsTransientCorrectedAK8DeepFlavour.addTagInfos = cms.bool(True) 
  process.updatedPatJetsTransientCorrectedAK8DeepFlavour.addBTagInfo = cms.bool(True)
else:
  raise ValueError('I could not find updatedPatJetsTransientCorrectedAK8DeepFlavour to embed the tagInfos, please check the cfg')

# QGLikelihood
process.load("DeepNTuples.DeepNtuplizer.QGLikelihood_cfi")
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource", "QGPoolDBESSource")
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag("selectedUpdatedPatJetsDeepFlavour")
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

process.QGTaggerSubJet = process.QGTagger.clone(
    srcJets = cms.InputTag("selectedUpdatedPatJetsSoftDropSubjetsDeepFlavour"), 
    jetsLabel = cms.string('QGL_AK4PFchs'),
    )

process.QGTaggerAK8 = process.QGTagger.clone(
    srcJets = cms.InputTag("selectedUpdatedPatJetsAK8DeepFlavour"), 
    jetsLabel = cms.string('QGL_AK4PFchs'),
    )


from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsWithNu = ak4GenJets.clone(src = 'packedGenParticles')
 
 ## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))

 ## Define GenJets
process.ak4GenJetsRecluster = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')
 
process.patGenJetMatchWithNu = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR           
    src         = cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),      # RECO jets (any View<Jet> is ok) 
    matched     = cms.InputTag("ak4GenJetsWithNu"),        # GEN jets  (must be GenJetCollection)              
    mcPdgId     = cms.vint32(),                      # n/a   
    mcStatus    = cms.vint32(),                      # n/a   
    checkCharge = cms.bool(False),                   # n/a   
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match   
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)                     
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object 
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first          
)

process.patGenJetMatchRecluster = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR           
    src         = cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),      # RECO jets (any View<Jet> is ok) 
    matched     = cms.InputTag("ak4GenJetsRecluster"),        # GEN jets  (must be GenJetCollection)              
    mcPdgId     = cms.vint32(),                      # n/a   
    mcStatus    = cms.vint32(),                      # n/a   
    checkCharge = cms.bool(False),                   # n/a   
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match   
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)                     
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object 
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first          
)

process.ak8GenJetsRecluster = ak4GenJets.clone(
    src = 'packedGenParticlesForJetsNoNu',
    jetAlgorithm = 'AntiKt',
    rParam = cms.double(0.8), 
    )

process.ak8GenJetsWithNu = ak4GenJets.clone(
    src = 'packedGenParticles',
    jetAlgorithm = 'AntiKt',
    rParam = cms.double(0.8),
    )

process.patGenJetsAK8MatchRecluster = cms.EDProducer("GenJetMatcher",
    src                   = cms.InputTag("selectedUpdatedPatJetsAK8DeepFlavour"),
    matched               = cms.InputTag("ak8GenJetsRecluster"),
    mcPdgId               = cms.vint32(),
    mcStatus              = cms.vint32(),
    checkCharge           = cms.bool(False),
    maxDeltaR             = cms.double(0.8),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    )

process.patGenJetsAK8MatchWithNu = cms.EDProducer("GenJetMatcher",  
    src                   = cms.InputTag("selectedUpdatedPatJetsAK8DeepFlavour"),      
    matched               = cms.InputTag("ak8GenJetsWithNu"),
    mcPdgId               = cms.vint32(), 
    mcStatus              = cms.vint32(), 
    checkCharge           = cms.bool(False), 
    maxDeltaR             = cms.double(0.8), 
    resolveAmbiguities    = cms.bool(True), 
    resolveByMatchQuality = cms.bool(False), 
)

process.genFatJetsSoftDropWithNu = ak4GenJets.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticles"), 
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0 = cms.double(0.8),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

process.genFatJetsSoftDropRecluster = ak4GenJets.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu"), 
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0 = cms.double(0.8),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

process.patGenSubJetsMatchWithNu = cms.EDProducer("GenJetMatcher",  
    src                   = cms.InputTag("selectedUpdatedPatJetsSoftDropSubjetsDeepFlavour"),      
    matched               = cms.InputTag("genFatJetsSoftDropWithNu:SubJets"),
    mcPdgId               = cms.vint32(), 
    mcStatus              = cms.vint32(), 
    checkCharge           = cms.bool(False), 
    maxDeltaR             = cms.double(0.4), 
    resolveAmbiguities    = cms.bool(True), 
    resolveByMatchQuality = cms.bool(False), 
)

process.patGenSubJetsMatchRecluster = cms.EDProducer("GenJetMatcher",  
    src                   = cms.InputTag("selectedUpdatedPatJetsSoftDropSubjetsDeepFlavour"),      
    matched               = cms.InputTag("genFatJetsSoftDropRecluster:SubJets"),
    mcPdgId               = cms.vint32(), 
    mcStatus              = cms.vint32(), 
    checkCharge           = cms.bool(False), 
    maxDeltaR             = cms.double(0.4), 
    resolveAmbiguities    = cms.bool(True), 
    resolveByMatchQuality = cms.bool(False), 
)

process.genJetSequence = cms.Sequence(process.packedGenParticlesForJetsNoNu
    #* process.ak4GenJetsWithNu 
    #* process.ak4GenJetsRecluster
    #* process.patGenJetMatchWithNu 
    #* process.patGenJetMatchRecluster
    * process.ak8GenJetsRecluster
    * process.ak8GenJetsWithNu
    * process.genFatJetsSoftDropWithNu
    * process.genFatJetsSoftDropRecluster
    * process.patGenSubJetsMatchWithNu
    * process.patGenSubJetsMatchRecluster
    * process.patGenJetsAK8MatchWithNu
    * process.patGenJetsAK8MatchRecluster
    )

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


outFileName = options.outputFile + '_' + str(options.job) +  '.root'
print ('Using output file ' + outFileName)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(outFileName))

# DeepNtuplizer
process.load("DeepNTuples.DeepNtuplizer.DeepNtuplizer_cfi")
process.deepntuplizer.jets = cms.InputTag('selectedUpdatedPatJetsDeepFlavour');
process.deepntuplizer.bDiscriminators = bTagDiscriminators 
process.deepntuplizer.bDiscriminators.append('pfCombinedMVAV2BJetTags')
process.deepntuplizer.LooseSVs = cms.InputTag("looseIVFinclusiveCandidateSecondaryVertices")
process.deepntuplizer.gluonReduction  = cms.double(options.gluonReduction)

process.deepntuplizerSj = process.deepntuplizer.clone(
    jets = cms.InputTag('selectedUpdatedPatJetsSoftDropSubjetsDeepFlavour'),
    jetPtMin = cms.double(20.0),
    genJetMatchWithNu = cms.InputTag('patGenSubJetsMatchWithNu'), 
    genJetMatchRecluster = cms.InputTag('patGenSubJetsMatchRecluster'),
    qgtagger = cms.string("QGTaggerSubJet"),
    )

process.deepntuplizerAK8 = process.deepntuplizer.clone(
    jets = cms.InputTag('selectedUpdatedPatJetsAK8DeepFlavour'),
    jetPtMin = cms.double(200.0),
    genJetMatchWithNu = cms.InputTag('patGenJetsAK8MatchWithNu'),
    genJetMatchRecluster = cms.InputTag('patGenJetsAK8MatchRecluster'),
    qgtagger = cms.string("QGTaggerAK8"),
    )
process.deepntuplizerAK8.bDiscriminators.append('pfBoostedDoubleSecondaryVertexAK8BJetTags')

process.p = cms.Path(
    #  process.QGTagger 
    #* 
    process.QGTaggerSubJet
    * process.QGTaggerAK8 
    * process.genJetSequence 
    #* process.deepntuplizer
    * process.deepntuplizerSj
    * process.deepntuplizerAK8
    ) 

open('dump.py', 'w').write(process.dumpPython())
