import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
for i in range(1,500):
    readFiles.extend( [
        '/store/user/dwalter/TT_Dilepton_v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/tt/180115_123529/0000/output_0_'+str(i)+'.root'
    ])
#for i in range(500,983):
#    readFiles.extend([
#        '/store/user/dwalter/TT_Dilepton_v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/tt/180115_123529/0000/output_0_'+str(i)+'.root'
#    ])
#for i in range(1,500):
#    readFiles.extend( [
#        '/store/user/dwalter/TT_Dilepton_v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/tt_backup/180115_123615/0000/output_0_'+str(i)+'.root'
#    ])
#for i in range(500,1000):
#    readFiles.extend( [
#        '/store/user/dwalter/TT_Dilepton_v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/tt_backup/180115_123615/0000/output_0_'+str(i)+'.root'
#    ])
#for i in range(1000,1191):
#    readFiles.extend( [
#        '/store/user/dwalter/TT_Dilepton_v4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/tt_backup/180115_123615/0001/output_0_'+str(i)+'.root'
#    ])

deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                               vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               secVertices = cms.InputTag("slimmedSecondaryVertices"),
                               jets       = cms.InputTag("GoodJets"),
                               jetR       = cms.double(0.4),
                               runFatJet = cms.bool(False),
                               pupInfo = cms.InputTag("slimmedAddPileupInfo"),
                               lheInfo = cms.InputTag("externalLHEProducer"),
                               rhoInfo = cms.InputTag("fixedGridRhoFastjetAll"),
                               SVs  = cms.InputTag("slimmedSecondaryVertices"),
                               LooseSVs = cms.InputTag("inclusiveCandidateSecondaryVertices"),
                               genJetMatchWithNu = cms.InputTag("patGenJetMatchWithNu"),
                               genJetMatchRecluster = cms.InputTag("patGenJetMatchRecluster"),
                               pruned = cms.InputTag("prunedGenParticles"),
                               fatjets = cms.InputTag('slimmedJetsAK8'),
                               muons = cms.InputTag("slimmedMuons"),
                               electrons = cms.InputTag("slimmedElectrons"),
                               jetPtMin     = cms.double(20.0),
                               jetPtMax     = cms.double(1000),
                               jetAbsEtaMin = cms.double(0.0),
                               jetAbsEtaMax = cms.double(2.5),
                               gluonReduction = cms.double(0.0),
                               tagInfoName = cms.string('deepNN'),
                               tagInfoFName = cms.string('pfBoostedDoubleSVAK8'),
                               bDiscriminators = cms.vstring(),
                               qgtagger        = cms.string("QGTagger"),
                               candidates      = cms.InputTag("packedPFCandidates"),


                               useHerwigCompatible=cms.bool(False),
                               isHerwig=cms.bool(False),
                               useOffsets=cms.bool(True),
                               applySelection=cms.bool(True),
                               isData=cms.bool(False),
                               useLHEWeights=cms.bool(False),

                               #for computation of event weights
                               pupDataDir=cms.string(   #Directory of the data pileup distribution root file for pileup reweighting
                                   "/afs/desy.de/user/d/dwalter/CMSSW_8_0_25/src/DeepNTuples/DeepNtuplizer/data/pileupData.root"),

                               pupMCDir=cms.string(     # Directory of the data pileup distribution root file for pileup reweighting
                                   "/afs/desy.de/user/d/dwalter/CMSSW_8_0_25/src/DeepNTuples/DeepNtuplizer/data/pileupMC.root"),


                               crossSection=cms.double(831760),
                               luminosity = cms.double(35.9),
                               efficiency = cms.double(1./(77081156 + 77867738))  #1/((effective) number of events
                               )