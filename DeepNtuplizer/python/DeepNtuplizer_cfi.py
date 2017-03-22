import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                               vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               secVertices = cms.InputTag("slimmedSecondaryVertices"),
                               jets       = cms.InputTag("slimmedJets"),
                               SVs  = cms.InputTag("slimmedSecondaryVertices"),
                               genJetMatchWithNu = cms.InputTag("patGenJetMatchWithNu"),
                               genJetMatchRecluster = cms.InputTag("patGenJetMatchRecluster"),
                               pruned = cms.InputTag("prunedGenParticles"),
                               jetPtMin     = cms.double(10.0),
                               jetPtMax     = cms.double(1000),
                               jetAbsEtaMin = cms.double(0.0),
                               jetAbsEtaMax = cms.double(2.4),
                               gluonReduction = cms.double(0.0),
															 tagInfoName = cms.string('deepNN'),
															 bDiscriminators = cms.vstring(),
                               )
