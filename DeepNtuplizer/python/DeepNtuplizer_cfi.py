import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                               vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               jets       = cms.InputTag("slimmedJets"),
                               genJetMatchWithNu = cms.InputTag("patGenJetMatchWithNu"),
                               genJetMatchRecluster = cms.InputTag("patGenJetMatchRecluster"),
                               jetPtMin     = cms.double(0.0),
                               jetPtMax     = cms.double(10000),
                               jetAbsEtaMin = cms.double(0.0),
                               jetAbsEtaMax = cms.double(5),
                               gluonReduction = cms.double(0.0),

                               )
