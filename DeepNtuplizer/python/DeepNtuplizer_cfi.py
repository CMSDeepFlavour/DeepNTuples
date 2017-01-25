import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                               vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               jets       = cms.InputTag("slimmedJets")
                               )
