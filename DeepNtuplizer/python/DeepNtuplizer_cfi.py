import FWCore.ParameterSet.Config as cms
import os


deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                                vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                secVertices = cms.InputTag("slimmedSecondaryVertices"),
                                jets       = cms.InputTag("slimmedJets"),
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
                                useLHEWeights=cms.bool(True),

                                #leave empty string if you don't want to use pileup weights
                                pileupData=cms.string(""),
                                pileupMC=cms.string(""),

                                #scalefactor information
                                sfMuons = cms.InputTag("slimmedMuons"),
                                sfElectrons=cms.InputTag("slimmedElectrons"),

                                # leave an empty string for the root file if you don't want to use a scalefactor
                                sfMuonTrigger=cms.string(""),
                                sfMuonTriggerHist = cms.string(""),
                                sfMuonId = cms.string(""),
                                sfMuonIdHist = cms.string(""),
                                sfMuonIso=cms.string(""),
                                sfMuonIsoHist=cms.string(""),
                                sfMuonTracking=cms.string(""),
                                sfMuonTrackingHist=cms.string(""),
                                sfElIdAndIso=cms.string(""),
                                sfElIdAndIsoHist=cms.string(""),



                                )
