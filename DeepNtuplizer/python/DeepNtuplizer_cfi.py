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

                                removeUndefined=cms.bool(True),
                                isData=cms.bool(False),

                                #All following is for event weight computation
                                useLHEWeights=cms.bool(True),

                                #leave empty string if you don't want to use pileup weights
                                pileupData=cms.string(""),
                                pileupMC=cms.string(""),

                                periods=cms.vstring(),
                                lumis=cms.vdouble(),        #weights are weighted with the luminosity of the period

                                #to use different triggers for the periods
                                triggerToken=cms.InputTag("TriggerResults::HLT"),
                                triggers=cms.vstring(),

                                #scalefactor information
                                sfMuons = cms.InputTag("slimmedMuons"),
                                sfElectrons=cms.InputTag("slimmedElectrons"),

                                # The scalefactors can be computed on different ways:
                                # - empty vector to not use a specific scalefactor,
                                # - one string to use the same hist for all periods
                                # - number of strings equal to number of periods to use a different hist for each period
                                sfTrigger_mu=cms.vstring(),
                                sfTrigger_mu_Hist = cms.vstring(),

                                sfTrigger_emu=cms.vstring(),
                                sfTrigger_emu_Hist=cms.vstring(),
                                sfMuonId = cms.vstring(),
                                sfMuonId_Hist = cms.vstring(),
                                sfMuonIso=cms.vstring(),
                                sfMuonIso_Hist=cms.vstring(),
                                sfMuonTracking=cms.vstring(),
                                sfMuonTracking_Hist=cms.vstring(),
                                sfElIdAndIso=cms.vstring(),
                                sfElIdAndIso_Hist=cms.vstring(),



                                )
