import FWCore.ParameterSet.Config as cms
import os

datapath=os.environ['CMSSW_BASE']+'/src/DeepNTuples/DeepNtuplizer/data/'

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
                                pileupData=cms.string(datapath+"pileup_data_2016.root"),
                                pileupMC=cms.string(datapath+"pileup_MC_2016.root"),

                                #scalefactor information
                                sfMuons = cms.InputTag("goodMuons"),
                                    # leave an empty string if you don't want to use a scalefactor
                                sfMuonId = cms.string(datapath+"EfficienciesAndSF_ID_GH.root"),
                                sfMuonIdName = cms.string("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")

                                )
