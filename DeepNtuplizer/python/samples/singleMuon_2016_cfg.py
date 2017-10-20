import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/000C6E52-8BEC-E611-B3FF-0025905C42FE.root',   #-> singleMuon_0.root
'/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/001E3E7D-57EB-E611-8469-0CC47A7C35D2.root',   #-> singleMuon_1.root
'/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/002BBC83-57EB-E611-9701-0CC47A4D7692.root',
'/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/00359BA2-35EC-E611-AD92-0025905C3E68.root'
#'/store/mc/PhaseIFall16MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PhaseIFall16PUFlat20to50_PhaseIFall16_81X_upgrade2017_realistic_v26-v1/50000/08358A47-61E3-E611-8B77-001E677928AE.root'   #-> TTJetsPhase1_0.root
#'/store/mc/PhaseIFall16MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PhaseIFall16PUFlat20to50_PhaseIFall16_81X_upgrade2017_realistic_v26-v1/50000/0AC6BCC6-BBE9-E611-BC3B-02163E019BD4.root'   #-> TTJetsPhase1_1.root
#'/store/mc/PhaseIFall16MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PhaseIFall16PUFlat20to50_PhaseIFall16_81X_upgrade2017_realistic_v26-v1/50000/10525025-01E4-E611-B847-008CFA56D894.root'    #-> TTJetsPhase1_2.root
])

secFiles.extend( [
               ] )
