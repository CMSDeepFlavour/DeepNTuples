#taken from cmsDriver.py
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: pat --filein file:t3mu.root --fileout file:t3mu_MINIAODSIM.root --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions 90X_upgrade2023_realistic_v1 --step PAT --era Phase2C2_timing --geometry Extended2023D4 --python_filename MiniAODTEST_cfg.py --beamspot HLLHC14TeV --customise Configuration/DataProcessing/Utils.addMonitoring -n 100 --no_exec
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

def customize(process):
   process.load('Configuration.StandardSequences.Services_cff')
   process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
   process.load('FWCore.MessageService.MessageLogger_cfi')
   process.load('Configuration.EventContent.EventContent_cff')
   process.load('SimGeneral.MixingModule.mixNoPU_cfi')
   process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
   process.load('Configuration.StandardSequences.MagneticField_cff')
   process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
   process.load('Configuration.StandardSequences.EndOfProcess_cff')
   process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
   
   # Other statements
   from Configuration.AlCa.GlobalTag import GlobalTag
   process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2023_realistic_v1', '')
   
   # Path and EndPath definitions
   process.endjob_step = cms.EndPath(process.endOfProcess)

   # Schedule definition
   process.schedule = cms.Schedule(
      process.endjob_step
      )
   
   from FWCore.ParameterSet.Utilities import convertToUnscheduled
   process=convertToUnscheduled(process)
   process.load('Configuration.StandardSequences.PATMC_cff')
   from FWCore.ParameterSet.Utilities import cleanUnscheduled
   process=cleanUnscheduled(process)

   # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
   from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 
   process = miniAOD_customizeAllMC(process)

   # Add early deletion of temporary data products to reduce peak memory need
   from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
   process = customiseEarlyDelete(process)
   return process
