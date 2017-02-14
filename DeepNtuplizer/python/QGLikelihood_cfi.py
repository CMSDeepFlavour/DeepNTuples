import FWCore.ParameterSet.Config as cms
import os

qgDatabaseVersion = 'cmssw8020_v2'

databasepath=os.environ['CMSSW_BASE']+'/src/DeepNTuples/DeepNtuplizer/data/QGL_cmssw8020_v2.db'

from CondCore.CondDB.CondDB_cfi import *
QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      CondDB.DBParameters,
      toGet = cms.VPSet(),
      connect = cms.string('sqlite_file:'+databasepath),
)

for type in ['AK4PFchs','AK4PFchs_antib']:
  QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
    record = cms.string('QGLikelihoodRcd'),
    tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
    label  = cms.untracked.string('QGL_'+type)
)))
