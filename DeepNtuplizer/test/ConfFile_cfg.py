import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/e/ebols/public/cmsswPreProcessing_RECO_wTagInfo.root'
    )
)

process.Ntupler = cms.EDAnalyzer('Ntupler',
                              onJetsrc = cms.InputTag("hltPFJetForBtag"),
                              onJetCalosrc = cms.InputTag("hltAK4CaloJetsCorrected"),
                              offJetsrc = cms.InputTag("ak4PFJetsCHS"),
                              onShallowsrc = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsInfos"),
                              offShallowsrc = cms.InputTag("pfDeepCSVTagInfos"),
                              onBtagsrc = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsPF:probb"),
                              offBtagsrc = cms.InputTag("pfDeepCSVJetTags:probb"),
                              offBBtagsrc = cms.InputTag("pfDeepCSVJetTags:probbb"),
                              onCtagsrc = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsPF:probc"),
                              offCtagsrc = cms.InputTag("pfDeepCSVJetTags:probc"),
                              onUDSGtagsrc = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsPF:probudsg"),
                              offUDSGtagsrc = cms.InputTag("pfDeepCSVJetTags:probudsg"),
                              onCSVtagsrc = cms.InputTag("hltCombinedSecondaryVertexBJetTagsPF"),
                              offCSVtagsrc = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                              onCSVCalotagsrc = cms.InputTag("hltCombinedSecondaryVertexBJetTagsCalo"),
                              onBCalotagsrc = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsCalo:probb"),
                              onCCalotagsrc = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsCalo:probc"),
                              onUDSGCalotagsrc = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsCalo:probudsg")
                             )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('/afs/cern.ch/work/e/ebols/public/ntuple.root')
                                   )


process.p = cms.Path(process.Ntupler)
