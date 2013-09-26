import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('file:/home/users/rwkelley//Data/nfs-7/edm/TTWJets_8TeV-madgraph_AODSIM_PU_S10_START53_V7A-v1_0000_E67C6173-9BDA-E111-BAF3-0030487F1309.root'), ## mandatory
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(10),                            ## optional
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('plots/analyzeFWLiteHistograms.root'),  ## mandatory
)

process.muonAnalyzer = cms.PSet(
    ## input specific for this analyzer
    muons = cms.InputTag('muons')
)
