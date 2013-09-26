import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.TagAndProbeAnalysis = cms.PSet(
    baby_file   = cms.string('file:/home/users/rwkelley/Data/nfs-7/edm/TTWJets_8TeV-madgraph_AODSIM_PU_S10_START53_V7A-v1_0000_E67C6173-9BDA-E111-BAF3-0030487F1309.root'), ## mandatory
    plot_file   = cms.string('plots/tnp_test.root'), 
    maxEvents   = cms.int32(-1),                     ## optional
)
