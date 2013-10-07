import FWCore.ParameterSet.Config as cms

process = cms.PSet()

# used to make the plots
process.TagAndProbeAnalysis = cms.PSet(
    baby_file   = cms.string('/Users/rwk7t/Data/babies/leptonTrees/tnp_V00-00-00/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/merged_leptontree.root'),
    plot_file   = cms.string('plots/tnp_test.root'), 
    maxEvents   = cms.int32(-1),
)

# use to perform the comparison
