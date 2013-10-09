import FWCore.ParameterSet.Config as cms
import os

process = cms.PSet()

# used to make the plots
# dataset = cms.PSet(
# 	files    = cms.vstring(['/Users/rwk7t/Data/babies/leptonTrees/tnp_V00-00-00/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/merged_leptontree.root']),
# 	is_data  = cms.bool(False),
# 	run_list = cms.string('')
# )

# define the input datasets
lepton_tree_path = "/Users/rwk7t/Data/babies/leptonTrees/tnp_V00-00-00"
analysis_path = os.getenv("TNP")
process.datasets = cms.VPSet([
	cms.PSet(name = cms.string("dy_full"), files= cms.vstring([lepton_tree_path+'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/*.root']), is_data = cms.bool(False), run_list = cms.string('')),
	cms.PSet(name = cms.string("dy_fast"), files= cms.vstring([lepton_tree_path+'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-START52_V9_FSIM-v1_AODSIM/*.root'])         , is_data = cms.bool(False), run_list = cms.string('')),
	cms.PSet(name = cms.string("data_mu"), files= cms.vstring([lepton_tree_path+'/DoubleMu_Run2012*/*.root'      ]), is_data = cms.bool(True), run_list = cms.string(analysis_path+'/json/final_19p49fb.txt')),
	cms.PSet(name = cms.string("data_el"), files= cms.vstring([lepton_tree_path+'/DoubleElectron_Run2012*/*.root']), is_data = cms.bool(True), run_list = cms.string(analysis_path+'/json/final_19p49fb.txt'))
	]
)

process.tnp_make_plots = cms.PSet(
   max_events = cms.int32(100),
   verbose    = cms.bool(True)
)



# use to perform the comparison
