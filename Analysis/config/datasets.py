import FWCore.ParameterSet.Config as cms
import os

## path the analysis
analysis_path = os.getenv("CMSSW_BASE") + "/src/TagAndProbe/Analysis"

## add the configuration path
import sys
sys.path.append(analysis_path + "/config")

## ------------------------------------------------------------------------------------------------------------------- #
## define the input datasets
## ------------------------------------------------------------------------------------------------------------------- #

## path to the lepton trees
lepton_tree_tag  = "tnp_V00-00-00"
lepton_tree_path = "/nfs-7/userdata/rwkelley/lepton_trees/" + lepton_tree_tag

## good run list
run_list = cms.string(analysis_path + "/json/final_19p49fb.txt")

## DY fullsim
dy_full = cms.PSet(
	name     = cms.string("dy_full"),
	files    = cms.vstring([lepton_tree_path+'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/*.root']),
	is_data  = cms.bool(False),
	run_list = cms.string('')
)

## DY fastsim
dy_fast = cms.PSet(
	name     = cms.string("dy_fast"),
	files    = cms.vstring([lepton_tree_path+'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-START52_V9_FSIM-v1_AODSIM/*.root']),
	is_data  = cms.bool(False),
	run_list = cms.string('')
)

## muon triggered data 
double_mu = cms.PSet(
	name     = cms.string("data_double_mu"),
	files    = cms.vstring([lepton_tree_path+'/DoubleMu_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

single_mu = cms.PSet(
	name     = cms.string("data_single_mu"),
	files    = cms.vstring([lepton_tree_path+'/SingleMu_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

## electron triggered data 
double_el = cms.PSet(
	name     = cms.string("data_double_el"),
	files    = cms.vstring([lepton_tree_path+'/DoubleElectron_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

single_el = cms.PSet(
	name     = cms.string("data_single_el"),
	files    = cms.vstring([lepton_tree_path+'/SingleElectron_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)
