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
lepton_tree_path = "/tas/gzevi/TagAndProbe/crab/merged"
#+ lepton_tree_tag

## good run list
run_list = cms.string(analysis_path + "/json/final_19p49fb.txt")

## DY fullsim
dy_full = cms.PSet(
	name     = cms.string("dy_full"),
	title    = cms.string("DY fullsim"),
#	files    = cms.vstring([lepton_tree_path+'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/*.root']),
        files    = cms.vstring([lepton_tree_path+'/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/*.root']),
	is_data  = cms.bool(False),
	run_list = cms.string('')
)

dy_full2 = cms.PSet(
	name     = cms.string("dy_full2"),
	title    = cms.string("DY fullsim 2"),
	files    = cms.vstring([lepton_tree_path+'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/*.root']),
#        files    = cms.vstring([lepton_tree_path+'/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/*.root']),
	is_data  = cms.bool(False),
	run_list = cms.string('')
)

## DY fastsim
dy_fast = cms.PSet(
	name     = cms.string("dy_fast"),
	title    = cms.string("DY faststim"),
	files    = cms.vstring([lepton_tree_path+'/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-START52_V9_FSIM-v1_AODSIM/*.root']),
	is_data  = cms.bool(False),
	run_list = cms.string('')
)

## muon triggered data 
double_mu = cms.PSet(
	name     = cms.string("data_double_mu"),
	title    = cms.string("DoubleMu_Run2012"),
	files    = cms.vstring([lepton_tree_path+'/DoubleMu_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

single_mu = cms.PSet(
	name     = cms.string("data_single_mu"),
	title    = cms.string("SingleMu_Run2012"),
	files    = cms.vstring([lepton_tree_path+'/SingleMu_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

## electron triggered data 
double_el = cms.PSet(
	name     = cms.string("data_double_el"),
	title    = cms.string("DoubleElectron_Run2012"),
#	files    = cms.vstring([lepton_tree_path+'/DoubleElectron_Run2012*/*.root']),
        files    = cms.vstring([lepton_tree_path+'/DoubleElectron_Run2012*22Jan2013*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

single_el = cms.PSet(
	name     = cms.string("data_single_el"),
	title    = cms.string("SingleElectron_Run2012"),
#	files    = cms.vstring(['/nfs-7/userdata/rwkelley/lepton_trees/V00-02-09/SingleElectron/merged_Moriond.root']),
	#files    = cms.vstring([lepton_tree_path+'/SingleElectron_Run2012*/*.root']),
        files    = cms.vstring([lepton_tree_path+'/SingleElectron_Run2012*22Jan2013*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)
