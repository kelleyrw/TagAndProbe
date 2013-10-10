import FWCore.ParameterSet.Config as cms
import os

process = cms.PSet()

## ------------------------------------------------------------------------------------------------------------------- #
## define the input datasets
## ------------------------------------------------------------------------------------------------------------------- #

## path to the lepton trees
lepton_tree_path = "/Users/rwk7t/Data/babies/leptonTrees/tnp_V00-00-00"

## path to the analysis
analysis_path = os.getenv("CMSSW_BASE") + "/src/TagAndProbe/Analysis"

## good run list
run_list = cms.string(analysis_path+'/json/final_19p49fb.txt')


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
data_mu = cms.PSet(
	name     = cms.string("data_double_mu"),
	files    = cms.vstring([lepton_tree_path+'/DoubleMu_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

single_mu = cms.PSet(
	name     = cms.string("data_single_mu"),
	files    = cms.vstring(['/Users/rwk7t/Data/babies/LeptonTree/V00-02-09/SingleMu/merged_Moriond.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

## electron triggered data 
data_el = cms.PSet(
	name     = cms.string("data_double_el"),
	files    = cms.vstring([lepton_tree_path+'/DoubleElectron_Run2012*/*.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)

single_el = cms.PSet(
	name     = cms.string("data_single_el"),
	files    = cms.vstring(['/Users/rwk7t/Data/babies/LeptonTree/V00-02-09/SingleElectron/merged_Moriond.root']),
	is_data  = cms.bool(True),
	run_list = run_list 
)


## ------------------------------------------------------------------------------------------------------------------- #
## Parameters for the selection and Plot making 
## ------------------------------------------------------------------------------------------------------------------- #

## datasets to process
datasets_mu = cms.VPSet(dy_full, single_mu)  # muons
datasets_el = cms.VPSet(dy_full, single_el)  # electrons

## maximum number of events to run on
# max_events = cms.int64(1000)
max_events = cms.int64(-1)

## verbose print out for troubleshooting 
# verbose = cms.bool(True)
verbose = cms.bool(False)

## mass range for resonance window
mass_low       = cms.double(60.0)  # GeV
mass_high      = cms.double(120.0) # GeV
mass_bin_width = cms.double(2.0) # GeV

## path and name to histogram for Pile UP reweighting
pileup_hist_file = cms.string("")
pileup_hist_name = cms.string("")

## output label to give it a unique name
output_label = cms.string("ss2012_test")

## suffix to print the plots (before the fit)
## blank means do not print
## available options are: eps, png, pdf
suffix = cms.string("png")
 
## muons
## ------------------------------------------ #
ss2012_mu = cms.PSet(
	## type of lepton 
	lepton_type = cms.string("muon"),

	## output label to give it a unique name
	output_label = output_label,
 
	## suffix to print the plots (before the fit)
	## blank means do not print
	## available options are: eps, png, pdf
	suffix = suffix, 
 
	## max number of events to run on
	max_events = max_events,

	## verbosity (for trouble shooting)
	verbose = verbose,

	## mass range for resonance window
	mass_low       = mass_low,
	mass_high      = mass_high,
	mass_bin_width = mass_bin_width,
	
	# datasets to run on
	datasets = datasets_mu,

	## bins for the observables
	## supported pt, eta, phi, and # vertices
	## note: for eta and phi, no negative bins means use |eta| and |phi|, respectively
	pt_bins   = cms.vdouble(10, 15, 20, 30, 40, 50, 200),
	eta_bins  = cms.vdouble(0, 1.2, 2.5),
	phi_bins  = cms.vdouble(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.15),
	nvtx_bins = cms.vdouble(0, 5, 10, 15, 20, 25, 30, 35, 40),
	
	## selection 
	## note: numerator and denominator must be one-to-one (e.i. numerator[0] corresponds to denominator[0])
	numerator   = cms.vstring("SameSignMuNum"  , "SameSignMuNum"   , "SameSignMuNum"    ),
	denominator = cms.vstring("SameSignMuDenID", "SameSignMuDenIso", "SameSignMuDenBoth"),
)

## electrons
## ------------------------------------------ #
egamma_el = cms.PSet(
	## type of lepton 
	lepton_type = cms.string("electron"),

	## output label to give it a unique name
	output_label = output_label,
 
	## suffix to print the plots (before the fit)
	## blank means do not print
	## available options are: eps, png, pdf
	suffix = suffix, 
 
	## max number of events to run on
	max_events = max_events,

	## verbosity (for trouble shooting)
	verbose = verbose,

	## mass range for resonance window
	mass_low       = mass_low,
	mass_high      = mass_high,
	mass_bin_width = mass_bin_width,
	
	# datasets to run on
	datasets    = datasets_el,

	## bins for the observables
	## supported pt, eta, phi, and # vertices
	## note: for eta and phi, no negative bins means use |eta| and |phi|, respectively
	pt_bins   = cms.vdouble(10, 15, 20, 30, 40, 50, 200),
	eta_bins  = cms.vdouble(0, 0.8, 1.4442, 1.566, 2.0, 2.5),
	phi_bins  = cms.vdouble(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.15),
	nvtx_bins = cms.vdouble(0, 5, 10, 15, 20, 25, 30, 35, 40),
	
	## selection 
	## note: numerator and denominator must be one-to-one (e.i. numerator[0] corresponds to denominator[0])
	numerator   = cms.vstring("EGammaNum"  , "EGammaNum"   , "EGammaNum"    ),
	denominator = cms.vstring("EGammaDenID", "EGammaDenIso", "EGammaDenBoth"),
)

## process to run to make the plots
## will make a set of plots for each element of the cms.VPSet
## ------------------------------------------ #

process.tnp_make_plots = cms.VPSet(egamma_el, ss2012_mu)

## ------------------------------------------------------------------------------------------------------------------- #
## Parameters for the fitting 
## ------------------------------------------------------------------------------------------------------------------- #
