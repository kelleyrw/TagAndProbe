import FWCore.ParameterSet.Config as cms
import os

## path the analysis (THIS SHOULD NOT CHANGE)
analysis_path = os.getenv("CMSSW_BASE") + "/src/TagAndProbe/Analysis"

## add the configuration path (THIS SHOULD NOT CHANGE)
import sys
sys.path.append(analysis_path + "/config")

## process to parse (THIS SHOULD NOT CHANGE)
process = cms.PSet()

## ------------------------------------------------------------- #
## define the input datasets
## ------------------------------------------------------------- #

from datasets import *

## ------------------------------------------------------------- #
## Parameters for the selection, plot making, and fitting 
## ------------------------------------------------------------- #

tnp_parameters = cms.PSet(
	
	## Parameters for the selection and mass plots 
	## --------------------------------------------------------- #
	
	## type of lepton 
	lepton_type = cms.string("muon"),

	## output label to give it a unique name
	output_label = cms.string("SameSign"),
 
	## suffix to print the plots (before the fit)
	## blank means do not print
	## available options are: eps, png, pdf
	suffix = cms.string("png"), 

	## pile re-weighting histogram
	pileup_hist_file = cms.string(analysis_path + "/data/puWeights_Summer12_53x_Observed.root"),
	pileup_hist_name = cms.string("puWeights"),
 
	## max number of events to run on
	max_events = cms.int64(-1),

	## verbosity (for trouble shooting)
	verbose = cms.bool(False),

	## mass range for resonance window
	mass_low       = cms.double(60.0),  # GeV
	mass_high      = cms.double(120.0), # GeV
	mass_bin_width = cms.double(2.0),   # GeV
	
	# datasets to run on
	datasets = cms.VPSet(dy_full, single_mu),

	## bins for the observables
	## supported pt, eta, phi, and # vertices
	## note: for eta and phi, no negative bins means use |eta| and |phi|, respectively
	pt_bins   = cms.vdouble(10, 15, 20, 30, 40, 50, 200),
	eta_bins  = cms.vdouble(0, 1.1, 2.5),
	phi_bins  = cms.vdouble(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.15),
	nvtx_bins = cms.vdouble(0, 5, 10, 15, 20, 25, 30, 35, 40),
	
	## selection (from Measurements.h/cc) 
	numerator   = cms.string("SameSignNum"    ),
	denominator = cms.string("SameSignDenBoth"),

	## Parameters for the fitting 
	## --------------------------------------------------------- #

	## models for pt bins 
	pt_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # pt0
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # pt1
		"BreitWignerCB", "MCTemplate"   , "ErfExp"     , "ErfExp"     , # pt2
		"BreitWignerCB", "MCTemplate"   , "ErfExp"     , "ErfExp"     , # pt3
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # pt4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # pt5
	),
	
	## models for eta bins 
	eta_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"MCTemplate"   , "BreitWignerCB", "Chebychev", "Chebychev", # eta0
		"MCTemplate"   , "BreitWignerCB", "Chebychev", "Chebychev", # eta2
	),
	
	## models for phi bins 
	phi_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # phi0
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # phi1
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # phi2
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # phi3
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # phi4
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # phi5
		"BreitWignerCB", "MCTemplate"   , "Exponential", "Exponential", # phi6
	),
	
	## models for nvtx bins 
	nvtx_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx0
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx1
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx2
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx3
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx4
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx5
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx6
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx7
	),
	
	## models for pt vs eta bins 
	pt_vs_eta_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta0, pt0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # eta1, pt0
	,                                                                                     
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta0, pt1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # eta1, pt1
	,                                                                                     
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta0, pt2
		"BreitWignerCB", "BreitWignerCB", "ErfExp"     , "ErfExp"       # eta1, pt2 
	,                                                                                     
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta0, pt3
		"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"        # eta1, pt3
	,                                                                                     
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta0, pt4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"    # eta1, pt4
	,                                                                                     
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta0, pt5
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"    # eta1, pt5
	),
	
	## models for eta vs phi bins 
	eta_vs_phi_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi0, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi1, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi2, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi3, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi4, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta0
	, 
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi0, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi1, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi2, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi3, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi4, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta0
	),
)


## ------------------------------------------------------------- #
## hisogram files to use with MC Templates (if used) 
## ------------------------------------------------------------- #

tnp_parameters.mc_template_file = cms.string("%s/plots/%s/%s/%s_%s/%s.root" % (analysis_path, 
                                                                               tnp_parameters.output_label.value(),
                                                                               tnp_parameters.lepton_type.value(),
                                                                               tnp_parameters.denominator.value(), 
                                                                               tnp_parameters.numerator.value(), 
                                                                               dy_full.name.value()))
	
## ------------------------------------------------------------- #
## process to make the plots
## will make a set of plots for each element of the cms.VPSet
## ------------------------------------------------------------- #

process.tnp_make_plots = tnp_parameters

## ------------------------------------------------------------- #
## process to fit the plots and extract the efficiency
## ------------------------------------------------------------- #

process.tnp_eff_plots = tnp_parameters
