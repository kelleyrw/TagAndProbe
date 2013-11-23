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
	lepton_type = cms.string("electron"),

	## output label to give it a unique name
	output_label = cms.string("ElectronID_EGammaMediumWP"),
 
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
#	datasets = cms.VPSet(dy_full, single_el),
#       datasets = cms.VPSet(dy_full, double_el),                                     
        datasets = cms.VPSet(dy_full2, double_el),

	## bins for the observables
	## supported pt, eta, phi, and # vertices
	## note: for eta and phi, no negative bins means use |eta| and |phi|, respectively
	pt_bins   = cms.vdouble(10, 15, 20, 30, 40, 50, 200),
	eta_bins  = cms.vdouble(0, 0.8, 1.4442, 1.566, 2.0, 2.5),
#	phi_bins  = cms.vdouble(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.15),
#	nvtx_bins = cms.vdouble(0, 5, 10, 15, 20, 25, 30, 35, 40),
 	phi_bins  = cms.vdouble(),               
 	nvtx_bins = cms.vdouble(),                    
	
	## selection (from Measurements.h/cc) 
	numerator   = cms.string("EGammaMediumWPNum"  ),
#	denominator = cms.string("EGammaMediumWPDenID"),
	denominator = cms.string("EGammaMediumWPDenBoth"),

	## Parameters for the fitting 
	## --------------------------------------------------------- #

	## models for pt bins 
	pt_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"MCTemplate", "MCTemplate", "Exponential", "Exponential", # pt0
		"MCTemplate", "MCTemplate"   , "Exponential", "Exponential"     , # pt1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # pt2
		"MCTemplate"   , "MCTemplate"   , "Exponential"     , "Exponential"  , # pt3
		"MCTemplate"   , "MCTemplate"   , "Exponential"  , "Exponential"  , # pt4
		"MCTemplate"   , "MCTemplate"   , "Exponential"  , "Exponential"  , # pt5
	),
	
	## models for eta bins 
	eta_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"MCTemplate"   , "MCTemplate", "Exponential", "Exponential", # eta0
		"MCTemplate"   , "MCTemplate", "Exponential", "Exponential", # eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential"      , "Exponential"      , # eta2  # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplate"   , "MCTemplate"   , "Exponential"  , "Exponential"  , # eta3
		"MCTemplate"   , "MCTemplate"   , "Exponential"  , "Exponential"  , # eta4
	),
	
	## models for phi bins 
	phi_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi0
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi1
		"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # phi2
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi3
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi5
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi6
	),
	
	## models for nvtx bins 
	nvtx_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx0
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # nvtx1
		"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # nvtx2
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # nvtx3
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # nvtx4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # nvtx5
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # nvtx6
		"MCTemplate"   , "BreitWignerCB", "Chebychev"  , "Chebychev"  , # nvtx7
	),
	
	## models for pt vs eta bins 
#	pt_vs_eta_models = cms.vstring( 
#	#          sig pass,        sig fail,      bkg pass,      bkg fail
#		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta0, pt0
#		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta1, pt0
#		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta2, pt0 # crack electrons (1.4442 < |eta| < 1.566)
#		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta3, pt0
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "Exponential"  # eta4, pt0
#	, 
#		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta0, pt1
#		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta1, pt1
#		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta2, pt1 # crack electrons (1.4442 < |eta| < 1.566)
#		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta3, pt1
#		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # eta4, pt1
#	,
#		"MCTemplate"   , "MCTemplate"   , "Exponential", "ErfExp"     , # eta0, pt2
#		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta1, pt2 
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta2, pt2 # crack electrons (1.4442 < |eta| < 1.566)
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta3, pt2
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"       # eta4, pt2
#	,
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta0, pt3
#		"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # eta1, pt3
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta2, pt3 # crack electrons (1.4442 < |eta| < 1.566)
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "Chebychev"  , # eta3, pt3
#		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"       # eta4, pt3
#	,
#		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta0, pt4
#		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta1, pt4
#		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta2, pt4 # crack electrons (1.4442 < |eta| < 1.566)
#		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta3, pt4
#		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"    # eta4, pt4
#	,
#		"BreitWignerCB", "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta0, pt5
#		"BreitWignerCB", "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta1, pt5
#		"BreitWignerCB", "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta2, pt5 # crack electrons (1.4442 < |eta| < 1.566)
#		"BreitWignerCB", "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta3, pt5
#		"BreitWignerCB", "MCTemplate"   , "Chebychev"  , "Chebychev"    # eta4, pt5
#	),
#	
	pt_vs_eta_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta0, pt0
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta1, pt0
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta2, pt0 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta3, pt0
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential"  # eta4, pt0
	,
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta0, pt1
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta1, pt1
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta2, pt1 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta3, pt1
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential"  # eta4, pt1
	,
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta0, pt2
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta1, pt2 
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta2, pt2 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta3, pt2
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential"  # eta4, pt2
	,
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta0, pt3
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta1, pt3
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta2, pt3 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta3, pt3
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential"  # eta4, pt3
	,
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta0, pt4
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta1, pt4
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta2, pt4 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta3, pt4
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential"  # eta4, pt4
	,
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta0, pt5
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta1, pt5
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta2, pt5 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential", # eta3, pt5
		"MCTemplateCB"   , "MCTemplateCB"   , "Exponential", "Exponential"  # eta4, pt5
	),
	
	## models for eta vs phi bins 
	eta_vs_phi_models = cms.vstring( 
	#          sig pass,        sig fail,      bkg pass,      bkg fail
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi0, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi1, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi2, eta0
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi3, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi4, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta0
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta0
	, 
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi0, eta1
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi1, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi2, eta1
		"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi3, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi4, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta1
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta0
	,
		"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # phi0, eta2 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # phi1, eta2 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi2, eta2 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi3, eta2 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi4, eta2 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta2 # crack electrons (1.4442 < |eta| < 1.566)
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta2 # crack electrons (1.4442 < |eta| < 1.566)
	,
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi0, eta3
		"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Exponential", # phi1, eta3
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi2, eta3
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi3, eta3
		"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi4, eta3
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta3
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta3
	,
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi0, eta4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi1, eta4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi2, eta4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi3, eta4
		"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi4, eta4
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta4
		"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta4
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
process.tnp_make_plots.suffix = cms.string("") 

## ------------------------------------------------------------- #
## process to fit the plots and extract the efficiency
## ------------------------------------------------------------- #

process.tnp_eff_plots = tnp_parameters
process.tnp_make_plots.suffix = cms.string("png") 
