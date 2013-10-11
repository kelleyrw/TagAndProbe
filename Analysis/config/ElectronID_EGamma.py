import FWCore.ParameterSet.Config as cms
import os

## path the analysis
analysis_path = os.getenv("CMSSW_BASE") + "/src/TagAndProbe/Analysis"

## add the configuration path (if not already there) 
import sys
sys.path.append(analysis_path + "/config")

## process to parse
process = cms.PSet()

## ------------------------------------------------------------------------------------------------------------------- #
## define the input datasets
## ------------------------------------------------------------------------------------------------------------------- #

from datasets import *

## ------------------------------------------------------------------------------------------------------------------- #
## Parameters for the selection and Plot making 
## ------------------------------------------------------------------------------------------------------------------- #

## datasets to process
datasets_el = cms.VPSet(dy_full, single_el)

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
output_label = cms.string("ElectronID_EGamma")

## suffix to print the plots (before the fit)
## blank means do not print
## available options are: eps, png, pdf
suffix = cms.string("png")

## ------------------------------------------------------------------------------------------------------------------- #
## Parameters for the fitting 
## ------------------------------------------------------------------------------------------------------------------- #

## models for pt bins 
pt_models = cms.vstring( 
#          sig pass,        sig fail,      bkg pass,      bkg fail
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # pt0
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # pt1
   	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # pt2
   	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # pt3
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # pt4
	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # pt5
)

## models for eta bins 
eta_models = cms.vstring( 
#          sig pass,        sig fail,      bkg pass,      bkg fail
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta0
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta1
   	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # eta2
   	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta3
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta4
)

## models for phi bins 
phi_models = cms.vstring( 
#          sig pass,        sig fail,      bkg pass,      bkg fail
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi0
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # phi1
   	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # phi2
   	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi3
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi4
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi5
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # phi6
)

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
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # nvtx7
)


## models for pt vs eta bins 
pt_vs_eta_models = cms.vstring( 
#          sig pass,        sig fail,      bkg pass,      bkg fail
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta0, pt0
 	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta1, pt0
 	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta2, pt0
 	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta3, pt0
 	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # eta4, pt0
, 
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta0, pt1
 	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta1, pt1
 	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # eta2, pt1
 	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # eta3, pt1
 	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # eta4, pt1
,                                                                
   	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # eta0, pt2
	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # eta1, pt2 
	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta2, pt2
	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta3, pt2
	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # eta4, pt2
,                                                                
   	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta0, pt3
	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # eta1, pt3
	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta2, pt3
	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # eta3, pt3
	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"       # eta4, pt3
,                                                                
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta0, pt4
	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta1, pt4
	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta2, pt4
	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta3, pt4
	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"    # eta4, pt4
,                                                                
	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta0, pt5
 	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta1, pt5
 	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta2, pt5
 	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # eta3, pt5
 	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"    # eta4, pt5
)

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
   	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # phi0, eta2
	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # phi1, eta2 
	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi2, eta2
	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi3, eta2
	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi4, eta2
 	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential", # phi5, eta2
 	"MCTemplate"   , "MCTemplate"   , "Exponential", "Exponential"  # phi6, eta2
,                                                                
   	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # phi0, eta3
	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # phi1, eta3
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
)

## electrons
## ------------------------------------------ #
electron = cms.PSet(
	## type of lepton 
	lepton_type = cms.string("electron"),

	## output label to give it a unique name
	output_label = output_label,
 
	## suffix to print the plots (before the fit)
	## blank means do not print
	## available options are: eps, png, pdf
	suffix = suffix, 
 
	## max number of events to run on
	max_events = cms.int64(-1),

	## verbosity (for trouble shooting)
	verbose = cms.bool(False),

	## mass range for resonance window
	mass_low       = mass_low,
	mass_high      = mass_high,
	mass_bin_width = mass_bin_width,
	
	# datasets to run on
	datasets = datasets_el,

	## bins for the observables
	## supported pt, eta, phi, and # vertices
	## note: for eta and phi, no negative bins means use |eta| and |phi|, respectively
# 	pt_bins   = cms.vdouble(10, 15),
# 	eta_bins  = cms.vdouble(0, 0.8),
# 	phi_bins  = cms.vdouble(),
# 	nvtx_bins = cms.vdouble(),
	pt_bins   = cms.vdouble(10, 15, 20, 30, 40, 50, 200),
	eta_bins  = cms.vdouble(0, 0.8, 1.4442, 1.566, 2.0, 2.5),
	phi_bins  = cms.vdouble(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.15),
	nvtx_bins = cms.vdouble(0, 5, 10, 15, 20, 25, 30, 35, 40),
	
	## selection 
	## note: numerator and denominator must be one-to-one (e.i. numerator[0] corresponds to denominator[0])
	numerator   = cms.vstring("EGammaNum"  ),
	denominator = cms.vstring("EGammaDenID"),

	# models
	pt_models         = pt_models,
	eta_models        = eta_models,
	phi_models        = phi_models,
	nvtx_models       = nvtx_models,
	pt_vs_eta_models  = pt_vs_eta_models,
	eta_vs_phi_models = eta_vs_phi_models,
)

## process to run to make the plots
## will make a set of plots for each element of the cms.VPSet
## ------------------------------------------ #

process.tnp_make_plots = cms.VPSet(electron)

## ------------------------------------------------------------------------------------------------------------------- #
## Parameters for the fitting 
## ------------------------------------------------------------------------------------------------------------------- #

process.tnp_fit_plots = electron
