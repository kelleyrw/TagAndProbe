import FWCore.ParameterSet.Config as cms
import os
import sys
sys.path.append(os.getenv("CMSSW_BASE") + "/src/TagAndProbe/Analysis/config")

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

# models for electron ID
pt_models = cms.vstring( 
#          sig pass,        sig fail,      bkg pass,      bkg fail
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # pt0
	"BreitWignerCB", "BreitWignerCB", "Exponential", "Exponential", # pt1
   	"MCTemplate"   , "MCTemplate"   , "Poly3"      , "Poly3"      , # pt2
   	"MCTemplate"   , "MCTemplate"   , "ErfExp"     , "ErfExp"     , # pt3
    "MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # pt4
	"MCTemplate"   , "MCTemplate"   , "Chebychev"  , "Chebychev"  , # pt5
)

pt_vs_eta_models = cms.vstring( 
#              sig pass,        sig fail,      bkg pass,      bkg fail
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

#pt_vs_eta_models2 = cms.VPSet( 
#		bin_eta0_pt0=cms.PSet(sig_pass=cms.string("BreitWignerCB"), sig_fail=cms.string("BreitWignerCB"), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta1_pt0=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta2_pt0=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta3_pt0=cms.PSet(sig_pass=cms.string("BreitWignerCB"), sig_fail=cms.string("BreitWignerCB"), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta4_pt0=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),                                                                                                                                                               
#		bin_eta0_pt1=cms.PSet(sig_pass=cms.string("BreitWignerCB"), sig_fail=cms.string("BreitWignerCB"), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta1_pt1=cms.PSet(sig_pass=cms.string("BreitWignerCB"), sig_fail=cms.string("BreitWignerCB"), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta2_pt1=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta3_pt1=cms.PSet(sig_pass=cms.string("BreitWignerCB"), sig_fail=cms.string("BreitWignerCB"), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),
#		bin_eta4_pt1=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),                                                                                                                                                               
#		bin_eta0_pt2=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Poly3"      ), bkg_fail=cms.string("Poly3"      )),
#		bin_eta1_pt2=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Poly3"      ), bkg_fail=cms.string("Poly3"      )), 
#		bin_eta2_pt2=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("ErfExp"     ), bkg_fail=cms.string("ErfExp"     )),
#		bin_eta3_pt2=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("ErfExp"     ), bkg_fail=cms.string("ErfExp"     )),
#		bin_eta4_pt2=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Exponential"), bkg_fail=cms.string("Exponential")),                                                                                                                                                               
#		bin_eta0_pt3=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("ErfExp"     ), bkg_fail=cms.string("ErfExp"     )),
#		bin_eta1_pt3=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Poly3"      ), bkg_fail=cms.string("Poly3"      )),
#		bin_eta2_pt3=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("ErfExp"     ), bkg_fail=cms.string("ErfExp"     )),
#		bin_eta3_pt3=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("ErfExp"     ), bkg_fail=cms.string("ErfExp"     )),
#		bin_eta4_pt3=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("ErfExp"     ), bkg_fail=cms.string("ErfExp"     )),                                                                                                                                                               
#		bin_eta0_pt4=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta1_pt4=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta2_pt4=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta3_pt4=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta4_pt4=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),                                                                                                                                                               
#		bin_eta0_pt5=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta1_pt5=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta2_pt5=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta3_pt5=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )),
#		bin_eta4_pt5=cms.PSet(sig_pass=cms.string("MCTemplate"   ), sig_fail=cms.string("MCTemplate"   ), bkg_pass=cms.string("Chebychev"  ), bkg_fail=cms.string("Chebychev"  )) 
#		)
#
#print pt_vs_eta_models2
 
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
	phi_bins  = cms.vdouble(),
	nvtx_bins = cms.vdouble(),
	pt_bins   = cms.vdouble(10, 15, 20, 30, 40, 50, 200),
	eta_bins  = cms.vdouble(0, 0.8, 1.4442, 1.566, 2.0, 2.5),
# 	phi_bins  = cms.vdouble(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.15),
# 	nvtx_bins = cms.vdouble(0, 5, 10, 15, 20, 25, 30, 35, 40),
	
	## selection 
	## note: numerator and denominator must be one-to-one (e.i. numerator[0] corresponds to denominator[0])
	numerator   = cms.vstring("EGammaNum"  ),
	denominator = cms.vstring("EGammaDenID"),

	# models
	pt_vs_eta_models = pt_vs_eta_models,
	pt_models = pt_models,
)

## process to run to make the plots
## will make a set of plots for each element of the cms.VPSet
## ------------------------------------------ #

process.tnp_make_plots = cms.VPSet(electron)

## ------------------------------------------------------------------------------------------------------------------- #
## Parameters for the fitting 
## ------------------------------------------------------------------------------------------------------------------- #

#pt_vs_eta_models = [
#    [#          sig pass,        sig fail,      bkg pass,      bkg fail
#		[cms.string("BreitWignerCB"), cms.string("BreitWignerCB"), cms.string("Exponential"), cms.string("Exponential")], # eta0, pt0
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Exponential"), cms.string("Exponential")], # eta1, pt0
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Exponential"), cms.string("Exponential")], # eta2, pt0
#     	[cms.string("BreitWignerCB"), cms.string("BreitWignerCB"), cms.string("Exponential"), cms.string("Exponential")], # eta3, pt0
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Exponential"), cms.string("Exponential")]  # eta4, pt0
#	],                                                                                                                 
#    [                                                                                                                  
#		[cms.string("BreitWignerCB"), cms.string("BreitWignerCB"), cms.string("Exponential"), cms.string("Exponential")], # eta0, pt1
#     	[cms.string("BreitWignerCB"), cms.string("BreitWignerCB"), cms.string("Exponential"), cms.string("Exponential")], # eta1, pt1
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Exponential"), cms.string("Exponential")], # eta2, pt1
#     	[cms.string("BreitWignerCB"), cms.string("BreitWignerCB"), cms.string("Exponential"), cms.string("Exponential")], # eta3, pt1
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Exponential"), cms.string("Exponential")]  # eta4, pt1
#	],                                                                                                                 
#    [                                                                                                                  
#		[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Poly3"      ), cms.string("Poly3"      )], # eta0, pt2
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Poly3"      ), cms.string("Poly3"      )], # eta1, pt2 
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("ErfExp"     ), cms.string("ErfExp"     )], # eta2, pt2
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("ErfExp"     ), cms.string("ErfExp"     )], # eta3, pt2
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Exponential"), cms.string("Exponential")]  # eta4, pt2
#	],                                                                                                                 
#    [                                                                                                                  
#		[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("ErfExp"     ), cms.string("ErfExp"     )], # eta0, pt3
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Poly3"      ), cms.string("Poly3"      )], # eta1, pt3
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("ErfExp"     ), cms.string("ErfExp"     )], # eta2, pt3
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("ErfExp"     ), cms.string("ErfExp"     )], # eta3, pt3
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("ErfExp"     ), cms.string("ErfExp"     )]  # eta4, pt3
#	],                                                                                                                 
#    [                                                                                                                  
#		[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta0, pt4
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta1, pt4
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta2, pt4
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta3, pt4
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )]  # eta4, pt4
#	],                                                                                                                 
#    [                                                                                                                  
#		[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta0, pt5
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta1, pt5
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta2, pt5
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )], # eta3, pt5
#     	[cms.string("MCTemplate"   ), cms.string("MCTemplate"   ), cms.string("Chebychev"  ), cms.string("Chebychev"  )]  # eta4, pt5
#	]
#]

# print pt_vs_eta_models

process.tnp_fit_plots = electron
