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
## Parameters for the comparison 
## ------------------------------------------------------------- #

# electron ID comparison
el_id = cms.PSet(
	
	## verbosity (for trouble shooting)
	verbose = cms.bool(False),

	## label to give unique name
	output_label = cms.string("compare"),

	## path to the efficiency results to compare 
	eff_results_path = cms.string(analysis_path+"/plots/ElectronID_EGammaMediumWP/electron/EGammaMediumWPDenID_EGammaMediumWPNum"),

	## first result to compare
	dataset1 = single_el,

	## second result to compare
	dataset2 = dy_full,

	## suffix to print the plots (before the fit)
	## blank means do not print
	## available options are: eps, png, pdf
	suffix = cms.string("png"), 
)

# muon ID comparison
mu_id = cms.PSet(
	
	## verbosity (for trouble shooting)
	verbose = cms.bool(False),

	## label to give unique name
	output_label = cms.string("compare"),

	## path to the efficiency results to compare 
	eff_results_path = cms.string(analysis_path+"/plots/MuID_MuTightWP/muon/MuTightWPDenID_MuTightWPNum"),

	## first result to compare
	dataset1 = single_mu,

	## second result to compare
	dataset2 = dy_full,

	## suffix to print the plots (before the fit)
	## blank means do not print
	## available options are: eps, png, pdf
	suffix = cms.string("png"), 
)

# a vector of all the comparison PSets --> will do them all.
process.tnp_compare = cms.VPSet(el_id, mu_id)
