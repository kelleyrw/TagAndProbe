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
	verbose = cms.bool(True),

	## label to give unique name
	output_label = cms.string("compare"),

	## path to the efficiency results to compare 
	eff_results_path = cms.string(analysis_path+"/plots/ElectronID_EGammaMediumWP/electron/EGammaMediumWPDenBoth_EGammaMediumWPNum"),

	## first result to compare
	dataset1 = double_el,

	## second result to compare
	dataset2 = dy_full2,

	## suffix to print the plots (before the fit)
	## blank means do not print
	## available options are: eps, png, pdf
	suffix = cms.string("png"), 
)

# a vector of all the comparison PSets --> will do them all.
process.tnp_compare = cms.VPSet(el_id)
