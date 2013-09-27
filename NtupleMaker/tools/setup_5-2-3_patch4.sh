#!/bin/bash

# -------------------------------------------------------------------------- #
# Simple Script to checkout and compile the /Smurf/ProcessingAndSkimming code
# needed for the leptonTreeMaker used by the this Tag and Probe analysis.
# The tags used were taken from: 
#   https://github.com/drkovalskyi/Smurf/blob/head/ProcessingAndSkimming/README 
# -------------------------------------------------------------------------- #

# metadata
# -------------------------------------------------------------------------- #

# not used yet
tag=""

# checkouts
# -------------------------------------------------------------------------- #

pushd $CMSSW_BASE/src

# Smurf Code
git clone https://github.com/drkovalskyi/Smurf.git

# Egamma ID code (this has been moved to CMSSW but in later releases)
cvs co -r V00-00-06 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
cvs up -r CutBasedId_V00-00-05 EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleId.cc
cvs up -r CutBasedId_V00-00-05 EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h
cvs up -r 1.3 EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h

# for Si Xie thesis MVA
# using combined ID+Isolation
cvs co -r HWWID_V00-00-00 HiggsAnalysis/HiggsToWW2Leptons

# 2012 HZZ ID+Iso combined MVA
cvs co -r R2012A_V2 -d Muon/MuonAnalysisTools UserCode/sixie/Muon/MuonAnalysisTools

# needed for the jet id
cvs co -r V00-02-05 -d CMGTools/External UserCode/CMG/CMGTools/External

# emulate pf no pileup (these tags are for 52X - may run into problems in earlier releases)
cvs co -r V00-03-11 CommonTools/ParticleFlow 
cvs co -r V15-01-06 RecoParticleFlow/PFProducer 

if [ $? -eq 0 ]; then
    echo "\n[TNP] checkout complete!\n"
else
    echo "\n[TNP] checkout failed!\n"
fi

# compile
# -------------------------------------------------------------------------- #

cmsenv
scramv1 b -j20

if [ $? -eq 0 ]; then
    echo "\n[TNP] compiled successfully!\n"
else
    echo "\n[TNP] compile failure!\n"
fi

# done
popd
