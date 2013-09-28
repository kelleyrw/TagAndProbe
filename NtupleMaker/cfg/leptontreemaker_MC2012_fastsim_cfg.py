import FWCore.ParameterSet.Config as cms

process = cms.Process("TNP")

# general config
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START52_V9::All"
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# lepton maker
process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")
process.leptonTreeMaker.rhoIsoAllInputTag     = cms.InputTag("kt6PFJets","rho") # override hardcoded "RECO" process name
process.leptonTreeMaker2011.rhoIsoAllInputTag = cms.InputTag("kt6PFJets","rho") # override hardcoded "RECO" process name
process.leptonTreeMaker2012.rhoIsoAllInputTag = cms.InputTag("kt6PFJets","rho") # override hardcoded "RECO" process name

# input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		'file:/home/users/rwkelley/Data/nfs-7/edm/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_AODSIM_START52_V9_FSIM-v1_00000_001985EC-07EE-E111-8517-00261894386A.root'
# 		'root://xrootd.unl.edu//store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/START52_V9_FSIM-v1/00000/001985EC-07EE-E111-8517-00261894386A.root'
    )
)

process.leptonTreeMaker2012.runFR = cms.bool(False)
process.leptonTreeMaker2012.runGamma = cms.bool(False)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.leptonTreeMakerSequenceMC * process.leptonTreeMaker2012)

