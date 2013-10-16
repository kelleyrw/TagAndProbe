import FWCore.ParameterSet.Config as cms

process = cms.Process("TNP")

# general config
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_52_V7::All"
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10) )

# lepton maker
process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")

# input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd.unl.edu//store/data/Run2012C/SingleMu/AOD/24Aug2012-v1/00000/F475B31B-03EF-E111-BCC3-002590747DDC.root'
    )
)

process.leptonTreeMaker2012.runFR = cms.bool(False)
process.leptonTreeMaker2012.runGamma = cms.bool(False)
process.leptonTreeMaker2012.pfJetCorrectorL1FastL2L3 = cms.string('ak5PFL1FastL2L3Residual')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.muonFilters * process.leptonTreeMakerSequenceData2012 * process.leptonTreeMaker2012)
