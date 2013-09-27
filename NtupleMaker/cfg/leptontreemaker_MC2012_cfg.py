import FWCore.ParameterSet.Config as cms

process = cms.Process("TNP")

# general config
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START53_V7F::All"
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

# lepton maker
process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")

# input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		#'file:/home/users/rwkelley/Data/nfs-7/edm/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_AODSIM_PU_S10_START53_V7A-v1_0000_00037C53-AAD1-E111-B1BE-003048D45F38.root'
		'root://xrootd.unl.edu//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root'
    )
)

process.leptonTreeMaker2012.runFR = cms.bool(False)
process.leptonTreeMaker2012.runGamma = cms.bool(False)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.leptonTreeMakerSequenceMC * process.leptonTreeMaker2012)

