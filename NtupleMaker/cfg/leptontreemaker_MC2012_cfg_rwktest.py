import FWCore.ParameterSet.Config as cms

process = cms.Process("SMURF")

#
# general config
#
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START53_V7F::All"
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

#
# lepton maker
#

process.load("Smurf.ProcessingAndSkimming.leptontreemaker_cff")

#
# input
#
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'rfio:/castor/cern.ch/user/e/emanuele/AODSummer12/00514868-B47A-E111-803E-001D0967DDC3.root'
#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/G_Pt-30to50_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/FCF7C5AE-B61A-E211-B10F-00266CF33100.root'
'file:/home/users/rwkelley/Data/nfs-7/edm/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_AODSIM_PU_S10_START53_V7A-v1_0000_00037C53-AAD1-E111-B1BE-003048D45F38.root'
    )
)

process.leptonTreeMaker2012.runFR = cms.bool(False)
process.leptonTreeMaker2012.runGamma = cms.bool(False)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.p = cms.Path(process.leptonTreeMakerSequenceMC * process.leptonTreeMaker2012)

