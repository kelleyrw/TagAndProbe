// -*- C++ -*-
#ifndef LeptonTree_H
#define LeptonTree_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class LeptonTree {
private: 
protected: 
	unsigned int index;
	unsigned int	event_;
	TBranch *event_branch;
	bool event_isLoaded;
	unsigned int	run_;
	TBranch *run_branch;
	bool run_isLoaded;
	unsigned int	lumi_;
	TBranch *lumi_branch;
	bool lumi_isLoaded;
	float	rnd_;
	TBranch *rnd_branch;
	bool rnd_isLoaded;
	unsigned int	nvtx_;
	TBranch *nvtx_branch;
	bool nvtx_isLoaded;
	unsigned int	npu_;
	TBranch *npu_branch;
	bool npu_isLoaded;
	unsigned int	npuPlusOne_;
	TBranch *npuPlusOne_branch;
	bool npuPlusOne_isLoaded;
	unsigned int	npuMinusOne_;
	TBranch *npuMinusOne_branch;
	bool npuMinusOne_isLoaded;
	unsigned int	leptonSelection_;
	TBranch *leptonSelection_branch;
	bool leptonSelection_isLoaded;
	unsigned int	eventSelection_;
	TBranch *eventSelection_branch;
	bool eventSelection_isLoaded;
	float	rhoIsoAll_;
	TBranch *rhoIsoAll_branch;
	bool rhoIsoAll_isLoaded;
	float	rhoIsoAllCentral_;
	TBranch *rhoIsoAllCentral_branch;
	bool rhoIsoAllCentral_isLoaded;
	float	rhoIsoNeutral_;
	TBranch *rhoIsoNeutral_branch;
	bool rhoIsoNeutral_isLoaded;
	float	tagAndProbeMass_;
	TBranch *tagAndProbeMass_branch;
	bool tagAndProbeMass_isLoaded;
	bool	tagAndProbeIsRandom_;
	TBranch *tagAndProbeIsRandom_branch;
	bool tagAndProbeIsRandom_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *tag_;
	TBranch *tag_branch;
	bool tag_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *probe_;
	TBranch *probe_branch;
	bool probe_isLoaded;
	int	qTag_;
	TBranch *qTag_branch;
	bool qTag_isLoaded;
	int	qProbe_;
	TBranch *qProbe_branch;
	bool qProbe_isLoaded;
	float	scale1fb_;
	TBranch *scale1fb_branch;
	bool scale1fb_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *jet1_;
	TBranch *jet1_branch;
	bool jet1_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *jet2_;
	TBranch *jet2_branch;
	bool jet2_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *jet3_;
	TBranch *jet3_branch;
	bool jet3_isLoaded;
	float	met_;
	TBranch *met_branch;
	bool met_isLoaded;
	float	metPhi_;
	TBranch *metPhi_branch;
	bool metPhi_isLoaded;
	float	trackMet_;
	TBranch *trackMet_branch;
	bool trackMet_isLoaded;
	float	trackMetPhi_;
	TBranch *trackMetPhi_branch;
	bool trackMetPhi_isLoaded;
	unsigned int	njets_;
	TBranch *njets_branch;
	bool njets_isLoaded;
	unsigned int	nleps_;
	TBranch *nleps_branch;
	bool nleps_isLoaded;
	float	hltPrescale_;
	TBranch *hltPrescale_branch;
	bool hltPrescale_isLoaded;
	float	sumet_;
	TBranch *sumet_branch;
	bool sumet_isLoaded;
	float	metSig_;
	TBranch *metSig_branch;
	bool metSig_isLoaded;
	float	mt_;
	TBranch *mt_branch;
	bool mt_isLoaded;
	float	dPhiProbeJet1_;
	TBranch *dPhiProbeJet1_branch;
	bool dPhiProbeJet1_isLoaded;
	float	chargedEmFracJet1_;
	TBranch *chargedEmFracJet1_branch;
	bool chargedEmFracJet1_isLoaded;
	float	neutralEmFracJet1_;
	TBranch *neutralEmFracJet1_branch;
	bool neutralEmFracJet1_isLoaded;
	float	chargedMuFracJet1_;
	TBranch *chargedMuFracJet1_branch;
	bool chargedMuFracJet1_isLoaded;
	float	electronHWW2011MVA_;
	TBranch *electronHWW2011MVA_branch;
	bool electronHWW2011MVA_isLoaded;
	float	egammaPOG2012MVA_;
	TBranch *egammaPOG2012MVA_branch;
	bool egammaPOG2012MVA_isLoaded;
	float	muonHZZ2012IsoRingsMVA_;
	TBranch *muonHZZ2012IsoRingsMVA_branch;
	bool muonHZZ2012IsoRingsMVA_isLoaded;
	unsigned int	vetoId_;
	TBranch *vetoId_branch;
	bool vetoId_isLoaded;
	unsigned int	looseId_;
	TBranch *looseId_branch;
	bool looseId_isLoaded;
	unsigned int	mediumId_;
	TBranch *mediumId_branch;
	bool mediumId_isLoaded;
	unsigned int	tightId_;
	TBranch *tightId_branch;
	bool tightId_isLoaded;
	float	pfmva_;
	TBranch *pfmva_branch;
	bool pfmva_isLoaded;
	float	sceta_;
	TBranch *sceta_branch;
	bool sceta_isLoaded;
	float	scphi_;
	TBranch *scphi_branch;
	bool scphi_isLoaded;
	float	scenergy_;
	TBranch *scenergy_branch;
	bool scenergy_isLoaded;
	bool	chargesAgree_;
	TBranch *chargesAgree_branch;
	bool chargesAgree_isLoaded;
	float	eopin_;
	TBranch *eopin_branch;
	bool eopin_isLoaded;
	float	ooemoop_;
	TBranch *ooemoop_branch;
	bool ooemoop_isLoaded;
	float	fbrem_;
	TBranch *fbrem_branch;
	bool fbrem_isLoaded;
	float	detain_;
	TBranch *detain_branch;
	bool detain_isLoaded;
	float	dphiin_;
	TBranch *dphiin_branch;
	bool dphiin_isLoaded;
	float	hoe_;
	TBranch *hoe_branch;
	bool hoe_isLoaded;
	float	hoetow_;
	TBranch *hoetow_branch;
	bool hoetow_isLoaded;
	float	sieie_;
	TBranch *sieie_branch;
	bool sieie_isLoaded;
	float	d0vtx_;
	TBranch *d0vtx_branch;
	bool d0vtx_isLoaded;
	float	dzvtx_;
	TBranch *dzvtx_branch;
	bool dzvtx_isLoaded;
	bool	vfitprob_;
	TBranch *vfitprob_branch;
	bool vfitprob_isLoaded;
	float	mhit_;
	TBranch *mhit_branch;
	bool mhit_isLoaded;
	float	ecaliso_;
	TBranch *ecaliso_branch;
	bool ecaliso_isLoaded;
	float	hcaliso_;
	TBranch *hcaliso_branch;
	bool hcaliso_isLoaded;
	float	trkiso_;
	TBranch *trkiso_branch;
	bool trkiso_isLoaded;
	float	pfemiso03_;
	TBranch *pfemiso03_branch;
	bool pfemiso03_isLoaded;
	float	pfchiso03_;
	TBranch *pfchiso03_branch;
	bool pfchiso03_isLoaded;
	float	pfnhiso03_;
	TBranch *pfnhiso03_branch;
	bool pfnhiso03_isLoaded;
	float	pfemiso04_;
	TBranch *pfemiso04_branch;
	bool pfemiso04_isLoaded;
	float	pfchiso04_;
	TBranch *pfchiso04_branch;
	bool pfchiso04_isLoaded;
	float	pfnhiso04_;
	TBranch *pfnhiso04_branch;
	bool pfnhiso04_isLoaded;
	float	radiso03_;
	TBranch *radiso03_branch;
	bool radiso03_isLoaded;
	float	radiso04_;
	TBranch *radiso04_branch;
	bool radiso04_isLoaded;
	float	iso2011_;
	TBranch *iso2011_branch;
	bool iso2011_isLoaded;
	float	ea04_;
	TBranch *ea04_branch;
	bool ea04_isLoaded;
	float	ea03_;
	TBranch *ea03_branch;
	bool ea03_isLoaded;
	float	dbeta03_;
	TBranch *dbeta03_branch;
	bool dbeta03_isLoaded;
	float	dbeta04_;
	TBranch *dbeta04_branch;
	bool dbeta04_isLoaded;
	float	el_test_pfchiso04_trkveto_;
	TBranch *el_test_pfchiso04_trkveto_branch;
	bool el_test_pfchiso04_trkveto_isLoaded;
	float	el_test_pfchiso04_dzcut_;
	TBranch *el_test_pfchiso04_dzcut_branch;
	bool el_test_pfchiso04_dzcut_isLoaded;
	float	el_test_pfchiso04_ebveto_;
	TBranch *el_test_pfchiso04_ebveto_branch;
	bool el_test_pfchiso04_ebveto_isLoaded;
	float	el_test_pfemiso04_ebveto_;
	TBranch *el_test_pfemiso04_ebveto_branch;
	bool el_test_pfemiso04_ebveto_isLoaded;
	float	eaem04_;
	TBranch *eaem04_branch;
	bool eaem04_isLoaded;
	float	eanh04_;
	TBranch *eanh04_branch;
	bool eanh04_isLoaded;
	float	gen_drs1_;
	TBranch *gen_drs1_branch;
	bool gen_drs1_isLoaded;
	float	gen_drs3_;
	TBranch *gen_drs3_branch;
	bool gen_drs3_isLoaded;
	unsigned int	HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_;
	TBranch *HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch;
	bool HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_isLoaded;
	unsigned int	HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_;
	TBranch *HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch;
	bool HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_isLoaded;
	unsigned int	HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_;
	TBranch *HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch;
	bool HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_isLoaded;
	unsigned int	HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_;
	TBranch *HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch;
	bool HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_isLoaded;
	unsigned int	HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_;
	TBranch *HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch;
	bool HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_isLoaded;
	unsigned int	HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_;
	TBranch *HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch;
	bool HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_isLoaded;
	unsigned int	HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_;
	TBranch *HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch;
	bool HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_isLoaded;
	unsigned int	HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_;
	TBranch *HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch;
	bool HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_isLoaded;
	unsigned int	HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_;
	TBranch *HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch;
	bool HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_isLoaded;
	unsigned int	HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_;
	TBranch *HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch;
	bool HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_isLoaded;
	unsigned int	HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_;
	TBranch *HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch;
	bool HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_isLoaded;
	unsigned int	HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_;
	TBranch *HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch;
	bool HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_isLoaded;
	unsigned int	HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_;
	TBranch *HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch;
	bool HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded;
	unsigned int	HLT_Mu17_Mu8_TrailingLeg_tag_;
	TBranch *HLT_Mu17_Mu8_TrailingLeg_tag_branch;
	bool HLT_Mu17_Mu8_TrailingLeg_tag_isLoaded;
	unsigned int	HLT_Mu17_Mu8_TrailingLeg_probe_;
	TBranch *HLT_Mu17_Mu8_TrailingLeg_probe_branch;
	bool HLT_Mu17_Mu8_TrailingLeg_probe_isLoaded;
	unsigned int	HLT_Mu17_Mu8_TrailingLeg_version_;
	TBranch *HLT_Mu17_Mu8_TrailingLeg_version_branch;
	bool HLT_Mu17_Mu8_TrailingLeg_version_isLoaded;
	unsigned int	HLT_Mu17_Mu8_LeadingLeg_tag_;
	TBranch *HLT_Mu17_Mu8_LeadingLeg_tag_branch;
	bool HLT_Mu17_Mu8_LeadingLeg_tag_isLoaded;
	unsigned int	HLT_Mu17_Mu8_LeadingLeg_probe_;
	TBranch *HLT_Mu17_Mu8_LeadingLeg_probe_branch;
	bool HLT_Mu17_Mu8_LeadingLeg_probe_isLoaded;
	unsigned int	HLT_Mu17_Mu8_LeadingLeg_version_;
	TBranch *HLT_Mu17_Mu8_LeadingLeg_version_branch;
	bool HLT_Mu17_Mu8_LeadingLeg_version_isLoaded;
	unsigned int	HLT_Mu17_Mu8_tag_;
	TBranch *HLT_Mu17_Mu8_tag_branch;
	bool HLT_Mu17_Mu8_tag_isLoaded;
	unsigned int	HLT_Mu17_Mu8_probe_;
	TBranch *HLT_Mu17_Mu8_probe_branch;
	bool HLT_Mu17_Mu8_probe_isLoaded;
	unsigned int	HLT_Mu17_Mu8_version_;
	TBranch *HLT_Mu17_Mu8_version_branch;
	bool HLT_Mu17_Mu8_version_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_TrailingLeg_tag_;
	TBranch *HLT_Mu17_TkMu8_TrailingLeg_tag_branch;
	bool HLT_Mu17_TkMu8_TrailingLeg_tag_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_TrailingLeg_probe_;
	TBranch *HLT_Mu17_TkMu8_TrailingLeg_probe_branch;
	bool HLT_Mu17_TkMu8_TrailingLeg_probe_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_TrailingLeg_version_;
	TBranch *HLT_Mu17_TkMu8_TrailingLeg_version_branch;
	bool HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_LeadingLeg_tag_;
	TBranch *HLT_Mu17_TkMu8_LeadingLeg_tag_branch;
	bool HLT_Mu17_TkMu8_LeadingLeg_tag_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_LeadingLeg_probe_;
	TBranch *HLT_Mu17_TkMu8_LeadingLeg_probe_branch;
	bool HLT_Mu17_TkMu8_LeadingLeg_probe_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_LeadingLeg_version_;
	TBranch *HLT_Mu17_TkMu8_LeadingLeg_version_branch;
	bool HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_;
	TBranch *HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch;
	bool HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_;
	TBranch *HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch;
	bool HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_;
	TBranch *HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch;
	bool HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_tag_;
	TBranch *HLT_Mu17_TkMu8_tag_branch;
	bool HLT_Mu17_TkMu8_tag_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_probe_;
	TBranch *HLT_Mu17_TkMu8_probe_branch;
	bool HLT_Mu17_TkMu8_probe_isLoaded;
	unsigned int	HLT_Mu17_TkMu8_version_;
	TBranch *HLT_Mu17_TkMu8_version_branch;
	bool HLT_Mu17_TkMu8_version_isLoaded;
	unsigned int	HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_;
	TBranch *HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch;
	bool HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_isLoaded;
	unsigned int	HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_;
	TBranch *HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch;
	bool HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_isLoaded;
	unsigned int	HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_;
	TBranch *HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch;
	bool HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded;
	unsigned int	HLT_IsoMu24_eta2p1_tag_;
	TBranch *HLT_IsoMu24_eta2p1_tag_branch;
	bool HLT_IsoMu24_eta2p1_tag_isLoaded;
	unsigned int	HLT_IsoMu24_eta2p1_probe_;
	TBranch *HLT_IsoMu24_eta2p1_probe_branch;
	bool HLT_IsoMu24_eta2p1_probe_isLoaded;
	unsigned int	HLT_IsoMu24_eta2p1_version_;
	TBranch *HLT_IsoMu24_eta2p1_version_branch;
	bool HLT_IsoMu24_eta2p1_version_isLoaded;
	unsigned int	HLT_Mu8_tag_;
	TBranch *HLT_Mu8_tag_branch;
	bool HLT_Mu8_tag_isLoaded;
	unsigned int	HLT_Mu8_probe_;
	TBranch *HLT_Mu8_probe_branch;
	bool HLT_Mu8_probe_isLoaded;
	unsigned int	HLT_Mu8_version_;
	TBranch *HLT_Mu8_version_branch;
	bool HLT_Mu8_version_isLoaded;
	unsigned int	HLT_Mu17_tag_;
	TBranch *HLT_Mu17_tag_branch;
	bool HLT_Mu17_tag_isLoaded;
	unsigned int	HLT_Mu17_probe_;
	TBranch *HLT_Mu17_probe_branch;
	bool HLT_Mu17_probe_isLoaded;
	unsigned int	HLT_Mu17_version_;
	TBranch *HLT_Mu17_version_branch;
	bool HLT_Mu17_version_isLoaded;
	unsigned int	HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_;
	TBranch *HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch;
	bool HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_isLoaded;
	unsigned int	HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_;
	TBranch *HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch;
	bool HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_isLoaded;
	unsigned int	HLT_Ele17_Ele8_L1sL1DoubleEG137_version_;
	TBranch *HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch;
	bool HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded;
	unsigned int	HLT_Ele17_Ele8_LeadingLeg_tag_;
	TBranch *HLT_Ele17_Ele8_LeadingLeg_tag_branch;
	bool HLT_Ele17_Ele8_LeadingLeg_tag_isLoaded;
	unsigned int	HLT_Ele17_Ele8_LeadingLeg_probe_;
	TBranch *HLT_Ele17_Ele8_LeadingLeg_probe_branch;
	bool HLT_Ele17_Ele8_LeadingLeg_probe_isLoaded;
	unsigned int	HLT_Ele17_Ele8_LeadingLeg_version_;
	TBranch *HLT_Ele17_Ele8_LeadingLeg_version_branch;
	bool HLT_Ele17_Ele8_LeadingLeg_version_isLoaded;
	unsigned int	HLT_Ele17_Ele8_TrailingLeg_tag_;
	TBranch *HLT_Ele17_Ele8_TrailingLeg_tag_branch;
	bool HLT_Ele17_Ele8_TrailingLeg_tag_isLoaded;
	unsigned int	HLT_Ele17_Ele8_TrailingLeg_probe_;
	TBranch *HLT_Ele17_Ele8_TrailingLeg_probe_branch;
	bool HLT_Ele17_Ele8_TrailingLeg_probe_isLoaded;
	unsigned int	HLT_Ele17_Ele8_TrailingLeg_version_;
	TBranch *HLT_Ele17_Ele8_TrailingLeg_version_branch;
	bool HLT_Ele17_Ele8_TrailingLeg_version_isLoaded;
	unsigned int	HLT_Ele17_Ele8_tag_;
	TBranch *HLT_Ele17_Ele8_tag_branch;
	bool HLT_Ele17_Ele8_tag_isLoaded;
	unsigned int	HLT_Ele17_Ele8_probe_;
	TBranch *HLT_Ele17_Ele8_probe_branch;
	bool HLT_Ele17_Ele8_probe_isLoaded;
	unsigned int	HLT_Ele17_Ele8_version_;
	TBranch *HLT_Ele17_Ele8_version_branch;
	bool HLT_Ele17_Ele8_version_isLoaded;
	unsigned int	HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_;
	TBranch *HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch;
	bool HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_isLoaded;
	unsigned int	HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_;
	TBranch *HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch;
	bool HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_isLoaded;
	unsigned int	HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_;
	TBranch *HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch;
	bool HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded;
	unsigned int	HLT_Ele27_WP80_tag_;
	TBranch *HLT_Ele27_WP80_tag_branch;
	bool HLT_Ele27_WP80_tag_isLoaded;
	unsigned int	HLT_Ele27_WP80_probe_;
	TBranch *HLT_Ele27_WP80_probe_branch;
	bool HLT_Ele27_WP80_probe_isLoaded;
	unsigned int	HLT_Ele27_WP80_version_;
	TBranch *HLT_Ele27_WP80_version_branch;
	bool HLT_Ele27_WP80_version_isLoaded;
	unsigned int	HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_;
	TBranch *HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch;
	bool HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_isLoaded;
	unsigned int	HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_;
	TBranch *HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch;
	bool HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_isLoaded;
	unsigned int	HLT_Ele17_Ele8_Mass50_LeadingLeg_version_;
	TBranch *HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch;
	bool HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded;
	unsigned int	HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_;
	TBranch *HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch;
	bool HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_isLoaded;
	unsigned int	HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_;
	TBranch *HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch;
	bool HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_isLoaded;
	unsigned int	HLT_Ele17_Ele8_Mass50_TrailingLeg_version_;
	TBranch *HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch;
	bool HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded;
	unsigned int	HLT_Ele20_SC4_Mass50_LeadingLeg_tag_;
	TBranch *HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch;
	bool HLT_Ele20_SC4_Mass50_LeadingLeg_tag_isLoaded;
	unsigned int	HLT_Ele20_SC4_Mass50_LeadingLeg_probe_;
	TBranch *HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch;
	bool HLT_Ele20_SC4_Mass50_LeadingLeg_probe_isLoaded;
	unsigned int	HLT_Ele20_SC4_Mass50_LeadingLeg_version_;
	TBranch *HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch;
	bool HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded;
	unsigned int	HLT_Ele20_SC4_Mass50_TrailingLeg_tag_;
	TBranch *HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch;
	bool HLT_Ele20_SC4_Mass50_TrailingLeg_tag_isLoaded;
	unsigned int	HLT_Ele20_SC4_Mass50_TrailingLeg_probe_;
	TBranch *HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch;
	bool HLT_Ele20_SC4_Mass50_TrailingLeg_probe_isLoaded;
	unsigned int	HLT_Ele20_SC4_Mass50_TrailingLeg_version_;
	TBranch *HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch;
	bool HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded;
	unsigned int	HLT_Ele32_SC17_Mass50_LeadingLeg_tag_;
	TBranch *HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch;
	bool HLT_Ele32_SC17_Mass50_LeadingLeg_tag_isLoaded;
	unsigned int	HLT_Ele32_SC17_Mass50_LeadingLeg_probe_;
	TBranch *HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch;
	bool HLT_Ele32_SC17_Mass50_LeadingLeg_probe_isLoaded;
	unsigned int	HLT_Ele32_SC17_Mass50_LeadingLeg_version_;
	TBranch *HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch;
	bool HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded;
	unsigned int	HLT_Ele32_SC17_Mass50_TrailingLeg_tag_;
	TBranch *HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch;
	bool HLT_Ele32_SC17_Mass50_TrailingLeg_tag_isLoaded;
	unsigned int	HLT_Ele32_SC17_Mass50_TrailingLeg_probe_;
	TBranch *HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch;
	bool HLT_Ele32_SC17_Mass50_TrailingLeg_probe_isLoaded;
	unsigned int	HLT_Ele32_SC17_Mass50_TrailingLeg_version_;
	TBranch *HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch;
	bool HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_;
	TBranch *HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch;
	bool HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_;
	TBranch *HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch;
	bool HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;
	TBranch *HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch;
	bool HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_;
	TBranch *HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch;
	bool HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_;
	TBranch *HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch;
	bool HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;
	TBranch *HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch;
	bool HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded;
	unsigned int	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_;
	TBranch *HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch;
	bool HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded;
	unsigned int	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_;
	TBranch *HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch;
	bool HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded;
	unsigned int	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;
	TBranch *HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch;
	bool HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded;
	unsigned int	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_;
	TBranch *HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch;
	bool HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded;
	unsigned int	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_;
	TBranch *HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch;
	bool HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded;
	unsigned int	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;
	TBranch *HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch;
	bool HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_TrkIdVL_tag_;
	TBranch *HLT_Ele8_CaloIdT_TrkIdVL_tag_branch;
	bool HLT_Ele8_CaloIdT_TrkIdVL_tag_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_TrkIdVL_probe_;
	TBranch *HLT_Ele8_CaloIdT_TrkIdVL_probe_branch;
	bool HLT_Ele8_CaloIdT_TrkIdVL_probe_isLoaded;
	unsigned int	HLT_Ele8_CaloIdT_TrkIdVL_version_;
	TBranch *HLT_Ele8_CaloIdT_TrkIdVL_version_branch;
	bool HLT_Ele8_CaloIdT_TrkIdVL_version_isLoaded;
public: 
void Init(TTree *tree) {
	tag_branch = 0;
	if (tree->GetBranch("tag") != 0) {
		tag_branch = tree->GetBranch("tag");
		if (tag_branch) {tag_branch->SetAddress(&tag_);}
	}
	probe_branch = 0;
	if (tree->GetBranch("probe") != 0) {
		probe_branch = tree->GetBranch("probe");
		if (probe_branch) {probe_branch->SetAddress(&probe_);}
	}
	jet1_branch = 0;
	if (tree->GetBranch("jet1") != 0) {
		jet1_branch = tree->GetBranch("jet1");
		if (jet1_branch) {jet1_branch->SetAddress(&jet1_);}
	}
	jet2_branch = 0;
	if (tree->GetBranch("jet2") != 0) {
		jet2_branch = tree->GetBranch("jet2");
		if (jet2_branch) {jet2_branch->SetAddress(&jet2_);}
	}
	jet3_branch = 0;
	if (tree->GetBranch("jet3") != 0) {
		jet3_branch = tree->GetBranch("jet3");
		if (jet3_branch) {jet3_branch->SetAddress(&jet3_);}
	}
  tree->SetMakeClass(1);
	event_branch = 0;
	if (tree->GetBranch("event") != 0) {
		event_branch = tree->GetBranch("event");
		if (event_branch) {event_branch->SetAddress(&event_);}
	}
	run_branch = 0;
	if (tree->GetBranch("run") != 0) {
		run_branch = tree->GetBranch("run");
		if (run_branch) {run_branch->SetAddress(&run_);}
	}
	lumi_branch = 0;
	if (tree->GetBranch("lumi") != 0) {
		lumi_branch = tree->GetBranch("lumi");
		if (lumi_branch) {lumi_branch->SetAddress(&lumi_);}
	}
	rnd_branch = 0;
	if (tree->GetBranch("rnd") != 0) {
		rnd_branch = tree->GetBranch("rnd");
		if (rnd_branch) {rnd_branch->SetAddress(&rnd_);}
	}
	nvtx_branch = 0;
	if (tree->GetBranch("nvtx") != 0) {
		nvtx_branch = tree->GetBranch("nvtx");
		if (nvtx_branch) {nvtx_branch->SetAddress(&nvtx_);}
	}
	npu_branch = 0;
	if (tree->GetBranch("npu") != 0) {
		npu_branch = tree->GetBranch("npu");
		if (npu_branch) {npu_branch->SetAddress(&npu_);}
	}
	npuPlusOne_branch = 0;
	if (tree->GetBranch("npuPlusOne") != 0) {
		npuPlusOne_branch = tree->GetBranch("npuPlusOne");
		if (npuPlusOne_branch) {npuPlusOne_branch->SetAddress(&npuPlusOne_);}
	}
	npuMinusOne_branch = 0;
	if (tree->GetBranch("npuMinusOne") != 0) {
		npuMinusOne_branch = tree->GetBranch("npuMinusOne");
		if (npuMinusOne_branch) {npuMinusOne_branch->SetAddress(&npuMinusOne_);}
	}
	leptonSelection_branch = 0;
	if (tree->GetBranch("leptonSelection") != 0) {
		leptonSelection_branch = tree->GetBranch("leptonSelection");
		if (leptonSelection_branch) {leptonSelection_branch->SetAddress(&leptonSelection_);}
	}
	eventSelection_branch = 0;
	if (tree->GetBranch("eventSelection") != 0) {
		eventSelection_branch = tree->GetBranch("eventSelection");
		if (eventSelection_branch) {eventSelection_branch->SetAddress(&eventSelection_);}
	}
	rhoIsoAll_branch = 0;
	if (tree->GetBranch("rhoIsoAll") != 0) {
		rhoIsoAll_branch = tree->GetBranch("rhoIsoAll");
		if (rhoIsoAll_branch) {rhoIsoAll_branch->SetAddress(&rhoIsoAll_);}
	}
	rhoIsoAllCentral_branch = 0;
	if (tree->GetBranch("rhoIsoAllCentral") != 0) {
		rhoIsoAllCentral_branch = tree->GetBranch("rhoIsoAllCentral");
		if (rhoIsoAllCentral_branch) {rhoIsoAllCentral_branch->SetAddress(&rhoIsoAllCentral_);}
	}
	rhoIsoNeutral_branch = 0;
	if (tree->GetBranch("rhoIsoNeutral") != 0) {
		rhoIsoNeutral_branch = tree->GetBranch("rhoIsoNeutral");
		if (rhoIsoNeutral_branch) {rhoIsoNeutral_branch->SetAddress(&rhoIsoNeutral_);}
	}
	tagAndProbeMass_branch = 0;
	if (tree->GetBranch("tagAndProbeMass") != 0) {
		tagAndProbeMass_branch = tree->GetBranch("tagAndProbeMass");
		if (tagAndProbeMass_branch) {tagAndProbeMass_branch->SetAddress(&tagAndProbeMass_);}
	}
	tagAndProbeIsRandom_branch = 0;
	if (tree->GetBranch("tagAndProbeIsRandom") != 0) {
		tagAndProbeIsRandom_branch = tree->GetBranch("tagAndProbeIsRandom");
		if (tagAndProbeIsRandom_branch) {tagAndProbeIsRandom_branch->SetAddress(&tagAndProbeIsRandom_);}
	}
	qTag_branch = 0;
	if (tree->GetBranch("qTag") != 0) {
		qTag_branch = tree->GetBranch("qTag");
		if (qTag_branch) {qTag_branch->SetAddress(&qTag_);}
	}
	qProbe_branch = 0;
	if (tree->GetBranch("qProbe") != 0) {
		qProbe_branch = tree->GetBranch("qProbe");
		if (qProbe_branch) {qProbe_branch->SetAddress(&qProbe_);}
	}
	scale1fb_branch = 0;
	if (tree->GetBranch("scale1fb") != 0) {
		scale1fb_branch = tree->GetBranch("scale1fb");
		if (scale1fb_branch) {scale1fb_branch->SetAddress(&scale1fb_);}
	}
	met_branch = 0;
	if (tree->GetBranch("met") != 0) {
		met_branch = tree->GetBranch("met");
		if (met_branch) {met_branch->SetAddress(&met_);}
	}
	metPhi_branch = 0;
	if (tree->GetBranch("metPhi") != 0) {
		metPhi_branch = tree->GetBranch("metPhi");
		if (metPhi_branch) {metPhi_branch->SetAddress(&metPhi_);}
	}
	trackMet_branch = 0;
	if (tree->GetBranch("trackMet") != 0) {
		trackMet_branch = tree->GetBranch("trackMet");
		if (trackMet_branch) {trackMet_branch->SetAddress(&trackMet_);}
	}
	trackMetPhi_branch = 0;
	if (tree->GetBranch("trackMetPhi") != 0) {
		trackMetPhi_branch = tree->GetBranch("trackMetPhi");
		if (trackMetPhi_branch) {trackMetPhi_branch->SetAddress(&trackMetPhi_);}
	}
	njets_branch = 0;
	if (tree->GetBranch("njets") != 0) {
		njets_branch = tree->GetBranch("njets");
		if (njets_branch) {njets_branch->SetAddress(&njets_);}
	}
	nleps_branch = 0;
	if (tree->GetBranch("nleps") != 0) {
		nleps_branch = tree->GetBranch("nleps");
		if (nleps_branch) {nleps_branch->SetAddress(&nleps_);}
	}
	hltPrescale_branch = 0;
	if (tree->GetBranch("hltPrescale") != 0) {
		hltPrescale_branch = tree->GetBranch("hltPrescale");
		if (hltPrescale_branch) {hltPrescale_branch->SetAddress(&hltPrescale_);}
	}
	sumet_branch = 0;
	if (tree->GetBranch("sumet") != 0) {
		sumet_branch = tree->GetBranch("sumet");
		if (sumet_branch) {sumet_branch->SetAddress(&sumet_);}
	}
	metSig_branch = 0;
	if (tree->GetBranch("metSig") != 0) {
		metSig_branch = tree->GetBranch("metSig");
		if (metSig_branch) {metSig_branch->SetAddress(&metSig_);}
	}
	mt_branch = 0;
	if (tree->GetBranch("mt") != 0) {
		mt_branch = tree->GetBranch("mt");
		if (mt_branch) {mt_branch->SetAddress(&mt_);}
	}
	dPhiProbeJet1_branch = 0;
	if (tree->GetBranch("dPhiProbeJet1") != 0) {
		dPhiProbeJet1_branch = tree->GetBranch("dPhiProbeJet1");
		if (dPhiProbeJet1_branch) {dPhiProbeJet1_branch->SetAddress(&dPhiProbeJet1_);}
	}
	chargedEmFracJet1_branch = 0;
	if (tree->GetBranch("chargedEmFracJet1") != 0) {
		chargedEmFracJet1_branch = tree->GetBranch("chargedEmFracJet1");
		if (chargedEmFracJet1_branch) {chargedEmFracJet1_branch->SetAddress(&chargedEmFracJet1_);}
	}
	neutralEmFracJet1_branch = 0;
	if (tree->GetBranch("neutralEmFracJet1") != 0) {
		neutralEmFracJet1_branch = tree->GetBranch("neutralEmFracJet1");
		if (neutralEmFracJet1_branch) {neutralEmFracJet1_branch->SetAddress(&neutralEmFracJet1_);}
	}
	chargedMuFracJet1_branch = 0;
	if (tree->GetBranch("chargedMuFracJet1") != 0) {
		chargedMuFracJet1_branch = tree->GetBranch("chargedMuFracJet1");
		if (chargedMuFracJet1_branch) {chargedMuFracJet1_branch->SetAddress(&chargedMuFracJet1_);}
	}
	electronHWW2011MVA_branch = 0;
	if (tree->GetBranch("electronHWW2011MVA") != 0) {
		electronHWW2011MVA_branch = tree->GetBranch("electronHWW2011MVA");
		if (electronHWW2011MVA_branch) {electronHWW2011MVA_branch->SetAddress(&electronHWW2011MVA_);}
	}
	egammaPOG2012MVA_branch = 0;
	if (tree->GetBranch("egammaPOG2012MVA") != 0) {
		egammaPOG2012MVA_branch = tree->GetBranch("egammaPOG2012MVA");
		if (egammaPOG2012MVA_branch) {egammaPOG2012MVA_branch->SetAddress(&egammaPOG2012MVA_);}
	}
	muonHZZ2012IsoRingsMVA_branch = 0;
	if (tree->GetBranch("muonHZZ2012IsoRingsMVA") != 0) {
		muonHZZ2012IsoRingsMVA_branch = tree->GetBranch("muonHZZ2012IsoRingsMVA");
		if (muonHZZ2012IsoRingsMVA_branch) {muonHZZ2012IsoRingsMVA_branch->SetAddress(&muonHZZ2012IsoRingsMVA_);}
	}
	vetoId_branch = 0;
	if (tree->GetBranch("vetoId") != 0) {
		vetoId_branch = tree->GetBranch("vetoId");
		if (vetoId_branch) {vetoId_branch->SetAddress(&vetoId_);}
	}
	looseId_branch = 0;
	if (tree->GetBranch("looseId") != 0) {
		looseId_branch = tree->GetBranch("looseId");
		if (looseId_branch) {looseId_branch->SetAddress(&looseId_);}
	}
	mediumId_branch = 0;
	if (tree->GetBranch("mediumId") != 0) {
		mediumId_branch = tree->GetBranch("mediumId");
		if (mediumId_branch) {mediumId_branch->SetAddress(&mediumId_);}
	}
	tightId_branch = 0;
	if (tree->GetBranch("tightId") != 0) {
		tightId_branch = tree->GetBranch("tightId");
		if (tightId_branch) {tightId_branch->SetAddress(&tightId_);}
	}
	pfmva_branch = 0;
	if (tree->GetBranch("pfmva") != 0) {
		pfmva_branch = tree->GetBranch("pfmva");
		if (pfmva_branch) {pfmva_branch->SetAddress(&pfmva_);}
	}
	sceta_branch = 0;
	if (tree->GetBranch("sceta") != 0) {
		sceta_branch = tree->GetBranch("sceta");
		if (sceta_branch) {sceta_branch->SetAddress(&sceta_);}
	}
	scphi_branch = 0;
	if (tree->GetBranch("scphi") != 0) {
		scphi_branch = tree->GetBranch("scphi");
		if (scphi_branch) {scphi_branch->SetAddress(&scphi_);}
	}
	scenergy_branch = 0;
	if (tree->GetBranch("scenergy") != 0) {
		scenergy_branch = tree->GetBranch("scenergy");
		if (scenergy_branch) {scenergy_branch->SetAddress(&scenergy_);}
	}
	chargesAgree_branch = 0;
	if (tree->GetBranch("chargesAgree") != 0) {
		chargesAgree_branch = tree->GetBranch("chargesAgree");
		if (chargesAgree_branch) {chargesAgree_branch->SetAddress(&chargesAgree_);}
	}
	eopin_branch = 0;
	if (tree->GetBranch("eopin") != 0) {
		eopin_branch = tree->GetBranch("eopin");
		if (eopin_branch) {eopin_branch->SetAddress(&eopin_);}
	}
	ooemoop_branch = 0;
	if (tree->GetBranch("ooemoop") != 0) {
		ooemoop_branch = tree->GetBranch("ooemoop");
		if (ooemoop_branch) {ooemoop_branch->SetAddress(&ooemoop_);}
	}
	fbrem_branch = 0;
	if (tree->GetBranch("fbrem") != 0) {
		fbrem_branch = tree->GetBranch("fbrem");
		if (fbrem_branch) {fbrem_branch->SetAddress(&fbrem_);}
	}
	detain_branch = 0;
	if (tree->GetBranch("detain") != 0) {
		detain_branch = tree->GetBranch("detain");
		if (detain_branch) {detain_branch->SetAddress(&detain_);}
	}
	dphiin_branch = 0;
	if (tree->GetBranch("dphiin") != 0) {
		dphiin_branch = tree->GetBranch("dphiin");
		if (dphiin_branch) {dphiin_branch->SetAddress(&dphiin_);}
	}
	hoe_branch = 0;
	if (tree->GetBranch("hoe") != 0) {
		hoe_branch = tree->GetBranch("hoe");
		if (hoe_branch) {hoe_branch->SetAddress(&hoe_);}
	}
	hoetow_branch = 0;
	if (tree->GetBranch("hoetow") != 0) {
		hoetow_branch = tree->GetBranch("hoetow");
		if (hoetow_branch) {hoetow_branch->SetAddress(&hoetow_);}
	}
	sieie_branch = 0;
	if (tree->GetBranch("sieie") != 0) {
		sieie_branch = tree->GetBranch("sieie");
		if (sieie_branch) {sieie_branch->SetAddress(&sieie_);}
	}
	d0vtx_branch = 0;
	if (tree->GetBranch("d0vtx") != 0) {
		d0vtx_branch = tree->GetBranch("d0vtx");
		if (d0vtx_branch) {d0vtx_branch->SetAddress(&d0vtx_);}
	}
	dzvtx_branch = 0;
	if (tree->GetBranch("dzvtx") != 0) {
		dzvtx_branch = tree->GetBranch("dzvtx");
		if (dzvtx_branch) {dzvtx_branch->SetAddress(&dzvtx_);}
	}
	vfitprob_branch = 0;
	if (tree->GetBranch("vfitprob") != 0) {
		vfitprob_branch = tree->GetBranch("vfitprob");
		if (vfitprob_branch) {vfitprob_branch->SetAddress(&vfitprob_);}
	}
	mhit_branch = 0;
	if (tree->GetBranch("mhit") != 0) {
		mhit_branch = tree->GetBranch("mhit");
		if (mhit_branch) {mhit_branch->SetAddress(&mhit_);}
	}
	ecaliso_branch = 0;
	if (tree->GetBranch("ecaliso") != 0) {
		ecaliso_branch = tree->GetBranch("ecaliso");
		if (ecaliso_branch) {ecaliso_branch->SetAddress(&ecaliso_);}
	}
	hcaliso_branch = 0;
	if (tree->GetBranch("hcaliso") != 0) {
		hcaliso_branch = tree->GetBranch("hcaliso");
		if (hcaliso_branch) {hcaliso_branch->SetAddress(&hcaliso_);}
	}
	trkiso_branch = 0;
	if (tree->GetBranch("trkiso") != 0) {
		trkiso_branch = tree->GetBranch("trkiso");
		if (trkiso_branch) {trkiso_branch->SetAddress(&trkiso_);}
	}
	pfemiso03_branch = 0;
	if (tree->GetBranch("pfemiso03") != 0) {
		pfemiso03_branch = tree->GetBranch("pfemiso03");
		if (pfemiso03_branch) {pfemiso03_branch->SetAddress(&pfemiso03_);}
	}
	pfchiso03_branch = 0;
	if (tree->GetBranch("pfchiso03") != 0) {
		pfchiso03_branch = tree->GetBranch("pfchiso03");
		if (pfchiso03_branch) {pfchiso03_branch->SetAddress(&pfchiso03_);}
	}
	pfnhiso03_branch = 0;
	if (tree->GetBranch("pfnhiso03") != 0) {
		pfnhiso03_branch = tree->GetBranch("pfnhiso03");
		if (pfnhiso03_branch) {pfnhiso03_branch->SetAddress(&pfnhiso03_);}
	}
	pfemiso04_branch = 0;
	if (tree->GetBranch("pfemiso04") != 0) {
		pfemiso04_branch = tree->GetBranch("pfemiso04");
		if (pfemiso04_branch) {pfemiso04_branch->SetAddress(&pfemiso04_);}
	}
	pfchiso04_branch = 0;
	if (tree->GetBranch("pfchiso04") != 0) {
		pfchiso04_branch = tree->GetBranch("pfchiso04");
		if (pfchiso04_branch) {pfchiso04_branch->SetAddress(&pfchiso04_);}
	}
	pfnhiso04_branch = 0;
	if (tree->GetBranch("pfnhiso04") != 0) {
		pfnhiso04_branch = tree->GetBranch("pfnhiso04");
		if (pfnhiso04_branch) {pfnhiso04_branch->SetAddress(&pfnhiso04_);}
	}
	radiso03_branch = 0;
	if (tree->GetBranch("radiso03") != 0) {
		radiso03_branch = tree->GetBranch("radiso03");
		if (radiso03_branch) {radiso03_branch->SetAddress(&radiso03_);}
	}
	radiso04_branch = 0;
	if (tree->GetBranch("radiso04") != 0) {
		radiso04_branch = tree->GetBranch("radiso04");
		if (radiso04_branch) {radiso04_branch->SetAddress(&radiso04_);}
	}
	iso2011_branch = 0;
	if (tree->GetBranch("iso2011") != 0) {
		iso2011_branch = tree->GetBranch("iso2011");
		if (iso2011_branch) {iso2011_branch->SetAddress(&iso2011_);}
	}
	ea04_branch = 0;
	if (tree->GetBranch("ea04") != 0) {
		ea04_branch = tree->GetBranch("ea04");
		if (ea04_branch) {ea04_branch->SetAddress(&ea04_);}
	}
	ea03_branch = 0;
	if (tree->GetBranch("ea03") != 0) {
		ea03_branch = tree->GetBranch("ea03");
		if (ea03_branch) {ea03_branch->SetAddress(&ea03_);}
	}
	dbeta03_branch = 0;
	if (tree->GetBranch("dbeta03") != 0) {
		dbeta03_branch = tree->GetBranch("dbeta03");
		if (dbeta03_branch) {dbeta03_branch->SetAddress(&dbeta03_);}
	}
	dbeta04_branch = 0;
	if (tree->GetBranch("dbeta04") != 0) {
		dbeta04_branch = tree->GetBranch("dbeta04");
		if (dbeta04_branch) {dbeta04_branch->SetAddress(&dbeta04_);}
	}
	el_test_pfchiso04_trkveto_branch = 0;
	if (tree->GetBranch("el_test_pfchiso04_trkveto") != 0) {
		el_test_pfchiso04_trkveto_branch = tree->GetBranch("el_test_pfchiso04_trkveto");
		if (el_test_pfchiso04_trkveto_branch) {el_test_pfchiso04_trkveto_branch->SetAddress(&el_test_pfchiso04_trkveto_);}
	}
	el_test_pfchiso04_dzcut_branch = 0;
	if (tree->GetBranch("el_test_pfchiso04_dzcut") != 0) {
		el_test_pfchiso04_dzcut_branch = tree->GetBranch("el_test_pfchiso04_dzcut");
		if (el_test_pfchiso04_dzcut_branch) {el_test_pfchiso04_dzcut_branch->SetAddress(&el_test_pfchiso04_dzcut_);}
	}
	el_test_pfchiso04_ebveto_branch = 0;
	if (tree->GetBranch("el_test_pfchiso04_ebveto") != 0) {
		el_test_pfchiso04_ebveto_branch = tree->GetBranch("el_test_pfchiso04_ebveto");
		if (el_test_pfchiso04_ebveto_branch) {el_test_pfchiso04_ebveto_branch->SetAddress(&el_test_pfchiso04_ebveto_);}
	}
	el_test_pfemiso04_ebveto_branch = 0;
	if (tree->GetBranch("el_test_pfemiso04_ebveto") != 0) {
		el_test_pfemiso04_ebveto_branch = tree->GetBranch("el_test_pfemiso04_ebveto");
		if (el_test_pfemiso04_ebveto_branch) {el_test_pfemiso04_ebveto_branch->SetAddress(&el_test_pfemiso04_ebveto_);}
	}
	eaem04_branch = 0;
	if (tree->GetBranch("eaem04") != 0) {
		eaem04_branch = tree->GetBranch("eaem04");
		if (eaem04_branch) {eaem04_branch->SetAddress(&eaem04_);}
	}
	eanh04_branch = 0;
	if (tree->GetBranch("eanh04") != 0) {
		eanh04_branch = tree->GetBranch("eanh04");
		if (eanh04_branch) {eanh04_branch->SetAddress(&eanh04_);}
	}
	gen_drs1_branch = 0;
	if (tree->GetBranch("gen_drs1") != 0) {
		gen_drs1_branch = tree->GetBranch("gen_drs1");
		if (gen_drs1_branch) {gen_drs1_branch->SetAddress(&gen_drs1_);}
	}
	gen_drs3_branch = 0;
	if (tree->GetBranch("gen_drs3") != 0) {
		gen_drs3_branch = tree->GetBranch("gen_drs3");
		if (gen_drs3_branch) {gen_drs3_branch->SetAddress(&gen_drs3_);}
	}
	HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch = 0;
	if (tree->GetBranch("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly") != 0) {
		HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch = tree->GetBranch("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly");
		if (HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch) {HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch->SetAddress(&HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_);}
	}
	HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch = 0;
	if (tree->GetBranch("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version") != 0) {
		HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch = tree->GetBranch("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version");
		if (HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch) {HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch->SetAddress(&HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_);}
	}
	HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch = 0;
	if (tree->GetBranch("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly") != 0) {
		HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch = tree->GetBranch("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly");
		if (HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch) {HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch->SetAddress(&HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_);}
	}
	HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch = 0;
	if (tree->GetBranch("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version") != 0) {
		HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch = tree->GetBranch("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version");
		if (HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch) {HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch->SetAddress(&HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_);}
	}
	HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch = 0;
	if (tree->GetBranch("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly") != 0) {
		HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch = tree->GetBranch("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly");
		if (HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch) {HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch->SetAddress(&HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_);}
	}
	HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch = 0;
	if (tree->GetBranch("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version") != 0) {
		HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch = tree->GetBranch("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version");
		if (HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch) {HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch->SetAddress(&HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_);}
	}
	HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch = 0;
	if (tree->GetBranch("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly") != 0) {
		HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch = tree->GetBranch("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly");
		if (HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch) {HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch->SetAddress(&HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_);}
	}
	HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch = 0;
	if (tree->GetBranch("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version") != 0) {
		HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch = tree->GetBranch("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version");
		if (HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch) {HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch->SetAddress(&HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_);}
	}
	HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch = 0;
	if (tree->GetBranch("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly") != 0) {
		HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch = tree->GetBranch("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly");
		if (HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch) {HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch->SetAddress(&HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_);}
	}
	HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch = 0;
	if (tree->GetBranch("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version") != 0) {
		HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch = tree->GetBranch("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version");
		if (HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch) {HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch->SetAddress(&HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_);}
	}
	HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag") != 0) {
		HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch = tree->GetBranch("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag");
		if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch) {HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch->SetAddress(&HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_);}
	}
	HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe") != 0) {
		HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch = tree->GetBranch("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe");
		if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch) {HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch->SetAddress(&HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_);}
	}
	HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version") != 0) {
		HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch = tree->GetBranch("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version");
		if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch) {HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch->SetAddress(&HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_);}
	}
	HLT_Mu17_Mu8_TrailingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_TrailingLeg_tag") != 0) {
		HLT_Mu17_Mu8_TrailingLeg_tag_branch = tree->GetBranch("HLT_Mu17_Mu8_TrailingLeg_tag");
		if (HLT_Mu17_Mu8_TrailingLeg_tag_branch) {HLT_Mu17_Mu8_TrailingLeg_tag_branch->SetAddress(&HLT_Mu17_Mu8_TrailingLeg_tag_);}
	}
	HLT_Mu17_Mu8_TrailingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_TrailingLeg_probe") != 0) {
		HLT_Mu17_Mu8_TrailingLeg_probe_branch = tree->GetBranch("HLT_Mu17_Mu8_TrailingLeg_probe");
		if (HLT_Mu17_Mu8_TrailingLeg_probe_branch) {HLT_Mu17_Mu8_TrailingLeg_probe_branch->SetAddress(&HLT_Mu17_Mu8_TrailingLeg_probe_);}
	}
	HLT_Mu17_Mu8_TrailingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_TrailingLeg_version") != 0) {
		HLT_Mu17_Mu8_TrailingLeg_version_branch = tree->GetBranch("HLT_Mu17_Mu8_TrailingLeg_version");
		if (HLT_Mu17_Mu8_TrailingLeg_version_branch) {HLT_Mu17_Mu8_TrailingLeg_version_branch->SetAddress(&HLT_Mu17_Mu8_TrailingLeg_version_);}
	}
	HLT_Mu17_Mu8_LeadingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_LeadingLeg_tag") != 0) {
		HLT_Mu17_Mu8_LeadingLeg_tag_branch = tree->GetBranch("HLT_Mu17_Mu8_LeadingLeg_tag");
		if (HLT_Mu17_Mu8_LeadingLeg_tag_branch) {HLT_Mu17_Mu8_LeadingLeg_tag_branch->SetAddress(&HLT_Mu17_Mu8_LeadingLeg_tag_);}
	}
	HLT_Mu17_Mu8_LeadingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_LeadingLeg_probe") != 0) {
		HLT_Mu17_Mu8_LeadingLeg_probe_branch = tree->GetBranch("HLT_Mu17_Mu8_LeadingLeg_probe");
		if (HLT_Mu17_Mu8_LeadingLeg_probe_branch) {HLT_Mu17_Mu8_LeadingLeg_probe_branch->SetAddress(&HLT_Mu17_Mu8_LeadingLeg_probe_);}
	}
	HLT_Mu17_Mu8_LeadingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_LeadingLeg_version") != 0) {
		HLT_Mu17_Mu8_LeadingLeg_version_branch = tree->GetBranch("HLT_Mu17_Mu8_LeadingLeg_version");
		if (HLT_Mu17_Mu8_LeadingLeg_version_branch) {HLT_Mu17_Mu8_LeadingLeg_version_branch->SetAddress(&HLT_Mu17_Mu8_LeadingLeg_version_);}
	}
	HLT_Mu17_Mu8_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_tag") != 0) {
		HLT_Mu17_Mu8_tag_branch = tree->GetBranch("HLT_Mu17_Mu8_tag");
		if (HLT_Mu17_Mu8_tag_branch) {HLT_Mu17_Mu8_tag_branch->SetAddress(&HLT_Mu17_Mu8_tag_);}
	}
	HLT_Mu17_Mu8_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_probe") != 0) {
		HLT_Mu17_Mu8_probe_branch = tree->GetBranch("HLT_Mu17_Mu8_probe");
		if (HLT_Mu17_Mu8_probe_branch) {HLT_Mu17_Mu8_probe_branch->SetAddress(&HLT_Mu17_Mu8_probe_);}
	}
	HLT_Mu17_Mu8_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_Mu8_version") != 0) {
		HLT_Mu17_Mu8_version_branch = tree->GetBranch("HLT_Mu17_Mu8_version");
		if (HLT_Mu17_Mu8_version_branch) {HLT_Mu17_Mu8_version_branch->SetAddress(&HLT_Mu17_Mu8_version_);}
	}
	HLT_Mu17_TkMu8_TrailingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_TrailingLeg_tag") != 0) {
		HLT_Mu17_TkMu8_TrailingLeg_tag_branch = tree->GetBranch("HLT_Mu17_TkMu8_TrailingLeg_tag");
		if (HLT_Mu17_TkMu8_TrailingLeg_tag_branch) {HLT_Mu17_TkMu8_TrailingLeg_tag_branch->SetAddress(&HLT_Mu17_TkMu8_TrailingLeg_tag_);}
	}
	HLT_Mu17_TkMu8_TrailingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_TrailingLeg_probe") != 0) {
		HLT_Mu17_TkMu8_TrailingLeg_probe_branch = tree->GetBranch("HLT_Mu17_TkMu8_TrailingLeg_probe");
		if (HLT_Mu17_TkMu8_TrailingLeg_probe_branch) {HLT_Mu17_TkMu8_TrailingLeg_probe_branch->SetAddress(&HLT_Mu17_TkMu8_TrailingLeg_probe_);}
	}
	HLT_Mu17_TkMu8_TrailingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_TrailingLeg_version") != 0) {
		HLT_Mu17_TkMu8_TrailingLeg_version_branch = tree->GetBranch("HLT_Mu17_TkMu8_TrailingLeg_version");
		if (HLT_Mu17_TkMu8_TrailingLeg_version_branch) {HLT_Mu17_TkMu8_TrailingLeg_version_branch->SetAddress(&HLT_Mu17_TkMu8_TrailingLeg_version_);}
	}
	HLT_Mu17_TkMu8_LeadingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_LeadingLeg_tag") != 0) {
		HLT_Mu17_TkMu8_LeadingLeg_tag_branch = tree->GetBranch("HLT_Mu17_TkMu8_LeadingLeg_tag");
		if (HLT_Mu17_TkMu8_LeadingLeg_tag_branch) {HLT_Mu17_TkMu8_LeadingLeg_tag_branch->SetAddress(&HLT_Mu17_TkMu8_LeadingLeg_tag_);}
	}
	HLT_Mu17_TkMu8_LeadingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_LeadingLeg_probe") != 0) {
		HLT_Mu17_TkMu8_LeadingLeg_probe_branch = tree->GetBranch("HLT_Mu17_TkMu8_LeadingLeg_probe");
		if (HLT_Mu17_TkMu8_LeadingLeg_probe_branch) {HLT_Mu17_TkMu8_LeadingLeg_probe_branch->SetAddress(&HLT_Mu17_TkMu8_LeadingLeg_probe_);}
	}
	HLT_Mu17_TkMu8_LeadingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_LeadingLeg_version") != 0) {
		HLT_Mu17_TkMu8_LeadingLeg_version_branch = tree->GetBranch("HLT_Mu17_TkMu8_LeadingLeg_version");
		if (HLT_Mu17_TkMu8_LeadingLeg_version_branch) {HLT_Mu17_TkMu8_LeadingLeg_version_branch->SetAddress(&HLT_Mu17_TkMu8_LeadingLeg_version_);}
	}
	HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag") != 0) {
		HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch = tree->GetBranch("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag");
		if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch) {HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch->SetAddress(&HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_);}
	}
	HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe") != 0) {
		HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch = tree->GetBranch("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe");
		if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch) {HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch->SetAddress(&HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_);}
	}
	HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version") != 0) {
		HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch = tree->GetBranch("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version");
		if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch) {HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch->SetAddress(&HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_);}
	}
	HLT_Mu17_TkMu8_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_tag") != 0) {
		HLT_Mu17_TkMu8_tag_branch = tree->GetBranch("HLT_Mu17_TkMu8_tag");
		if (HLT_Mu17_TkMu8_tag_branch) {HLT_Mu17_TkMu8_tag_branch->SetAddress(&HLT_Mu17_TkMu8_tag_);}
	}
	HLT_Mu17_TkMu8_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_probe") != 0) {
		HLT_Mu17_TkMu8_probe_branch = tree->GetBranch("HLT_Mu17_TkMu8_probe");
		if (HLT_Mu17_TkMu8_probe_branch) {HLT_Mu17_TkMu8_probe_branch->SetAddress(&HLT_Mu17_TkMu8_probe_);}
	}
	HLT_Mu17_TkMu8_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_TkMu8_version") != 0) {
		HLT_Mu17_TkMu8_version_branch = tree->GetBranch("HLT_Mu17_TkMu8_version");
		if (HLT_Mu17_TkMu8_version_branch) {HLT_Mu17_TkMu8_version_branch->SetAddress(&HLT_Mu17_TkMu8_version_);}
	}
	HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch = 0;
	if (tree->GetBranch("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag") != 0) {
		HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch = tree->GetBranch("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag");
		if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch) {HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch->SetAddress(&HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_);}
	}
	HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch = 0;
	if (tree->GetBranch("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe") != 0) {
		HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch = tree->GetBranch("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe");
		if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch) {HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch->SetAddress(&HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_);}
	}
	HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch = 0;
	if (tree->GetBranch("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version") != 0) {
		HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch = tree->GetBranch("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version");
		if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch) {HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch->SetAddress(&HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_);}
	}
	HLT_IsoMu24_eta2p1_tag_branch = 0;
	if (tree->GetBranch("HLT_IsoMu24_eta2p1_tag") != 0) {
		HLT_IsoMu24_eta2p1_tag_branch = tree->GetBranch("HLT_IsoMu24_eta2p1_tag");
		if (HLT_IsoMu24_eta2p1_tag_branch) {HLT_IsoMu24_eta2p1_tag_branch->SetAddress(&HLT_IsoMu24_eta2p1_tag_);}
	}
	HLT_IsoMu24_eta2p1_probe_branch = 0;
	if (tree->GetBranch("HLT_IsoMu24_eta2p1_probe") != 0) {
		HLT_IsoMu24_eta2p1_probe_branch = tree->GetBranch("HLT_IsoMu24_eta2p1_probe");
		if (HLT_IsoMu24_eta2p1_probe_branch) {HLT_IsoMu24_eta2p1_probe_branch->SetAddress(&HLT_IsoMu24_eta2p1_probe_);}
	}
	HLT_IsoMu24_eta2p1_version_branch = 0;
	if (tree->GetBranch("HLT_IsoMu24_eta2p1_version") != 0) {
		HLT_IsoMu24_eta2p1_version_branch = tree->GetBranch("HLT_IsoMu24_eta2p1_version");
		if (HLT_IsoMu24_eta2p1_version_branch) {HLT_IsoMu24_eta2p1_version_branch->SetAddress(&HLT_IsoMu24_eta2p1_version_);}
	}
	HLT_Mu8_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu8_tag") != 0) {
		HLT_Mu8_tag_branch = tree->GetBranch("HLT_Mu8_tag");
		if (HLT_Mu8_tag_branch) {HLT_Mu8_tag_branch->SetAddress(&HLT_Mu8_tag_);}
	}
	HLT_Mu8_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu8_probe") != 0) {
		HLT_Mu8_probe_branch = tree->GetBranch("HLT_Mu8_probe");
		if (HLT_Mu8_probe_branch) {HLT_Mu8_probe_branch->SetAddress(&HLT_Mu8_probe_);}
	}
	HLT_Mu8_version_branch = 0;
	if (tree->GetBranch("HLT_Mu8_version") != 0) {
		HLT_Mu8_version_branch = tree->GetBranch("HLT_Mu8_version");
		if (HLT_Mu8_version_branch) {HLT_Mu8_version_branch->SetAddress(&HLT_Mu8_version_);}
	}
	HLT_Mu17_tag_branch = 0;
	if (tree->GetBranch("HLT_Mu17_tag") != 0) {
		HLT_Mu17_tag_branch = tree->GetBranch("HLT_Mu17_tag");
		if (HLT_Mu17_tag_branch) {HLT_Mu17_tag_branch->SetAddress(&HLT_Mu17_tag_);}
	}
	HLT_Mu17_probe_branch = 0;
	if (tree->GetBranch("HLT_Mu17_probe") != 0) {
		HLT_Mu17_probe_branch = tree->GetBranch("HLT_Mu17_probe");
		if (HLT_Mu17_probe_branch) {HLT_Mu17_probe_branch->SetAddress(&HLT_Mu17_probe_);}
	}
	HLT_Mu17_version_branch = 0;
	if (tree->GetBranch("HLT_Mu17_version") != 0) {
		HLT_Mu17_version_branch = tree->GetBranch("HLT_Mu17_version");
		if (HLT_Mu17_version_branch) {HLT_Mu17_version_branch->SetAddress(&HLT_Mu17_version_);}
	}
	HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_L1sL1DoubleEG137_tag") != 0) {
		HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch = tree->GetBranch("HLT_Ele17_Ele8_L1sL1DoubleEG137_tag");
		if (HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch) {HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch->SetAddress(&HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_);}
	}
	HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_L1sL1DoubleEG137_probe") != 0) {
		HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch = tree->GetBranch("HLT_Ele17_Ele8_L1sL1DoubleEG137_probe");
		if (HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch) {HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch->SetAddress(&HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_);}
	}
	HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_L1sL1DoubleEG137_version") != 0) {
		HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch = tree->GetBranch("HLT_Ele17_Ele8_L1sL1DoubleEG137_version");
		if (HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch) {HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch->SetAddress(&HLT_Ele17_Ele8_L1sL1DoubleEG137_version_);}
	}
	HLT_Ele17_Ele8_LeadingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_LeadingLeg_tag") != 0) {
		HLT_Ele17_Ele8_LeadingLeg_tag_branch = tree->GetBranch("HLT_Ele17_Ele8_LeadingLeg_tag");
		if (HLT_Ele17_Ele8_LeadingLeg_tag_branch) {HLT_Ele17_Ele8_LeadingLeg_tag_branch->SetAddress(&HLT_Ele17_Ele8_LeadingLeg_tag_);}
	}
	HLT_Ele17_Ele8_LeadingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_LeadingLeg_probe") != 0) {
		HLT_Ele17_Ele8_LeadingLeg_probe_branch = tree->GetBranch("HLT_Ele17_Ele8_LeadingLeg_probe");
		if (HLT_Ele17_Ele8_LeadingLeg_probe_branch) {HLT_Ele17_Ele8_LeadingLeg_probe_branch->SetAddress(&HLT_Ele17_Ele8_LeadingLeg_probe_);}
	}
	HLT_Ele17_Ele8_LeadingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_LeadingLeg_version") != 0) {
		HLT_Ele17_Ele8_LeadingLeg_version_branch = tree->GetBranch("HLT_Ele17_Ele8_LeadingLeg_version");
		if (HLT_Ele17_Ele8_LeadingLeg_version_branch) {HLT_Ele17_Ele8_LeadingLeg_version_branch->SetAddress(&HLT_Ele17_Ele8_LeadingLeg_version_);}
	}
	HLT_Ele17_Ele8_TrailingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_TrailingLeg_tag") != 0) {
		HLT_Ele17_Ele8_TrailingLeg_tag_branch = tree->GetBranch("HLT_Ele17_Ele8_TrailingLeg_tag");
		if (HLT_Ele17_Ele8_TrailingLeg_tag_branch) {HLT_Ele17_Ele8_TrailingLeg_tag_branch->SetAddress(&HLT_Ele17_Ele8_TrailingLeg_tag_);}
	}
	HLT_Ele17_Ele8_TrailingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_TrailingLeg_probe") != 0) {
		HLT_Ele17_Ele8_TrailingLeg_probe_branch = tree->GetBranch("HLT_Ele17_Ele8_TrailingLeg_probe");
		if (HLT_Ele17_Ele8_TrailingLeg_probe_branch) {HLT_Ele17_Ele8_TrailingLeg_probe_branch->SetAddress(&HLT_Ele17_Ele8_TrailingLeg_probe_);}
	}
	HLT_Ele17_Ele8_TrailingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_TrailingLeg_version") != 0) {
		HLT_Ele17_Ele8_TrailingLeg_version_branch = tree->GetBranch("HLT_Ele17_Ele8_TrailingLeg_version");
		if (HLT_Ele17_Ele8_TrailingLeg_version_branch) {HLT_Ele17_Ele8_TrailingLeg_version_branch->SetAddress(&HLT_Ele17_Ele8_TrailingLeg_version_);}
	}
	HLT_Ele17_Ele8_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_tag") != 0) {
		HLT_Ele17_Ele8_tag_branch = tree->GetBranch("HLT_Ele17_Ele8_tag");
		if (HLT_Ele17_Ele8_tag_branch) {HLT_Ele17_Ele8_tag_branch->SetAddress(&HLT_Ele17_Ele8_tag_);}
	}
	HLT_Ele17_Ele8_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_probe") != 0) {
		HLT_Ele17_Ele8_probe_branch = tree->GetBranch("HLT_Ele17_Ele8_probe");
		if (HLT_Ele17_Ele8_probe_branch) {HLT_Ele17_Ele8_probe_branch->SetAddress(&HLT_Ele17_Ele8_probe_);}
	}
	HLT_Ele17_Ele8_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_version") != 0) {
		HLT_Ele17_Ele8_version_branch = tree->GetBranch("HLT_Ele17_Ele8_version");
		if (HLT_Ele17_Ele8_version_branch) {HLT_Ele17_Ele8_version_branch->SetAddress(&HLT_Ele17_Ele8_version_);}
	}
	HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag") != 0) {
		HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch = tree->GetBranch("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag");
		if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch) {HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch->SetAddress(&HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_);}
	}
	HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe") != 0) {
		HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch = tree->GetBranch("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe");
		if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch) {HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch->SetAddress(&HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_);}
	}
	HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch = 0;
	if (tree->GetBranch("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version") != 0) {
		HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch = tree->GetBranch("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version");
		if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch) {HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch->SetAddress(&HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_);}
	}
	HLT_Ele27_WP80_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele27_WP80_tag") != 0) {
		HLT_Ele27_WP80_tag_branch = tree->GetBranch("HLT_Ele27_WP80_tag");
		if (HLT_Ele27_WP80_tag_branch) {HLT_Ele27_WP80_tag_branch->SetAddress(&HLT_Ele27_WP80_tag_);}
	}
	HLT_Ele27_WP80_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele27_WP80_probe") != 0) {
		HLT_Ele27_WP80_probe_branch = tree->GetBranch("HLT_Ele27_WP80_probe");
		if (HLT_Ele27_WP80_probe_branch) {HLT_Ele27_WP80_probe_branch->SetAddress(&HLT_Ele27_WP80_probe_);}
	}
	HLT_Ele27_WP80_version_branch = 0;
	if (tree->GetBranch("HLT_Ele27_WP80_version") != 0) {
		HLT_Ele27_WP80_version_branch = tree->GetBranch("HLT_Ele27_WP80_version");
		if (HLT_Ele27_WP80_version_branch) {HLT_Ele27_WP80_version_branch->SetAddress(&HLT_Ele27_WP80_version_);}
	}
	HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_Mass50_LeadingLeg_tag") != 0) {
		HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch = tree->GetBranch("HLT_Ele17_Ele8_Mass50_LeadingLeg_tag");
		if (HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch) {HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch->SetAddress(&HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_);}
	}
	HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_Mass50_LeadingLeg_probe") != 0) {
		HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch = tree->GetBranch("HLT_Ele17_Ele8_Mass50_LeadingLeg_probe");
		if (HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch) {HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch->SetAddress(&HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_);}
	}
	HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_Mass50_LeadingLeg_version") != 0) {
		HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch = tree->GetBranch("HLT_Ele17_Ele8_Mass50_LeadingLeg_version");
		if (HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch) {HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch->SetAddress(&HLT_Ele17_Ele8_Mass50_LeadingLeg_version_);}
	}
	HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_Mass50_TrailingLeg_tag") != 0) {
		HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch = tree->GetBranch("HLT_Ele17_Ele8_Mass50_TrailingLeg_tag");
		if (HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch) {HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch->SetAddress(&HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_);}
	}
	HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_Mass50_TrailingLeg_probe") != 0) {
		HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch = tree->GetBranch("HLT_Ele17_Ele8_Mass50_TrailingLeg_probe");
		if (HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch) {HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch->SetAddress(&HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_);}
	}
	HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_Ele8_Mass50_TrailingLeg_version") != 0) {
		HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch = tree->GetBranch("HLT_Ele17_Ele8_Mass50_TrailingLeg_version");
		if (HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch) {HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch->SetAddress(&HLT_Ele17_Ele8_Mass50_TrailingLeg_version_);}
	}
	HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele20_SC4_Mass50_LeadingLeg_tag") != 0) {
		HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch = tree->GetBranch("HLT_Ele20_SC4_Mass50_LeadingLeg_tag");
		if (HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch) {HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch->SetAddress(&HLT_Ele20_SC4_Mass50_LeadingLeg_tag_);}
	}
	HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele20_SC4_Mass50_LeadingLeg_probe") != 0) {
		HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch = tree->GetBranch("HLT_Ele20_SC4_Mass50_LeadingLeg_probe");
		if (HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch) {HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch->SetAddress(&HLT_Ele20_SC4_Mass50_LeadingLeg_probe_);}
	}
	HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele20_SC4_Mass50_LeadingLeg_version") != 0) {
		HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch = tree->GetBranch("HLT_Ele20_SC4_Mass50_LeadingLeg_version");
		if (HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch) {HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch->SetAddress(&HLT_Ele20_SC4_Mass50_LeadingLeg_version_);}
	}
	HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele20_SC4_Mass50_TrailingLeg_tag") != 0) {
		HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch = tree->GetBranch("HLT_Ele20_SC4_Mass50_TrailingLeg_tag");
		if (HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch) {HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch->SetAddress(&HLT_Ele20_SC4_Mass50_TrailingLeg_tag_);}
	}
	HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele20_SC4_Mass50_TrailingLeg_probe") != 0) {
		HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch = tree->GetBranch("HLT_Ele20_SC4_Mass50_TrailingLeg_probe");
		if (HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch) {HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch->SetAddress(&HLT_Ele20_SC4_Mass50_TrailingLeg_probe_);}
	}
	HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele20_SC4_Mass50_TrailingLeg_version") != 0) {
		HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch = tree->GetBranch("HLT_Ele20_SC4_Mass50_TrailingLeg_version");
		if (HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch) {HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch->SetAddress(&HLT_Ele20_SC4_Mass50_TrailingLeg_version_);}
	}
	HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele32_SC17_Mass50_LeadingLeg_tag") != 0) {
		HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch = tree->GetBranch("HLT_Ele32_SC17_Mass50_LeadingLeg_tag");
		if (HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch) {HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch->SetAddress(&HLT_Ele32_SC17_Mass50_LeadingLeg_tag_);}
	}
	HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele32_SC17_Mass50_LeadingLeg_probe") != 0) {
		HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch = tree->GetBranch("HLT_Ele32_SC17_Mass50_LeadingLeg_probe");
		if (HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch) {HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch->SetAddress(&HLT_Ele32_SC17_Mass50_LeadingLeg_probe_);}
	}
	HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele32_SC17_Mass50_LeadingLeg_version") != 0) {
		HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch = tree->GetBranch("HLT_Ele32_SC17_Mass50_LeadingLeg_version");
		if (HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch) {HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch->SetAddress(&HLT_Ele32_SC17_Mass50_LeadingLeg_version_);}
	}
	HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele32_SC17_Mass50_TrailingLeg_tag") != 0) {
		HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch = tree->GetBranch("HLT_Ele32_SC17_Mass50_TrailingLeg_tag");
		if (HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch) {HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch->SetAddress(&HLT_Ele32_SC17_Mass50_TrailingLeg_tag_);}
	}
	HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele32_SC17_Mass50_TrailingLeg_probe") != 0) {
		HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch = tree->GetBranch("HLT_Ele32_SC17_Mass50_TrailingLeg_probe");
		if (HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch) {HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch->SetAddress(&HLT_Ele32_SC17_Mass50_TrailingLeg_probe_);}
	}
	HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch = 0;
	if (tree->GetBranch("HLT_Ele32_SC17_Mass50_TrailingLeg_version") != 0) {
		HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch = tree->GetBranch("HLT_Ele32_SC17_Mass50_TrailingLeg_version");
		if (HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch) {HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch->SetAddress(&HLT_Ele32_SC17_Mass50_TrailingLeg_version_);}
	}
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag") != 0) {
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch = tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag");
		if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch) {HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch->SetAddress(&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_);}
	}
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe") != 0) {
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch = tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe");
		if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch) {HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch->SetAddress(&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_);}
	}
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version") != 0) {
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch = tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version");
		if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch) {HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch->SetAddress(&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_);}
	}
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag") != 0) {
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch = tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag");
		if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch) {HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch->SetAddress(&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_);}
	}
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe") != 0) {
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch = tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe");
		if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch) {HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch->SetAddress(&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_);}
	}
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version") != 0) {
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch = tree->GetBranch("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version");
		if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch) {HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch->SetAddress(&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_);}
	}
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag") != 0) {
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch = tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag");
		if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch) {HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch->SetAddress(&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_);}
	}
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe") != 0) {
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch = tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe");
		if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch) {HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch->SetAddress(&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_);}
	}
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version") != 0) {
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch = tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version");
		if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch) {HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch->SetAddress(&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_);}
	}
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag") != 0) {
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch = tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag");
		if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch) {HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch->SetAddress(&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_);}
	}
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe") != 0) {
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch = tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe");
		if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch) {HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch->SetAddress(&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_);}
	}
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch = 0;
	if (tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version") != 0) {
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch = tree->GetBranch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version");
		if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch) {HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch->SetAddress(&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_);}
	}
	HLT_Ele8_CaloIdT_TrkIdVL_tag_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_TrkIdVL_tag") != 0) {
		HLT_Ele8_CaloIdT_TrkIdVL_tag_branch = tree->GetBranch("HLT_Ele8_CaloIdT_TrkIdVL_tag");
		if (HLT_Ele8_CaloIdT_TrkIdVL_tag_branch) {HLT_Ele8_CaloIdT_TrkIdVL_tag_branch->SetAddress(&HLT_Ele8_CaloIdT_TrkIdVL_tag_);}
	}
	HLT_Ele8_CaloIdT_TrkIdVL_probe_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_TrkIdVL_probe") != 0) {
		HLT_Ele8_CaloIdT_TrkIdVL_probe_branch = tree->GetBranch("HLT_Ele8_CaloIdT_TrkIdVL_probe");
		if (HLT_Ele8_CaloIdT_TrkIdVL_probe_branch) {HLT_Ele8_CaloIdT_TrkIdVL_probe_branch->SetAddress(&HLT_Ele8_CaloIdT_TrkIdVL_probe_);}
	}
	HLT_Ele8_CaloIdT_TrkIdVL_version_branch = 0;
	if (tree->GetBranch("HLT_Ele8_CaloIdT_TrkIdVL_version") != 0) {
		HLT_Ele8_CaloIdT_TrkIdVL_version_branch = tree->GetBranch("HLT_Ele8_CaloIdT_TrkIdVL_version");
		if (HLT_Ele8_CaloIdT_TrkIdVL_version_branch) {HLT_Ele8_CaloIdT_TrkIdVL_version_branch->SetAddress(&HLT_Ele8_CaloIdT_TrkIdVL_version_);}
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		event_isLoaded = false;
		run_isLoaded = false;
		lumi_isLoaded = false;
		rnd_isLoaded = false;
		nvtx_isLoaded = false;
		npu_isLoaded = false;
		npuPlusOne_isLoaded = false;
		npuMinusOne_isLoaded = false;
		leptonSelection_isLoaded = false;
		eventSelection_isLoaded = false;
		rhoIsoAll_isLoaded = false;
		rhoIsoAllCentral_isLoaded = false;
		rhoIsoNeutral_isLoaded = false;
		tagAndProbeMass_isLoaded = false;
		tagAndProbeIsRandom_isLoaded = false;
		tag_isLoaded = false;
		probe_isLoaded = false;
		qTag_isLoaded = false;
		qProbe_isLoaded = false;
		scale1fb_isLoaded = false;
		jet1_isLoaded = false;
		jet2_isLoaded = false;
		jet3_isLoaded = false;
		met_isLoaded = false;
		metPhi_isLoaded = false;
		trackMet_isLoaded = false;
		trackMetPhi_isLoaded = false;
		njets_isLoaded = false;
		nleps_isLoaded = false;
		hltPrescale_isLoaded = false;
		sumet_isLoaded = false;
		metSig_isLoaded = false;
		mt_isLoaded = false;
		dPhiProbeJet1_isLoaded = false;
		chargedEmFracJet1_isLoaded = false;
		neutralEmFracJet1_isLoaded = false;
		chargedMuFracJet1_isLoaded = false;
		electronHWW2011MVA_isLoaded = false;
		egammaPOG2012MVA_isLoaded = false;
		muonHZZ2012IsoRingsMVA_isLoaded = false;
		vetoId_isLoaded = false;
		looseId_isLoaded = false;
		mediumId_isLoaded = false;
		tightId_isLoaded = false;
		pfmva_isLoaded = false;
		sceta_isLoaded = false;
		scphi_isLoaded = false;
		scenergy_isLoaded = false;
		chargesAgree_isLoaded = false;
		eopin_isLoaded = false;
		ooemoop_isLoaded = false;
		fbrem_isLoaded = false;
		detain_isLoaded = false;
		dphiin_isLoaded = false;
		hoe_isLoaded = false;
		hoetow_isLoaded = false;
		sieie_isLoaded = false;
		d0vtx_isLoaded = false;
		dzvtx_isLoaded = false;
		vfitprob_isLoaded = false;
		mhit_isLoaded = false;
		ecaliso_isLoaded = false;
		hcaliso_isLoaded = false;
		trkiso_isLoaded = false;
		pfemiso03_isLoaded = false;
		pfchiso03_isLoaded = false;
		pfnhiso03_isLoaded = false;
		pfemiso04_isLoaded = false;
		pfchiso04_isLoaded = false;
		pfnhiso04_isLoaded = false;
		radiso03_isLoaded = false;
		radiso04_isLoaded = false;
		iso2011_isLoaded = false;
		ea04_isLoaded = false;
		ea03_isLoaded = false;
		dbeta03_isLoaded = false;
		dbeta04_isLoaded = false;
		el_test_pfchiso04_trkveto_isLoaded = false;
		el_test_pfchiso04_dzcut_isLoaded = false;
		el_test_pfchiso04_ebveto_isLoaded = false;
		el_test_pfemiso04_ebveto_isLoaded = false;
		eaem04_isLoaded = false;
		eanh04_isLoaded = false;
		gen_drs1_isLoaded = false;
		gen_drs3_isLoaded = false;
		HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_isLoaded = false;
		HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = false;
		HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_isLoaded = false;
		HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = false;
		HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_isLoaded = false;
		HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = false;
		HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_isLoaded = false;
		HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = false;
		HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_isLoaded = false;
		HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = false;
		HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_isLoaded = false;
		HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_isLoaded = false;
		HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded = false;
		HLT_Mu17_Mu8_TrailingLeg_tag_isLoaded = false;
		HLT_Mu17_Mu8_TrailingLeg_probe_isLoaded = false;
		HLT_Mu17_Mu8_TrailingLeg_version_isLoaded = false;
		HLT_Mu17_Mu8_LeadingLeg_tag_isLoaded = false;
		HLT_Mu17_Mu8_LeadingLeg_probe_isLoaded = false;
		HLT_Mu17_Mu8_LeadingLeg_version_isLoaded = false;
		HLT_Mu17_Mu8_tag_isLoaded = false;
		HLT_Mu17_Mu8_probe_isLoaded = false;
		HLT_Mu17_Mu8_version_isLoaded = false;
		HLT_Mu17_TkMu8_TrailingLeg_tag_isLoaded = false;
		HLT_Mu17_TkMu8_TrailingLeg_probe_isLoaded = false;
		HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded = false;
		HLT_Mu17_TkMu8_LeadingLeg_tag_isLoaded = false;
		HLT_Mu17_TkMu8_LeadingLeg_probe_isLoaded = false;
		HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded = false;
		HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_isLoaded = false;
		HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_isLoaded = false;
		HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded = false;
		HLT_Mu17_TkMu8_tag_isLoaded = false;
		HLT_Mu17_TkMu8_probe_isLoaded = false;
		HLT_Mu17_TkMu8_version_isLoaded = false;
		HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_isLoaded = false;
		HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_isLoaded = false;
		HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded = false;
		HLT_IsoMu24_eta2p1_tag_isLoaded = false;
		HLT_IsoMu24_eta2p1_probe_isLoaded = false;
		HLT_IsoMu24_eta2p1_version_isLoaded = false;
		HLT_Mu8_tag_isLoaded = false;
		HLT_Mu8_probe_isLoaded = false;
		HLT_Mu8_version_isLoaded = false;
		HLT_Mu17_tag_isLoaded = false;
		HLT_Mu17_probe_isLoaded = false;
		HLT_Mu17_version_isLoaded = false;
		HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_isLoaded = false;
		HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_isLoaded = false;
		HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded = false;
		HLT_Ele17_Ele8_LeadingLeg_tag_isLoaded = false;
		HLT_Ele17_Ele8_LeadingLeg_probe_isLoaded = false;
		HLT_Ele17_Ele8_LeadingLeg_version_isLoaded = false;
		HLT_Ele17_Ele8_TrailingLeg_tag_isLoaded = false;
		HLT_Ele17_Ele8_TrailingLeg_probe_isLoaded = false;
		HLT_Ele17_Ele8_TrailingLeg_version_isLoaded = false;
		HLT_Ele17_Ele8_tag_isLoaded = false;
		HLT_Ele17_Ele8_probe_isLoaded = false;
		HLT_Ele17_Ele8_version_isLoaded = false;
		HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_isLoaded = false;
		HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_isLoaded = false;
		HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded = false;
		HLT_Ele27_WP80_tag_isLoaded = false;
		HLT_Ele27_WP80_probe_isLoaded = false;
		HLT_Ele27_WP80_version_isLoaded = false;
		HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_isLoaded = false;
		HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_isLoaded = false;
		HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded = false;
		HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_isLoaded = false;
		HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_isLoaded = false;
		HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded = false;
		HLT_Ele20_SC4_Mass50_LeadingLeg_tag_isLoaded = false;
		HLT_Ele20_SC4_Mass50_LeadingLeg_probe_isLoaded = false;
		HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded = false;
		HLT_Ele20_SC4_Mass50_TrailingLeg_tag_isLoaded = false;
		HLT_Ele20_SC4_Mass50_TrailingLeg_probe_isLoaded = false;
		HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded = false;
		HLT_Ele32_SC17_Mass50_LeadingLeg_tag_isLoaded = false;
		HLT_Ele32_SC17_Mass50_LeadingLeg_probe_isLoaded = false;
		HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded = false;
		HLT_Ele32_SC17_Mass50_TrailingLeg_tag_isLoaded = false;
		HLT_Ele32_SC17_Mass50_TrailingLeg_probe_isLoaded = false;
		HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded = false;
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded = false;
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded = false;
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded = false;
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded = false;
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded = false;
		HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded = false;
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded = false;
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded = false;
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded = false;
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded = false;
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded = false;
		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded = false;
		HLT_Ele8_CaloIdT_TrkIdVL_tag_isLoaded = false;
		HLT_Ele8_CaloIdT_TrkIdVL_probe_isLoaded = false;
		HLT_Ele8_CaloIdT_TrkIdVL_version_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (event_branch != 0) event();
	if (run_branch != 0) run();
	if (lumi_branch != 0) lumi();
	if (rnd_branch != 0) rnd();
	if (nvtx_branch != 0) nvtx();
	if (npu_branch != 0) npu();
	if (npuPlusOne_branch != 0) npuPlusOne();
	if (npuMinusOne_branch != 0) npuMinusOne();
	if (leptonSelection_branch != 0) leptonSelection();
	if (eventSelection_branch != 0) eventSelection();
	if (rhoIsoAll_branch != 0) rhoIsoAll();
	if (rhoIsoAllCentral_branch != 0) rhoIsoAllCentral();
	if (rhoIsoNeutral_branch != 0) rhoIsoNeutral();
	if (tagAndProbeMass_branch != 0) tagAndProbeMass();
	if (tagAndProbeIsRandom_branch != 0) tagAndProbeIsRandom();
	if (tag_branch != 0) tag();
	if (probe_branch != 0) probe();
	if (qTag_branch != 0) qTag();
	if (qProbe_branch != 0) qProbe();
	if (scale1fb_branch != 0) scale1fb();
	if (jet1_branch != 0) jet1();
	if (jet2_branch != 0) jet2();
	if (jet3_branch != 0) jet3();
	if (met_branch != 0) met();
	if (metPhi_branch != 0) metPhi();
	if (trackMet_branch != 0) trackMet();
	if (trackMetPhi_branch != 0) trackMetPhi();
	if (njets_branch != 0) njets();
	if (nleps_branch != 0) nleps();
	if (hltPrescale_branch != 0) hltPrescale();
	if (sumet_branch != 0) sumet();
	if (metSig_branch != 0) metSig();
	if (mt_branch != 0) mt();
	if (dPhiProbeJet1_branch != 0) dPhiProbeJet1();
	if (chargedEmFracJet1_branch != 0) chargedEmFracJet1();
	if (neutralEmFracJet1_branch != 0) neutralEmFracJet1();
	if (chargedMuFracJet1_branch != 0) chargedMuFracJet1();
	if (electronHWW2011MVA_branch != 0) electronHWW2011MVA();
	if (egammaPOG2012MVA_branch != 0) egammaPOG2012MVA();
	if (muonHZZ2012IsoRingsMVA_branch != 0) muonHZZ2012IsoRingsMVA();
	if (vetoId_branch != 0) vetoId();
	if (looseId_branch != 0) looseId();
	if (mediumId_branch != 0) mediumId();
	if (tightId_branch != 0) tightId();
	if (pfmva_branch != 0) pfmva();
	if (sceta_branch != 0) sceta();
	if (scphi_branch != 0) scphi();
	if (scenergy_branch != 0) scenergy();
	if (chargesAgree_branch != 0) chargesAgree();
	if (eopin_branch != 0) eopin();
	if (ooemoop_branch != 0) ooemoop();
	if (fbrem_branch != 0) fbrem();
	if (detain_branch != 0) detain();
	if (dphiin_branch != 0) dphiin();
	if (hoe_branch != 0) hoe();
	if (hoetow_branch != 0) hoetow();
	if (sieie_branch != 0) sieie();
	if (d0vtx_branch != 0) d0vtx();
	if (dzvtx_branch != 0) dzvtx();
	if (vfitprob_branch != 0) vfitprob();
	if (mhit_branch != 0) mhit();
	if (ecaliso_branch != 0) ecaliso();
	if (hcaliso_branch != 0) hcaliso();
	if (trkiso_branch != 0) trkiso();
	if (pfemiso03_branch != 0) pfemiso03();
	if (pfchiso03_branch != 0) pfchiso03();
	if (pfnhiso03_branch != 0) pfnhiso03();
	if (pfemiso04_branch != 0) pfemiso04();
	if (pfchiso04_branch != 0) pfchiso04();
	if (pfnhiso04_branch != 0) pfnhiso04();
	if (radiso03_branch != 0) radiso03();
	if (radiso04_branch != 0) radiso04();
	if (iso2011_branch != 0) iso2011();
	if (ea04_branch != 0) ea04();
	if (ea03_branch != 0) ea03();
	if (dbeta03_branch != 0) dbeta03();
	if (dbeta04_branch != 0) dbeta04();
	if (el_test_pfchiso04_trkveto_branch != 0) el_test_pfchiso04_trkveto();
	if (el_test_pfchiso04_dzcut_branch != 0) el_test_pfchiso04_dzcut();
	if (el_test_pfchiso04_ebveto_branch != 0) el_test_pfchiso04_ebveto();
	if (el_test_pfemiso04_ebveto_branch != 0) el_test_pfemiso04_ebveto();
	if (eaem04_branch != 0) eaem04();
	if (eanh04_branch != 0) eanh04();
	if (gen_drs1_branch != 0) gen_drs1();
	if (gen_drs3_branch != 0) gen_drs3();
	if (HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch != 0) HLT_Photon22_R9Id90_HE10_Iso40_EBOnly();
	if (HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version();
	if (HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch != 0) HLT_Photon36_R9Id90_HE10_Iso40_EBOnly();
	if (HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version();
	if (HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch != 0) HLT_Photon50_R9Id90_HE10_Iso40_EBOnly();
	if (HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version();
	if (HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch != 0) HLT_Photon75_R9Id90_HE10_Iso40_EBOnly();
	if (HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version();
	if (HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch != 0) HLT_Photon90_R9Id90_HE10_Iso40_EBOnly();
	if (HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version();
	if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch != 0) HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag();
	if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch != 0) HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe();
	if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch != 0) HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version();
	if (HLT_Mu17_Mu8_TrailingLeg_tag_branch != 0) HLT_Mu17_Mu8_TrailingLeg_tag();
	if (HLT_Mu17_Mu8_TrailingLeg_probe_branch != 0) HLT_Mu17_Mu8_TrailingLeg_probe();
	if (HLT_Mu17_Mu8_TrailingLeg_version_branch != 0) HLT_Mu17_Mu8_TrailingLeg_version();
	if (HLT_Mu17_Mu8_LeadingLeg_tag_branch != 0) HLT_Mu17_Mu8_LeadingLeg_tag();
	if (HLT_Mu17_Mu8_LeadingLeg_probe_branch != 0) HLT_Mu17_Mu8_LeadingLeg_probe();
	if (HLT_Mu17_Mu8_LeadingLeg_version_branch != 0) HLT_Mu17_Mu8_LeadingLeg_version();
	if (HLT_Mu17_Mu8_tag_branch != 0) HLT_Mu17_Mu8_tag();
	if (HLT_Mu17_Mu8_probe_branch != 0) HLT_Mu17_Mu8_probe();
	if (HLT_Mu17_Mu8_version_branch != 0) HLT_Mu17_Mu8_version();
	if (HLT_Mu17_TkMu8_TrailingLeg_tag_branch != 0) HLT_Mu17_TkMu8_TrailingLeg_tag();
	if (HLT_Mu17_TkMu8_TrailingLeg_probe_branch != 0) HLT_Mu17_TkMu8_TrailingLeg_probe();
	if (HLT_Mu17_TkMu8_TrailingLeg_version_branch != 0) HLT_Mu17_TkMu8_TrailingLeg_version();
	if (HLT_Mu17_TkMu8_LeadingLeg_tag_branch != 0) HLT_Mu17_TkMu8_LeadingLeg_tag();
	if (HLT_Mu17_TkMu8_LeadingLeg_probe_branch != 0) HLT_Mu17_TkMu8_LeadingLeg_probe();
	if (HLT_Mu17_TkMu8_LeadingLeg_version_branch != 0) HLT_Mu17_TkMu8_LeadingLeg_version();
	if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch != 0) HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag();
	if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch != 0) HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe();
	if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch != 0) HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version();
	if (HLT_Mu17_TkMu8_tag_branch != 0) HLT_Mu17_TkMu8_tag();
	if (HLT_Mu17_TkMu8_probe_branch != 0) HLT_Mu17_TkMu8_probe();
	if (HLT_Mu17_TkMu8_version_branch != 0) HLT_Mu17_TkMu8_version();
	if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch != 0) HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag();
	if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch != 0) HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe();
	if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch != 0) HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version();
	if (HLT_IsoMu24_eta2p1_tag_branch != 0) HLT_IsoMu24_eta2p1_tag();
	if (HLT_IsoMu24_eta2p1_probe_branch != 0) HLT_IsoMu24_eta2p1_probe();
	if (HLT_IsoMu24_eta2p1_version_branch != 0) HLT_IsoMu24_eta2p1_version();
	if (HLT_Mu8_tag_branch != 0) HLT_Mu8_tag();
	if (HLT_Mu8_probe_branch != 0) HLT_Mu8_probe();
	if (HLT_Mu8_version_branch != 0) HLT_Mu8_version();
	if (HLT_Mu17_tag_branch != 0) HLT_Mu17_tag();
	if (HLT_Mu17_probe_branch != 0) HLT_Mu17_probe();
	if (HLT_Mu17_version_branch != 0) HLT_Mu17_version();
	if (HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch != 0) HLT_Ele17_Ele8_L1sL1DoubleEG137_tag();
	if (HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch != 0) HLT_Ele17_Ele8_L1sL1DoubleEG137_probe();
	if (HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch != 0) HLT_Ele17_Ele8_L1sL1DoubleEG137_version();
	if (HLT_Ele17_Ele8_LeadingLeg_tag_branch != 0) HLT_Ele17_Ele8_LeadingLeg_tag();
	if (HLT_Ele17_Ele8_LeadingLeg_probe_branch != 0) HLT_Ele17_Ele8_LeadingLeg_probe();
	if (HLT_Ele17_Ele8_LeadingLeg_version_branch != 0) HLT_Ele17_Ele8_LeadingLeg_version();
	if (HLT_Ele17_Ele8_TrailingLeg_tag_branch != 0) HLT_Ele17_Ele8_TrailingLeg_tag();
	if (HLT_Ele17_Ele8_TrailingLeg_probe_branch != 0) HLT_Ele17_Ele8_TrailingLeg_probe();
	if (HLT_Ele17_Ele8_TrailingLeg_version_branch != 0) HLT_Ele17_Ele8_TrailingLeg_version();
	if (HLT_Ele17_Ele8_tag_branch != 0) HLT_Ele17_Ele8_tag();
	if (HLT_Ele17_Ele8_probe_branch != 0) HLT_Ele17_Ele8_probe();
	if (HLT_Ele17_Ele8_version_branch != 0) HLT_Ele17_Ele8_version();
	if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch != 0) HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag();
	if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch != 0) HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe();
	if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch != 0) HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version();
	if (HLT_Ele27_WP80_tag_branch != 0) HLT_Ele27_WP80_tag();
	if (HLT_Ele27_WP80_probe_branch != 0) HLT_Ele27_WP80_probe();
	if (HLT_Ele27_WP80_version_branch != 0) HLT_Ele27_WP80_version();
	if (HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch != 0) HLT_Ele17_Ele8_Mass50_LeadingLeg_tag();
	if (HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch != 0) HLT_Ele17_Ele8_Mass50_LeadingLeg_probe();
	if (HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch != 0) HLT_Ele17_Ele8_Mass50_LeadingLeg_version();
	if (HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch != 0) HLT_Ele17_Ele8_Mass50_TrailingLeg_tag();
	if (HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch != 0) HLT_Ele17_Ele8_Mass50_TrailingLeg_probe();
	if (HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch != 0) HLT_Ele17_Ele8_Mass50_TrailingLeg_version();
	if (HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch != 0) HLT_Ele20_SC4_Mass50_LeadingLeg_tag();
	if (HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch != 0) HLT_Ele20_SC4_Mass50_LeadingLeg_probe();
	if (HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch != 0) HLT_Ele20_SC4_Mass50_LeadingLeg_version();
	if (HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch != 0) HLT_Ele20_SC4_Mass50_TrailingLeg_tag();
	if (HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch != 0) HLT_Ele20_SC4_Mass50_TrailingLeg_probe();
	if (HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch != 0) HLT_Ele20_SC4_Mass50_TrailingLeg_version();
	if (HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch != 0) HLT_Ele32_SC17_Mass50_LeadingLeg_tag();
	if (HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch != 0) HLT_Ele32_SC17_Mass50_LeadingLeg_probe();
	if (HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch != 0) HLT_Ele32_SC17_Mass50_LeadingLeg_version();
	if (HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch != 0) HLT_Ele32_SC17_Mass50_TrailingLeg_tag();
	if (HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch != 0) HLT_Ele32_SC17_Mass50_TrailingLeg_probe();
	if (HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch != 0) HLT_Ele32_SC17_Mass50_TrailingLeg_version();
	if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch != 0) HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag();
	if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch != 0) HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe();
	if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch != 0) HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version();
	if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch != 0) HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag();
	if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch != 0) HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe();
	if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch != 0) HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version();
	if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch != 0) HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag();
	if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch != 0) HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe();
	if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch != 0) HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version();
	if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch != 0) HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag();
	if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch != 0) HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe();
	if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch != 0) HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version();
	if (HLT_Ele8_CaloIdT_TrkIdVL_tag_branch != 0) HLT_Ele8_CaloIdT_TrkIdVL_tag();
	if (HLT_Ele8_CaloIdT_TrkIdVL_probe_branch != 0) HLT_Ele8_CaloIdT_TrkIdVL_probe();
	if (HLT_Ele8_CaloIdT_TrkIdVL_version_branch != 0) HLT_Ele8_CaloIdT_TrkIdVL_version();
}

	unsigned int &event()
	{
		if (not event_isLoaded) {
			if (event_branch != 0) {
				event_branch->GetEntry(index);
			} else { 
				printf("branch event_branch does not exist!\n");
				exit(1);
			}
			event_isLoaded = true;
		}
		return event_;
	}
	unsigned int &run()
	{
		if (not run_isLoaded) {
			if (run_branch != 0) {
				run_branch->GetEntry(index);
			} else { 
				printf("branch run_branch does not exist!\n");
				exit(1);
			}
			run_isLoaded = true;
		}
		return run_;
	}
	unsigned int &lumi()
	{
		if (not lumi_isLoaded) {
			if (lumi_branch != 0) {
				lumi_branch->GetEntry(index);
			} else { 
				printf("branch lumi_branch does not exist!\n");
				exit(1);
			}
			lumi_isLoaded = true;
		}
		return lumi_;
	}
	float &rnd()
	{
		if (not rnd_isLoaded) {
			if (rnd_branch != 0) {
				rnd_branch->GetEntry(index);
			} else { 
				printf("branch rnd_branch does not exist!\n");
				exit(1);
			}
			rnd_isLoaded = true;
		}
		return rnd_;
	}
	unsigned int &nvtx()
	{
		if (not nvtx_isLoaded) {
			if (nvtx_branch != 0) {
				nvtx_branch->GetEntry(index);
			} else { 
				printf("branch nvtx_branch does not exist!\n");
				exit(1);
			}
			nvtx_isLoaded = true;
		}
		return nvtx_;
	}
	unsigned int &npu()
	{
		if (not npu_isLoaded) {
			if (npu_branch != 0) {
				npu_branch->GetEntry(index);
			} else { 
				printf("branch npu_branch does not exist!\n");
				exit(1);
			}
			npu_isLoaded = true;
		}
		return npu_;
	}
	unsigned int &npuPlusOne()
	{
		if (not npuPlusOne_isLoaded) {
			if (npuPlusOne_branch != 0) {
				npuPlusOne_branch->GetEntry(index);
			} else { 
				printf("branch npuPlusOne_branch does not exist!\n");
				exit(1);
			}
			npuPlusOne_isLoaded = true;
		}
		return npuPlusOne_;
	}
	unsigned int &npuMinusOne()
	{
		if (not npuMinusOne_isLoaded) {
			if (npuMinusOne_branch != 0) {
				npuMinusOne_branch->GetEntry(index);
			} else { 
				printf("branch npuMinusOne_branch does not exist!\n");
				exit(1);
			}
			npuMinusOne_isLoaded = true;
		}
		return npuMinusOne_;
	}
	unsigned int &leptonSelection()
	{
		if (not leptonSelection_isLoaded) {
			if (leptonSelection_branch != 0) {
				leptonSelection_branch->GetEntry(index);
			} else { 
				printf("branch leptonSelection_branch does not exist!\n");
				exit(1);
			}
			leptonSelection_isLoaded = true;
		}
		return leptonSelection_;
	}
	unsigned int &eventSelection()
	{
		if (not eventSelection_isLoaded) {
			if (eventSelection_branch != 0) {
				eventSelection_branch->GetEntry(index);
			} else { 
				printf("branch eventSelection_branch does not exist!\n");
				exit(1);
			}
			eventSelection_isLoaded = true;
		}
		return eventSelection_;
	}
	float &rhoIsoAll()
	{
		if (not rhoIsoAll_isLoaded) {
			if (rhoIsoAll_branch != 0) {
				rhoIsoAll_branch->GetEntry(index);
			} else { 
				printf("branch rhoIsoAll_branch does not exist!\n");
				exit(1);
			}
			rhoIsoAll_isLoaded = true;
		}
		return rhoIsoAll_;
	}
	float &rhoIsoAllCentral()
	{
		if (not rhoIsoAllCentral_isLoaded) {
			if (rhoIsoAllCentral_branch != 0) {
				rhoIsoAllCentral_branch->GetEntry(index);
			} else { 
				printf("branch rhoIsoAllCentral_branch does not exist!\n");
				exit(1);
			}
			rhoIsoAllCentral_isLoaded = true;
		}
		return rhoIsoAllCentral_;
	}
	float &rhoIsoNeutral()
	{
		if (not rhoIsoNeutral_isLoaded) {
			if (rhoIsoNeutral_branch != 0) {
				rhoIsoNeutral_branch->GetEntry(index);
			} else { 
				printf("branch rhoIsoNeutral_branch does not exist!\n");
				exit(1);
			}
			rhoIsoNeutral_isLoaded = true;
		}
		return rhoIsoNeutral_;
	}
	float &tagAndProbeMass()
	{
		if (not tagAndProbeMass_isLoaded) {
			if (tagAndProbeMass_branch != 0) {
				tagAndProbeMass_branch->GetEntry(index);
			} else { 
				printf("branch tagAndProbeMass_branch does not exist!\n");
				exit(1);
			}
			tagAndProbeMass_isLoaded = true;
		}
		return tagAndProbeMass_;
	}
	bool &	tagAndProbeIsRandom()
	{
		if (not tagAndProbeIsRandom_isLoaded) {
			if (tagAndProbeIsRandom_branch != 0) {
				tagAndProbeIsRandom_branch->GetEntry(index);
			} else { 
				printf("branch tagAndProbeIsRandom_branch does not exist!\n");
				exit(1);
			}
			tagAndProbeIsRandom_isLoaded = true;
		}
		return tagAndProbeIsRandom_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &tag()
	{
		if (not tag_isLoaded) {
			if (tag_branch != 0) {
				tag_branch->GetEntry(index);
			} else { 
				printf("branch tag_branch does not exist!\n");
				exit(1);
			}
			tag_isLoaded = true;
		}
		return *tag_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &probe()
	{
		if (not probe_isLoaded) {
			if (probe_branch != 0) {
				probe_branch->GetEntry(index);
			} else { 
				printf("branch probe_branch does not exist!\n");
				exit(1);
			}
			probe_isLoaded = true;
		}
		return *probe_;
	}
	int &qTag()
	{
		if (not qTag_isLoaded) {
			if (qTag_branch != 0) {
				qTag_branch->GetEntry(index);
			} else { 
				printf("branch qTag_branch does not exist!\n");
				exit(1);
			}
			qTag_isLoaded = true;
		}
		return qTag_;
	}
	int &qProbe()
	{
		if (not qProbe_isLoaded) {
			if (qProbe_branch != 0) {
				qProbe_branch->GetEntry(index);
			} else { 
				printf("branch qProbe_branch does not exist!\n");
				exit(1);
			}
			qProbe_isLoaded = true;
		}
		return qProbe_;
	}
	float &scale1fb()
	{
		if (not scale1fb_isLoaded) {
			if (scale1fb_branch != 0) {
				scale1fb_branch->GetEntry(index);
			} else { 
				printf("branch scale1fb_branch does not exist!\n");
				exit(1);
			}
			scale1fb_isLoaded = true;
		}
		return scale1fb_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet1()
	{
		if (not jet1_isLoaded) {
			if (jet1_branch != 0) {
				jet1_branch->GetEntry(index);
			} else { 
				printf("branch jet1_branch does not exist!\n");
				exit(1);
			}
			jet1_isLoaded = true;
		}
		return *jet1_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet2()
	{
		if (not jet2_isLoaded) {
			if (jet2_branch != 0) {
				jet2_branch->GetEntry(index);
			} else { 
				printf("branch jet2_branch does not exist!\n");
				exit(1);
			}
			jet2_isLoaded = true;
		}
		return *jet2_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet3()
	{
		if (not jet3_isLoaded) {
			if (jet3_branch != 0) {
				jet3_branch->GetEntry(index);
			} else { 
				printf("branch jet3_branch does not exist!\n");
				exit(1);
			}
			jet3_isLoaded = true;
		}
		return *jet3_;
	}
	float &met()
	{
		if (not met_isLoaded) {
			if (met_branch != 0) {
				met_branch->GetEntry(index);
			} else { 
				printf("branch met_branch does not exist!\n");
				exit(1);
			}
			met_isLoaded = true;
		}
		return met_;
	}
	float &metPhi()
	{
		if (not metPhi_isLoaded) {
			if (metPhi_branch != 0) {
				metPhi_branch->GetEntry(index);
			} else { 
				printf("branch metPhi_branch does not exist!\n");
				exit(1);
			}
			metPhi_isLoaded = true;
		}
		return metPhi_;
	}
	float &trackMet()
	{
		if (not trackMet_isLoaded) {
			if (trackMet_branch != 0) {
				trackMet_branch->GetEntry(index);
			} else { 
				printf("branch trackMet_branch does not exist!\n");
				exit(1);
			}
			trackMet_isLoaded = true;
		}
		return trackMet_;
	}
	float &trackMetPhi()
	{
		if (not trackMetPhi_isLoaded) {
			if (trackMetPhi_branch != 0) {
				trackMetPhi_branch->GetEntry(index);
			} else { 
				printf("branch trackMetPhi_branch does not exist!\n");
				exit(1);
			}
			trackMetPhi_isLoaded = true;
		}
		return trackMetPhi_;
	}
	unsigned int &njets()
	{
		if (not njets_isLoaded) {
			if (njets_branch != 0) {
				njets_branch->GetEntry(index);
			} else { 
				printf("branch njets_branch does not exist!\n");
				exit(1);
			}
			njets_isLoaded = true;
		}
		return njets_;
	}
	unsigned int &nleps()
	{
		if (not nleps_isLoaded) {
			if (nleps_branch != 0) {
				nleps_branch->GetEntry(index);
			} else { 
				printf("branch nleps_branch does not exist!\n");
				exit(1);
			}
			nleps_isLoaded = true;
		}
		return nleps_;
	}
	float &hltPrescale()
	{
		if (not hltPrescale_isLoaded) {
			if (hltPrescale_branch != 0) {
				hltPrescale_branch->GetEntry(index);
			} else { 
				printf("branch hltPrescale_branch does not exist!\n");
				exit(1);
			}
			hltPrescale_isLoaded = true;
		}
		return hltPrescale_;
	}
	float &sumet()
	{
		if (not sumet_isLoaded) {
			if (sumet_branch != 0) {
				sumet_branch->GetEntry(index);
			} else { 
				printf("branch sumet_branch does not exist!\n");
				exit(1);
			}
			sumet_isLoaded = true;
		}
		return sumet_;
	}
	float &metSig()
	{
		if (not metSig_isLoaded) {
			if (metSig_branch != 0) {
				metSig_branch->GetEntry(index);
			} else { 
				printf("branch metSig_branch does not exist!\n");
				exit(1);
			}
			metSig_isLoaded = true;
		}
		return metSig_;
	}
	float &mt()
	{
		if (not mt_isLoaded) {
			if (mt_branch != 0) {
				mt_branch->GetEntry(index);
			} else { 
				printf("branch mt_branch does not exist!\n");
				exit(1);
			}
			mt_isLoaded = true;
		}
		return mt_;
	}
	float &dPhiProbeJet1()
	{
		if (not dPhiProbeJet1_isLoaded) {
			if (dPhiProbeJet1_branch != 0) {
				dPhiProbeJet1_branch->GetEntry(index);
			} else { 
				printf("branch dPhiProbeJet1_branch does not exist!\n");
				exit(1);
			}
			dPhiProbeJet1_isLoaded = true;
		}
		return dPhiProbeJet1_;
	}
	float &chargedEmFracJet1()
	{
		if (not chargedEmFracJet1_isLoaded) {
			if (chargedEmFracJet1_branch != 0) {
				chargedEmFracJet1_branch->GetEntry(index);
			} else { 
				printf("branch chargedEmFracJet1_branch does not exist!\n");
				exit(1);
			}
			chargedEmFracJet1_isLoaded = true;
		}
		return chargedEmFracJet1_;
	}
	float &neutralEmFracJet1()
	{
		if (not neutralEmFracJet1_isLoaded) {
			if (neutralEmFracJet1_branch != 0) {
				neutralEmFracJet1_branch->GetEntry(index);
			} else { 
				printf("branch neutralEmFracJet1_branch does not exist!\n");
				exit(1);
			}
			neutralEmFracJet1_isLoaded = true;
		}
		return neutralEmFracJet1_;
	}
	float &chargedMuFracJet1()
	{
		if (not chargedMuFracJet1_isLoaded) {
			if (chargedMuFracJet1_branch != 0) {
				chargedMuFracJet1_branch->GetEntry(index);
			} else { 
				printf("branch chargedMuFracJet1_branch does not exist!\n");
				exit(1);
			}
			chargedMuFracJet1_isLoaded = true;
		}
		return chargedMuFracJet1_;
	}
	float &electronHWW2011MVA()
	{
		if (not electronHWW2011MVA_isLoaded) {
			if (electronHWW2011MVA_branch != 0) {
				electronHWW2011MVA_branch->GetEntry(index);
			} else { 
				printf("branch electronHWW2011MVA_branch does not exist!\n");
				exit(1);
			}
			electronHWW2011MVA_isLoaded = true;
		}
		return electronHWW2011MVA_;
	}
	float &egammaPOG2012MVA()
	{
		if (not egammaPOG2012MVA_isLoaded) {
			if (egammaPOG2012MVA_branch != 0) {
				egammaPOG2012MVA_branch->GetEntry(index);
			} else { 
				printf("branch egammaPOG2012MVA_branch does not exist!\n");
				exit(1);
			}
			egammaPOG2012MVA_isLoaded = true;
		}
		return egammaPOG2012MVA_;
	}
	float &muonHZZ2012IsoRingsMVA()
	{
		if (not muonHZZ2012IsoRingsMVA_isLoaded) {
			if (muonHZZ2012IsoRingsMVA_branch != 0) {
				muonHZZ2012IsoRingsMVA_branch->GetEntry(index);
			} else { 
				printf("branch muonHZZ2012IsoRingsMVA_branch does not exist!\n");
				exit(1);
			}
			muonHZZ2012IsoRingsMVA_isLoaded = true;
		}
		return muonHZZ2012IsoRingsMVA_;
	}
	unsigned int &vetoId()
	{
		if (not vetoId_isLoaded) {
			if (vetoId_branch != 0) {
				vetoId_branch->GetEntry(index);
			} else { 
				printf("branch vetoId_branch does not exist!\n");
				exit(1);
			}
			vetoId_isLoaded = true;
		}
		return vetoId_;
	}
	unsigned int &looseId()
	{
		if (not looseId_isLoaded) {
			if (looseId_branch != 0) {
				looseId_branch->GetEntry(index);
			} else { 
				printf("branch looseId_branch does not exist!\n");
				exit(1);
			}
			looseId_isLoaded = true;
		}
		return looseId_;
	}
	unsigned int &mediumId()
	{
		if (not mediumId_isLoaded) {
			if (mediumId_branch != 0) {
				mediumId_branch->GetEntry(index);
			} else { 
				printf("branch mediumId_branch does not exist!\n");
				exit(1);
			}
			mediumId_isLoaded = true;
		}
		return mediumId_;
	}
	unsigned int &tightId()
	{
		if (not tightId_isLoaded) {
			if (tightId_branch != 0) {
				tightId_branch->GetEntry(index);
			} else { 
				printf("branch tightId_branch does not exist!\n");
				exit(1);
			}
			tightId_isLoaded = true;
		}
		return tightId_;
	}
	float &pfmva()
	{
		if (not pfmva_isLoaded) {
			if (pfmva_branch != 0) {
				pfmva_branch->GetEntry(index);
			} else { 
				printf("branch pfmva_branch does not exist!\n");
				exit(1);
			}
			pfmva_isLoaded = true;
		}
		return pfmva_;
	}
	float &sceta()
	{
		if (not sceta_isLoaded) {
			if (sceta_branch != 0) {
				sceta_branch->GetEntry(index);
			} else { 
				printf("branch sceta_branch does not exist!\n");
				exit(1);
			}
			sceta_isLoaded = true;
		}
		return sceta_;
	}
	float &scphi()
	{
		if (not scphi_isLoaded) {
			if (scphi_branch != 0) {
				scphi_branch->GetEntry(index);
			} else { 
				printf("branch scphi_branch does not exist!\n");
				exit(1);
			}
			scphi_isLoaded = true;
		}
		return scphi_;
	}
	float &scenergy()
	{
		if (not scenergy_isLoaded) {
			if (scenergy_branch != 0) {
				scenergy_branch->GetEntry(index);
			} else { 
				printf("branch scenergy_branch does not exist!\n");
				exit(1);
			}
			scenergy_isLoaded = true;
		}
		return scenergy_;
	}
	bool &	chargesAgree()
	{
		if (not chargesAgree_isLoaded) {
			if (chargesAgree_branch != 0) {
				chargesAgree_branch->GetEntry(index);
			} else { 
				printf("branch chargesAgree_branch does not exist!\n");
				exit(1);
			}
			chargesAgree_isLoaded = true;
		}
		return chargesAgree_;
	}
	float &eopin()
	{
		if (not eopin_isLoaded) {
			if (eopin_branch != 0) {
				eopin_branch->GetEntry(index);
			} else { 
				printf("branch eopin_branch does not exist!\n");
				exit(1);
			}
			eopin_isLoaded = true;
		}
		return eopin_;
	}
	float &ooemoop()
	{
		if (not ooemoop_isLoaded) {
			if (ooemoop_branch != 0) {
				ooemoop_branch->GetEntry(index);
			} else { 
				printf("branch ooemoop_branch does not exist!\n");
				exit(1);
			}
			ooemoop_isLoaded = true;
		}
		return ooemoop_;
	}
	float &fbrem()
	{
		if (not fbrem_isLoaded) {
			if (fbrem_branch != 0) {
				fbrem_branch->GetEntry(index);
			} else { 
				printf("branch fbrem_branch does not exist!\n");
				exit(1);
			}
			fbrem_isLoaded = true;
		}
		return fbrem_;
	}
	float &detain()
	{
		if (not detain_isLoaded) {
			if (detain_branch != 0) {
				detain_branch->GetEntry(index);
			} else { 
				printf("branch detain_branch does not exist!\n");
				exit(1);
			}
			detain_isLoaded = true;
		}
		return detain_;
	}
	float &dphiin()
	{
		if (not dphiin_isLoaded) {
			if (dphiin_branch != 0) {
				dphiin_branch->GetEntry(index);
			} else { 
				printf("branch dphiin_branch does not exist!\n");
				exit(1);
			}
			dphiin_isLoaded = true;
		}
		return dphiin_;
	}
	float &hoe()
	{
		if (not hoe_isLoaded) {
			if (hoe_branch != 0) {
				hoe_branch->GetEntry(index);
			} else { 
				printf("branch hoe_branch does not exist!\n");
				exit(1);
			}
			hoe_isLoaded = true;
		}
		return hoe_;
	}
	float &hoetow()
	{
		if (not hoetow_isLoaded) {
			if (hoetow_branch != 0) {
				hoetow_branch->GetEntry(index);
			} else { 
				printf("branch hoetow_branch does not exist!\n");
				exit(1);
			}
			hoetow_isLoaded = true;
		}
		return hoetow_;
	}
	float &sieie()
	{
		if (not sieie_isLoaded) {
			if (sieie_branch != 0) {
				sieie_branch->GetEntry(index);
			} else { 
				printf("branch sieie_branch does not exist!\n");
				exit(1);
			}
			sieie_isLoaded = true;
		}
		return sieie_;
	}
	float &d0vtx()
	{
		if (not d0vtx_isLoaded) {
			if (d0vtx_branch != 0) {
				d0vtx_branch->GetEntry(index);
			} else { 
				printf("branch d0vtx_branch does not exist!\n");
				exit(1);
			}
			d0vtx_isLoaded = true;
		}
		return d0vtx_;
	}
	float &dzvtx()
	{
		if (not dzvtx_isLoaded) {
			if (dzvtx_branch != 0) {
				dzvtx_branch->GetEntry(index);
			} else { 
				printf("branch dzvtx_branch does not exist!\n");
				exit(1);
			}
			dzvtx_isLoaded = true;
		}
		return dzvtx_;
	}
	bool &	vfitprob()
	{
		if (not vfitprob_isLoaded) {
			if (vfitprob_branch != 0) {
				vfitprob_branch->GetEntry(index);
			} else { 
				printf("branch vfitprob_branch does not exist!\n");
				exit(1);
			}
			vfitprob_isLoaded = true;
		}
		return vfitprob_;
	}
	float &mhit()
	{
		if (not mhit_isLoaded) {
			if (mhit_branch != 0) {
				mhit_branch->GetEntry(index);
			} else { 
				printf("branch mhit_branch does not exist!\n");
				exit(1);
			}
			mhit_isLoaded = true;
		}
		return mhit_;
	}
	float &ecaliso()
	{
		if (not ecaliso_isLoaded) {
			if (ecaliso_branch != 0) {
				ecaliso_branch->GetEntry(index);
			} else { 
				printf("branch ecaliso_branch does not exist!\n");
				exit(1);
			}
			ecaliso_isLoaded = true;
		}
		return ecaliso_;
	}
	float &hcaliso()
	{
		if (not hcaliso_isLoaded) {
			if (hcaliso_branch != 0) {
				hcaliso_branch->GetEntry(index);
			} else { 
				printf("branch hcaliso_branch does not exist!\n");
				exit(1);
			}
			hcaliso_isLoaded = true;
		}
		return hcaliso_;
	}
	float &trkiso()
	{
		if (not trkiso_isLoaded) {
			if (trkiso_branch != 0) {
				trkiso_branch->GetEntry(index);
			} else { 
				printf("branch trkiso_branch does not exist!\n");
				exit(1);
			}
			trkiso_isLoaded = true;
		}
		return trkiso_;
	}
	float &pfemiso03()
	{
		if (not pfemiso03_isLoaded) {
			if (pfemiso03_branch != 0) {
				pfemiso03_branch->GetEntry(index);
			} else { 
				printf("branch pfemiso03_branch does not exist!\n");
				exit(1);
			}
			pfemiso03_isLoaded = true;
		}
		return pfemiso03_;
	}
	float &pfchiso03()
	{
		if (not pfchiso03_isLoaded) {
			if (pfchiso03_branch != 0) {
				pfchiso03_branch->GetEntry(index);
			} else { 
				printf("branch pfchiso03_branch does not exist!\n");
				exit(1);
			}
			pfchiso03_isLoaded = true;
		}
		return pfchiso03_;
	}
	float &pfnhiso03()
	{
		if (not pfnhiso03_isLoaded) {
			if (pfnhiso03_branch != 0) {
				pfnhiso03_branch->GetEntry(index);
			} else { 
				printf("branch pfnhiso03_branch does not exist!\n");
				exit(1);
			}
			pfnhiso03_isLoaded = true;
		}
		return pfnhiso03_;
	}
	float &pfemiso04()
	{
		if (not pfemiso04_isLoaded) {
			if (pfemiso04_branch != 0) {
				pfemiso04_branch->GetEntry(index);
			} else { 
				printf("branch pfemiso04_branch does not exist!\n");
				exit(1);
			}
			pfemiso04_isLoaded = true;
		}
		return pfemiso04_;
	}
	float &pfchiso04()
	{
		if (not pfchiso04_isLoaded) {
			if (pfchiso04_branch != 0) {
				pfchiso04_branch->GetEntry(index);
			} else { 
				printf("branch pfchiso04_branch does not exist!\n");
				exit(1);
			}
			pfchiso04_isLoaded = true;
		}
		return pfchiso04_;
	}
	float &pfnhiso04()
	{
		if (not pfnhiso04_isLoaded) {
			if (pfnhiso04_branch != 0) {
				pfnhiso04_branch->GetEntry(index);
			} else { 
				printf("branch pfnhiso04_branch does not exist!\n");
				exit(1);
			}
			pfnhiso04_isLoaded = true;
		}
		return pfnhiso04_;
	}
	float &radiso03()
	{
		if (not radiso03_isLoaded) {
			if (radiso03_branch != 0) {
				radiso03_branch->GetEntry(index);
			} else { 
				printf("branch radiso03_branch does not exist!\n");
				exit(1);
			}
			radiso03_isLoaded = true;
		}
		return radiso03_;
	}
	float &radiso04()
	{
		if (not radiso04_isLoaded) {
			if (radiso04_branch != 0) {
				radiso04_branch->GetEntry(index);
			} else { 
				printf("branch radiso04_branch does not exist!\n");
				exit(1);
			}
			radiso04_isLoaded = true;
		}
		return radiso04_;
	}
	float &iso2011()
	{
		if (not iso2011_isLoaded) {
			if (iso2011_branch != 0) {
				iso2011_branch->GetEntry(index);
			} else { 
				printf("branch iso2011_branch does not exist!\n");
				exit(1);
			}
			iso2011_isLoaded = true;
		}
		return iso2011_;
	}
	float &ea04()
	{
		if (not ea04_isLoaded) {
			if (ea04_branch != 0) {
				ea04_branch->GetEntry(index);
			} else { 
				printf("branch ea04_branch does not exist!\n");
				exit(1);
			}
			ea04_isLoaded = true;
		}
		return ea04_;
	}
	float &ea03()
	{
		if (not ea03_isLoaded) {
			if (ea03_branch != 0) {
				ea03_branch->GetEntry(index);
			} else { 
				printf("branch ea03_branch does not exist!\n");
				exit(1);
			}
			ea03_isLoaded = true;
		}
		return ea03_;
	}
	float &dbeta03()
	{
		if (not dbeta03_isLoaded) {
			if (dbeta03_branch != 0) {
				dbeta03_branch->GetEntry(index);
			} else { 
				printf("branch dbeta03_branch does not exist!\n");
				exit(1);
			}
			dbeta03_isLoaded = true;
		}
		return dbeta03_;
	}
	float &dbeta04()
	{
		if (not dbeta04_isLoaded) {
			if (dbeta04_branch != 0) {
				dbeta04_branch->GetEntry(index);
			} else { 
				printf("branch dbeta04_branch does not exist!\n");
				exit(1);
			}
			dbeta04_isLoaded = true;
		}
		return dbeta04_;
	}
	float &el_test_pfchiso04_trkveto()
	{
		if (not el_test_pfchiso04_trkveto_isLoaded) {
			if (el_test_pfchiso04_trkveto_branch != 0) {
				el_test_pfchiso04_trkveto_branch->GetEntry(index);
			} else { 
				printf("branch el_test_pfchiso04_trkveto_branch does not exist!\n");
				exit(1);
			}
			el_test_pfchiso04_trkveto_isLoaded = true;
		}
		return el_test_pfchiso04_trkveto_;
	}
	float &el_test_pfchiso04_dzcut()
	{
		if (not el_test_pfchiso04_dzcut_isLoaded) {
			if (el_test_pfchiso04_dzcut_branch != 0) {
				el_test_pfchiso04_dzcut_branch->GetEntry(index);
			} else { 
				printf("branch el_test_pfchiso04_dzcut_branch does not exist!\n");
				exit(1);
			}
			el_test_pfchiso04_dzcut_isLoaded = true;
		}
		return el_test_pfchiso04_dzcut_;
	}
	float &el_test_pfchiso04_ebveto()
	{
		if (not el_test_pfchiso04_ebveto_isLoaded) {
			if (el_test_pfchiso04_ebveto_branch != 0) {
				el_test_pfchiso04_ebveto_branch->GetEntry(index);
			} else { 
				printf("branch el_test_pfchiso04_ebveto_branch does not exist!\n");
				exit(1);
			}
			el_test_pfchiso04_ebveto_isLoaded = true;
		}
		return el_test_pfchiso04_ebveto_;
	}
	float &el_test_pfemiso04_ebveto()
	{
		if (not el_test_pfemiso04_ebveto_isLoaded) {
			if (el_test_pfemiso04_ebveto_branch != 0) {
				el_test_pfemiso04_ebveto_branch->GetEntry(index);
			} else { 
				printf("branch el_test_pfemiso04_ebveto_branch does not exist!\n");
				exit(1);
			}
			el_test_pfemiso04_ebveto_isLoaded = true;
		}
		return el_test_pfemiso04_ebveto_;
	}
	float &eaem04()
	{
		if (not eaem04_isLoaded) {
			if (eaem04_branch != 0) {
				eaem04_branch->GetEntry(index);
			} else { 
				printf("branch eaem04_branch does not exist!\n");
				exit(1);
			}
			eaem04_isLoaded = true;
		}
		return eaem04_;
	}
	float &eanh04()
	{
		if (not eanh04_isLoaded) {
			if (eanh04_branch != 0) {
				eanh04_branch->GetEntry(index);
			} else { 
				printf("branch eanh04_branch does not exist!\n");
				exit(1);
			}
			eanh04_isLoaded = true;
		}
		return eanh04_;
	}
	float &gen_drs1()
	{
		if (not gen_drs1_isLoaded) {
			if (gen_drs1_branch != 0) {
				gen_drs1_branch->GetEntry(index);
			} else { 
				printf("branch gen_drs1_branch does not exist!\n");
				exit(1);
			}
			gen_drs1_isLoaded = true;
		}
		return gen_drs1_;
	}
	float &gen_drs3()
	{
		if (not gen_drs3_isLoaded) {
			if (gen_drs3_branch != 0) {
				gen_drs3_branch->GetEntry(index);
			} else { 
				printf("branch gen_drs3_branch does not exist!\n");
				exit(1);
			}
			gen_drs3_isLoaded = true;
		}
		return gen_drs3_;
	}
	unsigned int &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly()
	{
		if (not HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_isLoaded) {
			if (HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch != 0) {
				HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_isLoaded = true;
		}
		return HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_;
	}
	unsigned int &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version()
	{
		if (not HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_isLoaded) {
			if (HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) {
				HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = true;
		}
		return HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version_;
	}
	unsigned int &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly()
	{
		if (not HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_isLoaded) {
			if (HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch != 0) {
				HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_isLoaded = true;
		}
		return HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_;
	}
	unsigned int &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version()
	{
		if (not HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_isLoaded) {
			if (HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) {
				HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = true;
		}
		return HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version_;
	}
	unsigned int &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly()
	{
		if (not HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_isLoaded) {
			if (HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch != 0) {
				HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_isLoaded = true;
		}
		return HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_;
	}
	unsigned int &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version()
	{
		if (not HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_isLoaded) {
			if (HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) {
				HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = true;
		}
		return HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version_;
	}
	unsigned int &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly()
	{
		if (not HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_isLoaded) {
			if (HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch != 0) {
				HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_isLoaded = true;
		}
		return HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_;
	}
	unsigned int &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version()
	{
		if (not HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_isLoaded) {
			if (HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) {
				HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = true;
		}
		return HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version_;
	}
	unsigned int &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly()
	{
		if (not HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_isLoaded) {
			if (HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch != 0) {
				HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_isLoaded = true;
		}
		return HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_;
	}
	unsigned int &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version()
	{
		if (not HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_isLoaded) {
			if (HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch != 0) {
				HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_isLoaded = true;
		}
		return HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version_;
	}
	unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag()
	{
		if (not HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_isLoaded) {
			if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch != 0) {
				HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_isLoaded = true;
		}
		return HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_;
	}
	unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe()
	{
		if (not HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_isLoaded) {
			if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch != 0) {
				HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_isLoaded = true;
		}
		return HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_;
	}
	unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version()
	{
		if (not HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded) {
			if (HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch != 0) {
				HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded = true;
		}
		return HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_;
	}
	unsigned int &HLT_Mu17_Mu8_TrailingLeg_tag()
	{
		if (not HLT_Mu17_Mu8_TrailingLeg_tag_isLoaded) {
			if (HLT_Mu17_Mu8_TrailingLeg_tag_branch != 0) {
				HLT_Mu17_Mu8_TrailingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_TrailingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_TrailingLeg_tag_isLoaded = true;
		}
		return HLT_Mu17_Mu8_TrailingLeg_tag_;
	}
	unsigned int &HLT_Mu17_Mu8_TrailingLeg_probe()
	{
		if (not HLT_Mu17_Mu8_TrailingLeg_probe_isLoaded) {
			if (HLT_Mu17_Mu8_TrailingLeg_probe_branch != 0) {
				HLT_Mu17_Mu8_TrailingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_TrailingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_TrailingLeg_probe_isLoaded = true;
		}
		return HLT_Mu17_Mu8_TrailingLeg_probe_;
	}
	unsigned int &HLT_Mu17_Mu8_TrailingLeg_version()
	{
		if (not HLT_Mu17_Mu8_TrailingLeg_version_isLoaded) {
			if (HLT_Mu17_Mu8_TrailingLeg_version_branch != 0) {
				HLT_Mu17_Mu8_TrailingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_TrailingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_TrailingLeg_version_isLoaded = true;
		}
		return HLT_Mu17_Mu8_TrailingLeg_version_;
	}
	unsigned int &HLT_Mu17_Mu8_LeadingLeg_tag()
	{
		if (not HLT_Mu17_Mu8_LeadingLeg_tag_isLoaded) {
			if (HLT_Mu17_Mu8_LeadingLeg_tag_branch != 0) {
				HLT_Mu17_Mu8_LeadingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_LeadingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_LeadingLeg_tag_isLoaded = true;
		}
		return HLT_Mu17_Mu8_LeadingLeg_tag_;
	}
	unsigned int &HLT_Mu17_Mu8_LeadingLeg_probe()
	{
		if (not HLT_Mu17_Mu8_LeadingLeg_probe_isLoaded) {
			if (HLT_Mu17_Mu8_LeadingLeg_probe_branch != 0) {
				HLT_Mu17_Mu8_LeadingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_LeadingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_LeadingLeg_probe_isLoaded = true;
		}
		return HLT_Mu17_Mu8_LeadingLeg_probe_;
	}
	unsigned int &HLT_Mu17_Mu8_LeadingLeg_version()
	{
		if (not HLT_Mu17_Mu8_LeadingLeg_version_isLoaded) {
			if (HLT_Mu17_Mu8_LeadingLeg_version_branch != 0) {
				HLT_Mu17_Mu8_LeadingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_LeadingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_LeadingLeg_version_isLoaded = true;
		}
		return HLT_Mu17_Mu8_LeadingLeg_version_;
	}
	unsigned int &HLT_Mu17_Mu8_tag()
	{
		if (not HLT_Mu17_Mu8_tag_isLoaded) {
			if (HLT_Mu17_Mu8_tag_branch != 0) {
				HLT_Mu17_Mu8_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_tag_isLoaded = true;
		}
		return HLT_Mu17_Mu8_tag_;
	}
	unsigned int &HLT_Mu17_Mu8_probe()
	{
		if (not HLT_Mu17_Mu8_probe_isLoaded) {
			if (HLT_Mu17_Mu8_probe_branch != 0) {
				HLT_Mu17_Mu8_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_probe_isLoaded = true;
		}
		return HLT_Mu17_Mu8_probe_;
	}
	unsigned int &HLT_Mu17_Mu8_version()
	{
		if (not HLT_Mu17_Mu8_version_isLoaded) {
			if (HLT_Mu17_Mu8_version_branch != 0) {
				HLT_Mu17_Mu8_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_Mu8_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_Mu8_version_isLoaded = true;
		}
		return HLT_Mu17_Mu8_version_;
	}
	unsigned int &HLT_Mu17_TkMu8_TrailingLeg_tag()
	{
		if (not HLT_Mu17_TkMu8_TrailingLeg_tag_isLoaded) {
			if (HLT_Mu17_TkMu8_TrailingLeg_tag_branch != 0) {
				HLT_Mu17_TkMu8_TrailingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_TrailingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_TrailingLeg_tag_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_TrailingLeg_tag_;
	}
	unsigned int &HLT_Mu17_TkMu8_TrailingLeg_probe()
	{
		if (not HLT_Mu17_TkMu8_TrailingLeg_probe_isLoaded) {
			if (HLT_Mu17_TkMu8_TrailingLeg_probe_branch != 0) {
				HLT_Mu17_TkMu8_TrailingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_TrailingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_TrailingLeg_probe_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_TrailingLeg_probe_;
	}
	unsigned int &HLT_Mu17_TkMu8_TrailingLeg_version()
	{
		if (not HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded) {
			if (HLT_Mu17_TkMu8_TrailingLeg_version_branch != 0) {
				HLT_Mu17_TkMu8_TrailingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_TrailingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_TrailingLeg_version_;
	}
	unsigned int &HLT_Mu17_TkMu8_LeadingLeg_tag()
	{
		if (not HLT_Mu17_TkMu8_LeadingLeg_tag_isLoaded) {
			if (HLT_Mu17_TkMu8_LeadingLeg_tag_branch != 0) {
				HLT_Mu17_TkMu8_LeadingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_LeadingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_LeadingLeg_tag_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_LeadingLeg_tag_;
	}
	unsigned int &HLT_Mu17_TkMu8_LeadingLeg_probe()
	{
		if (not HLT_Mu17_TkMu8_LeadingLeg_probe_isLoaded) {
			if (HLT_Mu17_TkMu8_LeadingLeg_probe_branch != 0) {
				HLT_Mu17_TkMu8_LeadingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_LeadingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_LeadingLeg_probe_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_LeadingLeg_probe_;
	}
	unsigned int &HLT_Mu17_TkMu8_LeadingLeg_version()
	{
		if (not HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded) {
			if (HLT_Mu17_TkMu8_LeadingLeg_version_branch != 0) {
				HLT_Mu17_TkMu8_LeadingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_LeadingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_LeadingLeg_version_;
	}
	unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag()
	{
		if (not HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_isLoaded) {
			if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch != 0) {
				HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_;
	}
	unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe()
	{
		if (not HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_isLoaded) {
			if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch != 0) {
				HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_;
	}
	unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version()
	{
		if (not HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded) {
			if (HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch != 0) {
				HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_;
	}
	unsigned int &HLT_Mu17_TkMu8_tag()
	{
		if (not HLT_Mu17_TkMu8_tag_isLoaded) {
			if (HLT_Mu17_TkMu8_tag_branch != 0) {
				HLT_Mu17_TkMu8_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_tag_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_tag_;
	}
	unsigned int &HLT_Mu17_TkMu8_probe()
	{
		if (not HLT_Mu17_TkMu8_probe_isLoaded) {
			if (HLT_Mu17_TkMu8_probe_branch != 0) {
				HLT_Mu17_TkMu8_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_probe_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_probe_;
	}
	unsigned int &HLT_Mu17_TkMu8_version()
	{
		if (not HLT_Mu17_TkMu8_version_isLoaded) {
			if (HLT_Mu17_TkMu8_version_branch != 0) {
				HLT_Mu17_TkMu8_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_TkMu8_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_TkMu8_version_isLoaded = true;
		}
		return HLT_Mu17_TkMu8_version_;
	}
	unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag()
	{
		if (not HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_isLoaded) {
			if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch != 0) {
				HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_isLoaded = true;
		}
		return HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_;
	}
	unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe()
	{
		if (not HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_isLoaded) {
			if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch != 0) {
				HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_isLoaded = true;
		}
		return HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_;
	}
	unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version()
	{
		if (not HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded) {
			if (HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch != 0) {
				HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch does not exist!\n");
				exit(1);
			}
			HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded = true;
		}
		return HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_;
	}
	unsigned int &HLT_IsoMu24_eta2p1_tag()
	{
		if (not HLT_IsoMu24_eta2p1_tag_isLoaded) {
			if (HLT_IsoMu24_eta2p1_tag_branch != 0) {
				HLT_IsoMu24_eta2p1_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_IsoMu24_eta2p1_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_IsoMu24_eta2p1_tag_isLoaded = true;
		}
		return HLT_IsoMu24_eta2p1_tag_;
	}
	unsigned int &HLT_IsoMu24_eta2p1_probe()
	{
		if (not HLT_IsoMu24_eta2p1_probe_isLoaded) {
			if (HLT_IsoMu24_eta2p1_probe_branch != 0) {
				HLT_IsoMu24_eta2p1_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_IsoMu24_eta2p1_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_IsoMu24_eta2p1_probe_isLoaded = true;
		}
		return HLT_IsoMu24_eta2p1_probe_;
	}
	unsigned int &HLT_IsoMu24_eta2p1_version()
	{
		if (not HLT_IsoMu24_eta2p1_version_isLoaded) {
			if (HLT_IsoMu24_eta2p1_version_branch != 0) {
				HLT_IsoMu24_eta2p1_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_IsoMu24_eta2p1_version_branch does not exist!\n");
				exit(1);
			}
			HLT_IsoMu24_eta2p1_version_isLoaded = true;
		}
		return HLT_IsoMu24_eta2p1_version_;
	}
	unsigned int &HLT_Mu8_tag()
	{
		if (not HLT_Mu8_tag_isLoaded) {
			if (HLT_Mu8_tag_branch != 0) {
				HLT_Mu8_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu8_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu8_tag_isLoaded = true;
		}
		return HLT_Mu8_tag_;
	}
	unsigned int &HLT_Mu8_probe()
	{
		if (not HLT_Mu8_probe_isLoaded) {
			if (HLT_Mu8_probe_branch != 0) {
				HLT_Mu8_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu8_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu8_probe_isLoaded = true;
		}
		return HLT_Mu8_probe_;
	}
	unsigned int &HLT_Mu8_version()
	{
		if (not HLT_Mu8_version_isLoaded) {
			if (HLT_Mu8_version_branch != 0) {
				HLT_Mu8_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu8_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu8_version_isLoaded = true;
		}
		return HLT_Mu8_version_;
	}
	unsigned int &HLT_Mu17_tag()
	{
		if (not HLT_Mu17_tag_isLoaded) {
			if (HLT_Mu17_tag_branch != 0) {
				HLT_Mu17_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_tag_isLoaded = true;
		}
		return HLT_Mu17_tag_;
	}
	unsigned int &HLT_Mu17_probe()
	{
		if (not HLT_Mu17_probe_isLoaded) {
			if (HLT_Mu17_probe_branch != 0) {
				HLT_Mu17_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_probe_isLoaded = true;
		}
		return HLT_Mu17_probe_;
	}
	unsigned int &HLT_Mu17_version()
	{
		if (not HLT_Mu17_version_isLoaded) {
			if (HLT_Mu17_version_branch != 0) {
				HLT_Mu17_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Mu17_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Mu17_version_isLoaded = true;
		}
		return HLT_Mu17_version_;
	}
	unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_tag()
	{
		if (not HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_isLoaded) {
			if (HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch != 0) {
				HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_isLoaded = true;
		}
		return HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_;
	}
	unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_probe()
	{
		if (not HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_isLoaded) {
			if (HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch != 0) {
				HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_isLoaded = true;
		}
		return HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_;
	}
	unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_version()
	{
		if (not HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded) {
			if (HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch != 0) {
				HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded = true;
		}
		return HLT_Ele17_Ele8_L1sL1DoubleEG137_version_;
	}
	unsigned int &HLT_Ele17_Ele8_LeadingLeg_tag()
	{
		if (not HLT_Ele17_Ele8_LeadingLeg_tag_isLoaded) {
			if (HLT_Ele17_Ele8_LeadingLeg_tag_branch != 0) {
				HLT_Ele17_Ele8_LeadingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_LeadingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_LeadingLeg_tag_isLoaded = true;
		}
		return HLT_Ele17_Ele8_LeadingLeg_tag_;
	}
	unsigned int &HLT_Ele17_Ele8_LeadingLeg_probe()
	{
		if (not HLT_Ele17_Ele8_LeadingLeg_probe_isLoaded) {
			if (HLT_Ele17_Ele8_LeadingLeg_probe_branch != 0) {
				HLT_Ele17_Ele8_LeadingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_LeadingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_LeadingLeg_probe_isLoaded = true;
		}
		return HLT_Ele17_Ele8_LeadingLeg_probe_;
	}
	unsigned int &HLT_Ele17_Ele8_LeadingLeg_version()
	{
		if (not HLT_Ele17_Ele8_LeadingLeg_version_isLoaded) {
			if (HLT_Ele17_Ele8_LeadingLeg_version_branch != 0) {
				HLT_Ele17_Ele8_LeadingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_LeadingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_LeadingLeg_version_isLoaded = true;
		}
		return HLT_Ele17_Ele8_LeadingLeg_version_;
	}
	unsigned int &HLT_Ele17_Ele8_TrailingLeg_tag()
	{
		if (not HLT_Ele17_Ele8_TrailingLeg_tag_isLoaded) {
			if (HLT_Ele17_Ele8_TrailingLeg_tag_branch != 0) {
				HLT_Ele17_Ele8_TrailingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_TrailingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_TrailingLeg_tag_isLoaded = true;
		}
		return HLT_Ele17_Ele8_TrailingLeg_tag_;
	}
	unsigned int &HLT_Ele17_Ele8_TrailingLeg_probe()
	{
		if (not HLT_Ele17_Ele8_TrailingLeg_probe_isLoaded) {
			if (HLT_Ele17_Ele8_TrailingLeg_probe_branch != 0) {
				HLT_Ele17_Ele8_TrailingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_TrailingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_TrailingLeg_probe_isLoaded = true;
		}
		return HLT_Ele17_Ele8_TrailingLeg_probe_;
	}
	unsigned int &HLT_Ele17_Ele8_TrailingLeg_version()
	{
		if (not HLT_Ele17_Ele8_TrailingLeg_version_isLoaded) {
			if (HLT_Ele17_Ele8_TrailingLeg_version_branch != 0) {
				HLT_Ele17_Ele8_TrailingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_TrailingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_TrailingLeg_version_isLoaded = true;
		}
		return HLT_Ele17_Ele8_TrailingLeg_version_;
	}
	unsigned int &HLT_Ele17_Ele8_tag()
	{
		if (not HLT_Ele17_Ele8_tag_isLoaded) {
			if (HLT_Ele17_Ele8_tag_branch != 0) {
				HLT_Ele17_Ele8_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_tag_isLoaded = true;
		}
		return HLT_Ele17_Ele8_tag_;
	}
	unsigned int &HLT_Ele17_Ele8_probe()
	{
		if (not HLT_Ele17_Ele8_probe_isLoaded) {
			if (HLT_Ele17_Ele8_probe_branch != 0) {
				HLT_Ele17_Ele8_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_probe_isLoaded = true;
		}
		return HLT_Ele17_Ele8_probe_;
	}
	unsigned int &HLT_Ele17_Ele8_version()
	{
		if (not HLT_Ele17_Ele8_version_isLoaded) {
			if (HLT_Ele17_Ele8_version_branch != 0) {
				HLT_Ele17_Ele8_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_version_isLoaded = true;
		}
		return HLT_Ele17_Ele8_version_;
	}
	unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag()
	{
		if (not HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_isLoaded) {
			if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch != 0) {
				HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_isLoaded = true;
		}
		return HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_;
	}
	unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe()
	{
		if (not HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_isLoaded) {
			if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch != 0) {
				HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_isLoaded = true;
		}
		return HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_;
	}
	unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version()
	{
		if (not HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded) {
			if (HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch != 0) {
				HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded = true;
		}
		return HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_;
	}
	unsigned int &HLT_Ele27_WP80_tag()
	{
		if (not HLT_Ele27_WP80_tag_isLoaded) {
			if (HLT_Ele27_WP80_tag_branch != 0) {
				HLT_Ele27_WP80_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele27_WP80_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele27_WP80_tag_isLoaded = true;
		}
		return HLT_Ele27_WP80_tag_;
	}
	unsigned int &HLT_Ele27_WP80_probe()
	{
		if (not HLT_Ele27_WP80_probe_isLoaded) {
			if (HLT_Ele27_WP80_probe_branch != 0) {
				HLT_Ele27_WP80_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele27_WP80_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele27_WP80_probe_isLoaded = true;
		}
		return HLT_Ele27_WP80_probe_;
	}
	unsigned int &HLT_Ele27_WP80_version()
	{
		if (not HLT_Ele27_WP80_version_isLoaded) {
			if (HLT_Ele27_WP80_version_branch != 0) {
				HLT_Ele27_WP80_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele27_WP80_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele27_WP80_version_isLoaded = true;
		}
		return HLT_Ele27_WP80_version_;
	}
	unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_tag()
	{
		if (not HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_isLoaded) {
			if (HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch != 0) {
				HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_isLoaded = true;
		}
		return HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_;
	}
	unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_probe()
	{
		if (not HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_isLoaded) {
			if (HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch != 0) {
				HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_isLoaded = true;
		}
		return HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_;
	}
	unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_version()
	{
		if (not HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded) {
			if (HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch != 0) {
				HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded = true;
		}
		return HLT_Ele17_Ele8_Mass50_LeadingLeg_version_;
	}
	unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_tag()
	{
		if (not HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_isLoaded) {
			if (HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch != 0) {
				HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_isLoaded = true;
		}
		return HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_;
	}
	unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_probe()
	{
		if (not HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_isLoaded) {
			if (HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch != 0) {
				HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_isLoaded = true;
		}
		return HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_;
	}
	unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_version()
	{
		if (not HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded) {
			if (HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch != 0) {
				HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded = true;
		}
		return HLT_Ele17_Ele8_Mass50_TrailingLeg_version_;
	}
	unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_tag()
	{
		if (not HLT_Ele20_SC4_Mass50_LeadingLeg_tag_isLoaded) {
			if (HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch != 0) {
				HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele20_SC4_Mass50_LeadingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele20_SC4_Mass50_LeadingLeg_tag_isLoaded = true;
		}
		return HLT_Ele20_SC4_Mass50_LeadingLeg_tag_;
	}
	unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_probe()
	{
		if (not HLT_Ele20_SC4_Mass50_LeadingLeg_probe_isLoaded) {
			if (HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch != 0) {
				HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele20_SC4_Mass50_LeadingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele20_SC4_Mass50_LeadingLeg_probe_isLoaded = true;
		}
		return HLT_Ele20_SC4_Mass50_LeadingLeg_probe_;
	}
	unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_version()
	{
		if (not HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded) {
			if (HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch != 0) {
				HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded = true;
		}
		return HLT_Ele20_SC4_Mass50_LeadingLeg_version_;
	}
	unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_tag()
	{
		if (not HLT_Ele20_SC4_Mass50_TrailingLeg_tag_isLoaded) {
			if (HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch != 0) {
				HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele20_SC4_Mass50_TrailingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele20_SC4_Mass50_TrailingLeg_tag_isLoaded = true;
		}
		return HLT_Ele20_SC4_Mass50_TrailingLeg_tag_;
	}
	unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_probe()
	{
		if (not HLT_Ele20_SC4_Mass50_TrailingLeg_probe_isLoaded) {
			if (HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch != 0) {
				HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele20_SC4_Mass50_TrailingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele20_SC4_Mass50_TrailingLeg_probe_isLoaded = true;
		}
		return HLT_Ele20_SC4_Mass50_TrailingLeg_probe_;
	}
	unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_version()
	{
		if (not HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded) {
			if (HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch != 0) {
				HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded = true;
		}
		return HLT_Ele20_SC4_Mass50_TrailingLeg_version_;
	}
	unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_tag()
	{
		if (not HLT_Ele32_SC17_Mass50_LeadingLeg_tag_isLoaded) {
			if (HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch != 0) {
				HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele32_SC17_Mass50_LeadingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele32_SC17_Mass50_LeadingLeg_tag_isLoaded = true;
		}
		return HLT_Ele32_SC17_Mass50_LeadingLeg_tag_;
	}
	unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_probe()
	{
		if (not HLT_Ele32_SC17_Mass50_LeadingLeg_probe_isLoaded) {
			if (HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch != 0) {
				HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele32_SC17_Mass50_LeadingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele32_SC17_Mass50_LeadingLeg_probe_isLoaded = true;
		}
		return HLT_Ele32_SC17_Mass50_LeadingLeg_probe_;
	}
	unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_version()
	{
		if (not HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded) {
			if (HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch != 0) {
				HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded = true;
		}
		return HLT_Ele32_SC17_Mass50_LeadingLeg_version_;
	}
	unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_tag()
	{
		if (not HLT_Ele32_SC17_Mass50_TrailingLeg_tag_isLoaded) {
			if (HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch != 0) {
				HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele32_SC17_Mass50_TrailingLeg_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele32_SC17_Mass50_TrailingLeg_tag_isLoaded = true;
		}
		return HLT_Ele32_SC17_Mass50_TrailingLeg_tag_;
	}
	unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_probe()
	{
		if (not HLT_Ele32_SC17_Mass50_TrailingLeg_probe_isLoaded) {
			if (HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch != 0) {
				HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele32_SC17_Mass50_TrailingLeg_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele32_SC17_Mass50_TrailingLeg_probe_isLoaded = true;
		}
		return HLT_Ele32_SC17_Mass50_TrailingLeg_probe_;
	}
	unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_version()
	{
		if (not HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded) {
			if (HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch != 0) {
				HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded = true;
		}
		return HLT_Ele32_SC17_Mass50_TrailingLeg_version_;
	}
	unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag()
	{
		if (not HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded) {
			if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch != 0) {
				HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_;
	}
	unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe()
	{
		if (not HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded) {
			if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch != 0) {
				HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_;
	}
	unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version()
	{
		if (not HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded) {
			if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch != 0) {
				HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;
	}
	unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag()
	{
		if (not HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded) {
			if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch != 0) {
				HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_;
	}
	unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe()
	{
		if (not HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded) {
			if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch != 0) {
				HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_;
	}
	unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version()
	{
		if (not HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded) {
			if (HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch != 0) {
				HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;
	}
	unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag()
	{
		if (not HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded) {
			if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch != 0) {
				HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_isLoaded = true;
		}
		return HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag_;
	}
	unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe()
	{
		if (not HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded) {
			if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch != 0) {
				HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_isLoaded = true;
		}
		return HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_;
	}
	unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version()
	{
		if (not HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded) {
			if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch != 0) {
				HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_isLoaded = true;
		}
		return HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_;
	}
	unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag()
	{
		if (not HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded) {
			if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch != 0) {
				HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_isLoaded = true;
		}
		return HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag_;
	}
	unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe()
	{
		if (not HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded) {
			if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch != 0) {
				HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_isLoaded = true;
		}
		return HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_;
	}
	unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version()
	{
		if (not HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded) {
			if (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch != 0) {
				HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_isLoaded = true;
		}
		return HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_;
	}
	unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_tag()
	{
		if (not HLT_Ele8_CaloIdT_TrkIdVL_tag_isLoaded) {
			if (HLT_Ele8_CaloIdT_TrkIdVL_tag_branch != 0) {
				HLT_Ele8_CaloIdT_TrkIdVL_tag_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_TrkIdVL_tag_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_TrkIdVL_tag_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_TrkIdVL_tag_;
	}
	unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_probe()
	{
		if (not HLT_Ele8_CaloIdT_TrkIdVL_probe_isLoaded) {
			if (HLT_Ele8_CaloIdT_TrkIdVL_probe_branch != 0) {
				HLT_Ele8_CaloIdT_TrkIdVL_probe_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_TrkIdVL_probe_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_TrkIdVL_probe_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_TrkIdVL_probe_;
	}
	unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_version()
	{
		if (not HLT_Ele8_CaloIdT_TrkIdVL_version_isLoaded) {
			if (HLT_Ele8_CaloIdT_TrkIdVL_version_branch != 0) {
				HLT_Ele8_CaloIdT_TrkIdVL_version_branch->GetEntry(index);
			} else { 
				printf("branch HLT_Ele8_CaloIdT_TrkIdVL_version_branch does not exist!\n");
				exit(1);
			}
			HLT_Ele8_CaloIdT_TrkIdVL_version_isLoaded = true;
		}
		return HLT_Ele8_CaloIdT_TrkIdVL_version_;
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern LeptonTree lepton_tree_obj;
#endif

namespace lepton_tree {
	const unsigned int &event();
	const unsigned int &run();
	const unsigned int &lumi();
	const float &rnd();
	const unsigned int &nvtx();
	const unsigned int &npu();
	const unsigned int &npuPlusOne();
	const unsigned int &npuMinusOne();
	const unsigned int &leptonSelection();
	const unsigned int &eventSelection();
	const float &rhoIsoAll();
	const float &rhoIsoAllCentral();
	const float &rhoIsoNeutral();
	const float &tagAndProbeMass();
	const bool &tagAndProbeIsRandom();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &tag();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &probe();
	const int &qTag();
	const int &qProbe();
	const float &scale1fb();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet1();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet2();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet3();
	const float &met();
	const float &metPhi();
	const float &trackMet();
	const float &trackMetPhi();
	const unsigned int &njets();
	const unsigned int &nleps();
	const float &hltPrescale();
	const float &sumet();
	const float &metSig();
	const float &mt();
	const float &dPhiProbeJet1();
	const float &chargedEmFracJet1();
	const float &neutralEmFracJet1();
	const float &chargedMuFracJet1();
	const float &electronHWW2011MVA();
	const float &egammaPOG2012MVA();
	const float &muonHZZ2012IsoRingsMVA();
	const unsigned int &vetoId();
	const unsigned int &looseId();
	const unsigned int &mediumId();
	const unsigned int &tightId();
	const float &pfmva();
	const float &sceta();
	const float &scphi();
	const float &scenergy();
	const bool &chargesAgree();
	const float &eopin();
	const float &ooemoop();
	const float &fbrem();
	const float &detain();
	const float &dphiin();
	const float &hoe();
	const float &hoetow();
	const float &sieie();
	const float &d0vtx();
	const float &dzvtx();
	const bool &vfitprob();
	const float &mhit();
	const float &ecaliso();
	const float &hcaliso();
	const float &trkiso();
	const float &pfemiso03();
	const float &pfchiso03();
	const float &pfnhiso03();
	const float &pfemiso04();
	const float &pfchiso04();
	const float &pfnhiso04();
	const float &radiso03();
	const float &radiso04();
	const float &iso2011();
	const float &ea04();
	const float &ea03();
	const float &dbeta03();
	const float &dbeta04();
	const float &el_test_pfchiso04_trkveto();
	const float &el_test_pfchiso04_dzcut();
	const float &el_test_pfchiso04_ebveto();
	const float &el_test_pfemiso04_ebveto();
	const float &eaem04();
	const float &eanh04();
	const float &gen_drs1();
	const float &gen_drs3();
	const unsigned int &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly();
	const unsigned int &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version();
	const unsigned int &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly();
	const unsigned int &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version();
	const unsigned int &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly();
	const unsigned int &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version();
	const unsigned int &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly();
	const unsigned int &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version();
	const unsigned int &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly();
	const unsigned int &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version();
	const unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag();
	const unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe();
	const unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version();
	const unsigned int &HLT_Mu17_Mu8_TrailingLeg_tag();
	const unsigned int &HLT_Mu17_Mu8_TrailingLeg_probe();
	const unsigned int &HLT_Mu17_Mu8_TrailingLeg_version();
	const unsigned int &HLT_Mu17_Mu8_LeadingLeg_tag();
	const unsigned int &HLT_Mu17_Mu8_LeadingLeg_probe();
	const unsigned int &HLT_Mu17_Mu8_LeadingLeg_version();
	const unsigned int &HLT_Mu17_Mu8_tag();
	const unsigned int &HLT_Mu17_Mu8_probe();
	const unsigned int &HLT_Mu17_Mu8_version();
	const unsigned int &HLT_Mu17_TkMu8_TrailingLeg_tag();
	const unsigned int &HLT_Mu17_TkMu8_TrailingLeg_probe();
	const unsigned int &HLT_Mu17_TkMu8_TrailingLeg_version();
	const unsigned int &HLT_Mu17_TkMu8_LeadingLeg_tag();
	const unsigned int &HLT_Mu17_TkMu8_LeadingLeg_probe();
	const unsigned int &HLT_Mu17_TkMu8_LeadingLeg_version();
	const unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag();
	const unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe();
	const unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version();
	const unsigned int &HLT_Mu17_TkMu8_tag();
	const unsigned int &HLT_Mu17_TkMu8_probe();
	const unsigned int &HLT_Mu17_TkMu8_version();
	const unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag();
	const unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe();
	const unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version();
	const unsigned int &HLT_IsoMu24_eta2p1_tag();
	const unsigned int &HLT_IsoMu24_eta2p1_probe();
	const unsigned int &HLT_IsoMu24_eta2p1_version();
	const unsigned int &HLT_Mu8_tag();
	const unsigned int &HLT_Mu8_probe();
	const unsigned int &HLT_Mu8_version();
	const unsigned int &HLT_Mu17_tag();
	const unsigned int &HLT_Mu17_probe();
	const unsigned int &HLT_Mu17_version();
	const unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_tag();
	const unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_probe();
	const unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_version();
	const unsigned int &HLT_Ele17_Ele8_LeadingLeg_tag();
	const unsigned int &HLT_Ele17_Ele8_LeadingLeg_probe();
	const unsigned int &HLT_Ele17_Ele8_LeadingLeg_version();
	const unsigned int &HLT_Ele17_Ele8_TrailingLeg_tag();
	const unsigned int &HLT_Ele17_Ele8_TrailingLeg_probe();
	const unsigned int &HLT_Ele17_Ele8_TrailingLeg_version();
	const unsigned int &HLT_Ele17_Ele8_tag();
	const unsigned int &HLT_Ele17_Ele8_probe();
	const unsigned int &HLT_Ele17_Ele8_version();
	const unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag();
	const unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe();
	const unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version();
	const unsigned int &HLT_Ele27_WP80_tag();
	const unsigned int &HLT_Ele27_WP80_probe();
	const unsigned int &HLT_Ele27_WP80_version();
	const unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_tag();
	const unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_probe();
	const unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_version();
	const unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_tag();
	const unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_probe();
	const unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_version();
	const unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_tag();
	const unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_probe();
	const unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_version();
	const unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_tag();
	const unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_probe();
	const unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_version();
	const unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_tag();
	const unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_probe();
	const unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_version();
	const unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_tag();
	const unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_probe();
	const unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_version();
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag();
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe();
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version();
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag();
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe();
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version();
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag();
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe();
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version();
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag();
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe();
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version();
	const unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_tag();
	const unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_probe();
	const unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_version();
}
#endif
