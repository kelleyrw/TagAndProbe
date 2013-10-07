#include "TagAndProbe/Analysis/interface/LeptonTree.h"
LeptonTree lepton_tree_obj;
namespace lepton_tree {
	const unsigned int &event() { return lepton_tree_obj.event(); }
	const unsigned int &run() { return lepton_tree_obj.run(); }
	const unsigned int &lumi() { return lepton_tree_obj.lumi(); }
	const float &rnd() { return lepton_tree_obj.rnd(); }
	const unsigned int &nvtx() { return lepton_tree_obj.nvtx(); }
	const unsigned int &npu() { return lepton_tree_obj.npu(); }
	const unsigned int &npuPlusOne() { return lepton_tree_obj.npuPlusOne(); }
	const unsigned int &npuMinusOne() { return lepton_tree_obj.npuMinusOne(); }
	const unsigned int &leptonSelection() { return lepton_tree_obj.leptonSelection(); }
	const unsigned int &eventSelection() { return lepton_tree_obj.eventSelection(); }
	const float &rhoIsoAll() { return lepton_tree_obj.rhoIsoAll(); }
	const float &rhoIsoAllCentral() { return lepton_tree_obj.rhoIsoAllCentral(); }
	const float &rhoIsoNeutral() { return lepton_tree_obj.rhoIsoNeutral(); }
	const float &tagAndProbeMass() { return lepton_tree_obj.tagAndProbeMass(); }
	const bool &tagAndProbeIsRandom() { return lepton_tree_obj.tagAndProbeIsRandom(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &tag() { return lepton_tree_obj.tag(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &probe() { return lepton_tree_obj.probe(); }
	const int &qTag() { return lepton_tree_obj.qTag(); }
	const int &qProbe() { return lepton_tree_obj.qProbe(); }
	const float &scale1fb() { return lepton_tree_obj.scale1fb(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet1() { return lepton_tree_obj.jet1(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet2() { return lepton_tree_obj.jet2(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &jet3() { return lepton_tree_obj.jet3(); }
	const float &met() { return lepton_tree_obj.met(); }
	const float &metPhi() { return lepton_tree_obj.metPhi(); }
	const float &trackMet() { return lepton_tree_obj.trackMet(); }
	const float &trackMetPhi() { return lepton_tree_obj.trackMetPhi(); }
	const unsigned int &njets() { return lepton_tree_obj.njets(); }
	const unsigned int &nleps() { return lepton_tree_obj.nleps(); }
	const float &hltPrescale() { return lepton_tree_obj.hltPrescale(); }
	const float &sumet() { return lepton_tree_obj.sumet(); }
	const float &metSig() { return lepton_tree_obj.metSig(); }
	const float &mt() { return lepton_tree_obj.mt(); }
	const float &dPhiProbeJet1() { return lepton_tree_obj.dPhiProbeJet1(); }
	const float &chargedEmFracJet1() { return lepton_tree_obj.chargedEmFracJet1(); }
	const float &neutralEmFracJet1() { return lepton_tree_obj.neutralEmFracJet1(); }
	const float &chargedMuFracJet1() { return lepton_tree_obj.chargedMuFracJet1(); }
	const float &electronHWW2011MVA() { return lepton_tree_obj.electronHWW2011MVA(); }
	const float &egammaPOG2012MVA() { return lepton_tree_obj.egammaPOG2012MVA(); }
	const float &muonHZZ2012IsoRingsMVA() { return lepton_tree_obj.muonHZZ2012IsoRingsMVA(); }
	const unsigned int &vetoId() { return lepton_tree_obj.vetoId(); }
	const unsigned int &looseId() { return lepton_tree_obj.looseId(); }
	const unsigned int &mediumId() { return lepton_tree_obj.mediumId(); }
	const unsigned int &tightId() { return lepton_tree_obj.tightId(); }
	const float &pfmva() { return lepton_tree_obj.pfmva(); }
	const float &sceta() { return lepton_tree_obj.sceta(); }
	const float &scphi() { return lepton_tree_obj.scphi(); }
	const float &scenergy() { return lepton_tree_obj.scenergy(); }
	const bool &chargesAgree() { return lepton_tree_obj.chargesAgree(); }
	const float &eopin() { return lepton_tree_obj.eopin(); }
	const float &ooemoop() { return lepton_tree_obj.ooemoop(); }
	const float &fbrem() { return lepton_tree_obj.fbrem(); }
	const float &detain() { return lepton_tree_obj.detain(); }
	const float &dphiin() { return lepton_tree_obj.dphiin(); }
	const float &hoe() { return lepton_tree_obj.hoe(); }
	const float &hoetow() { return lepton_tree_obj.hoetow(); }
	const float &sieie() { return lepton_tree_obj.sieie(); }
	const float &d0vtx() { return lepton_tree_obj.d0vtx(); }
	const float &dzvtx() { return lepton_tree_obj.dzvtx(); }
	const bool &vfitprob() { return lepton_tree_obj.vfitprob(); }
	const float &mhit() { return lepton_tree_obj.mhit(); }
	const float &ecaliso() { return lepton_tree_obj.ecaliso(); }
	const float &hcaliso() { return lepton_tree_obj.hcaliso(); }
	const float &trkiso() { return lepton_tree_obj.trkiso(); }
	const float &pfemiso03() { return lepton_tree_obj.pfemiso03(); }
	const float &pfchiso03() { return lepton_tree_obj.pfchiso03(); }
	const float &pfnhiso03() { return lepton_tree_obj.pfnhiso03(); }
	const float &pfemiso04() { return lepton_tree_obj.pfemiso04(); }
	const float &pfchiso04() { return lepton_tree_obj.pfchiso04(); }
	const float &pfnhiso04() { return lepton_tree_obj.pfnhiso04(); }
	const float &radiso03() { return lepton_tree_obj.radiso03(); }
	const float &radiso04() { return lepton_tree_obj.radiso04(); }
	const float &iso2011() { return lepton_tree_obj.iso2011(); }
	const float &ea04() { return lepton_tree_obj.ea04(); }
	const float &ea03() { return lepton_tree_obj.ea03(); }
	const float &dbeta03() { return lepton_tree_obj.dbeta03(); }
	const float &dbeta04() { return lepton_tree_obj.dbeta04(); }
	const float &el_test_pfchiso04_trkveto() { return lepton_tree_obj.el_test_pfchiso04_trkveto(); }
	const float &el_test_pfchiso04_dzcut() { return lepton_tree_obj.el_test_pfchiso04_dzcut(); }
	const float &el_test_pfchiso04_ebveto() { return lepton_tree_obj.el_test_pfchiso04_ebveto(); }
	const float &el_test_pfemiso04_ebveto() { return lepton_tree_obj.el_test_pfemiso04_ebveto(); }
	const float &eaem04() { return lepton_tree_obj.eaem04(); }
	const float &eanh04() { return lepton_tree_obj.eanh04(); }
	const float &gen_drs1() { return lepton_tree_obj.gen_drs1(); }
	const float &gen_drs3() { return lepton_tree_obj.gen_drs3(); }
	const unsigned int &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly() { return lepton_tree_obj.HLT_Photon22_R9Id90_HE10_Iso40_EBOnly(); }
	const unsigned int &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version() { return lepton_tree_obj.HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_version(); }
	const unsigned int &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly() { return lepton_tree_obj.HLT_Photon36_R9Id90_HE10_Iso40_EBOnly(); }
	const unsigned int &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version() { return lepton_tree_obj.HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_version(); }
	const unsigned int &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly() { return lepton_tree_obj.HLT_Photon50_R9Id90_HE10_Iso40_EBOnly(); }
	const unsigned int &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version() { return lepton_tree_obj.HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_version(); }
	const unsigned int &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly() { return lepton_tree_obj.HLT_Photon75_R9Id90_HE10_Iso40_EBOnly(); }
	const unsigned int &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version() { return lepton_tree_obj.HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_version(); }
	const unsigned int &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly() { return lepton_tree_obj.HLT_Photon90_R9Id90_HE10_Iso40_EBOnly(); }
	const unsigned int &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version() { return lepton_tree_obj.HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_version(); }
	const unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag() { return lepton_tree_obj.HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag(); }
	const unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe() { return lepton_tree_obj.HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe(); }
	const unsigned int &HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version() { return lepton_tree_obj.HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version(); }
	const unsigned int &HLT_Mu17_Mu8_TrailingLeg_tag() { return lepton_tree_obj.HLT_Mu17_Mu8_TrailingLeg_tag(); }
	const unsigned int &HLT_Mu17_Mu8_TrailingLeg_probe() { return lepton_tree_obj.HLT_Mu17_Mu8_TrailingLeg_probe(); }
	const unsigned int &HLT_Mu17_Mu8_TrailingLeg_version() { return lepton_tree_obj.HLT_Mu17_Mu8_TrailingLeg_version(); }
	const unsigned int &HLT_Mu17_Mu8_LeadingLeg_tag() { return lepton_tree_obj.HLT_Mu17_Mu8_LeadingLeg_tag(); }
	const unsigned int &HLT_Mu17_Mu8_LeadingLeg_probe() { return lepton_tree_obj.HLT_Mu17_Mu8_LeadingLeg_probe(); }
	const unsigned int &HLT_Mu17_Mu8_LeadingLeg_version() { return lepton_tree_obj.HLT_Mu17_Mu8_LeadingLeg_version(); }
	const unsigned int &HLT_Mu17_Mu8_tag() { return lepton_tree_obj.HLT_Mu17_Mu8_tag(); }
	const unsigned int &HLT_Mu17_Mu8_probe() { return lepton_tree_obj.HLT_Mu17_Mu8_probe(); }
	const unsigned int &HLT_Mu17_Mu8_version() { return lepton_tree_obj.HLT_Mu17_Mu8_version(); }
	const unsigned int &HLT_Mu17_TkMu8_TrailingLeg_tag() { return lepton_tree_obj.HLT_Mu17_TkMu8_TrailingLeg_tag(); }
	const unsigned int &HLT_Mu17_TkMu8_TrailingLeg_probe() { return lepton_tree_obj.HLT_Mu17_TkMu8_TrailingLeg_probe(); }
	const unsigned int &HLT_Mu17_TkMu8_TrailingLeg_version() { return lepton_tree_obj.HLT_Mu17_TkMu8_TrailingLeg_version(); }
	const unsigned int &HLT_Mu17_TkMu8_LeadingLeg_tag() { return lepton_tree_obj.HLT_Mu17_TkMu8_LeadingLeg_tag(); }
	const unsigned int &HLT_Mu17_TkMu8_LeadingLeg_probe() { return lepton_tree_obj.HLT_Mu17_TkMu8_LeadingLeg_probe(); }
	const unsigned int &HLT_Mu17_TkMu8_LeadingLeg_version() { return lepton_tree_obj.HLT_Mu17_TkMu8_LeadingLeg_version(); }
	const unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag() { return lepton_tree_obj.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag(); }
	const unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe() { return lepton_tree_obj.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe(); }
	const unsigned int &HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version() { return lepton_tree_obj.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version(); }
	const unsigned int &HLT_Mu17_TkMu8_tag() { return lepton_tree_obj.HLT_Mu17_TkMu8_tag(); }
	const unsigned int &HLT_Mu17_TkMu8_probe() { return lepton_tree_obj.HLT_Mu17_TkMu8_probe(); }
	const unsigned int &HLT_Mu17_TkMu8_version() { return lepton_tree_obj.HLT_Mu17_TkMu8_version(); }
	const unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag() { return lepton_tree_obj.HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag(); }
	const unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe() { return lepton_tree_obj.HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe(); }
	const unsigned int &HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version() { return lepton_tree_obj.HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version(); }
	const unsigned int &HLT_IsoMu24_eta2p1_tag() { return lepton_tree_obj.HLT_IsoMu24_eta2p1_tag(); }
	const unsigned int &HLT_IsoMu24_eta2p1_probe() { return lepton_tree_obj.HLT_IsoMu24_eta2p1_probe(); }
	const unsigned int &HLT_IsoMu24_eta2p1_version() { return lepton_tree_obj.HLT_IsoMu24_eta2p1_version(); }
	const unsigned int &HLT_Mu8_tag() { return lepton_tree_obj.HLT_Mu8_tag(); }
	const unsigned int &HLT_Mu8_probe() { return lepton_tree_obj.HLT_Mu8_probe(); }
	const unsigned int &HLT_Mu8_version() { return lepton_tree_obj.HLT_Mu8_version(); }
	const unsigned int &HLT_Mu17_tag() { return lepton_tree_obj.HLT_Mu17_tag(); }
	const unsigned int &HLT_Mu17_probe() { return lepton_tree_obj.HLT_Mu17_probe(); }
	const unsigned int &HLT_Mu17_version() { return lepton_tree_obj.HLT_Mu17_version(); }
	const unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_tag() { return lepton_tree_obj.HLT_Ele17_Ele8_L1sL1DoubleEG137_tag(); }
	const unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_probe() { return lepton_tree_obj.HLT_Ele17_Ele8_L1sL1DoubleEG137_probe(); }
	const unsigned int &HLT_Ele17_Ele8_L1sL1DoubleEG137_version() { return lepton_tree_obj.HLT_Ele17_Ele8_L1sL1DoubleEG137_version(); }
	const unsigned int &HLT_Ele17_Ele8_LeadingLeg_tag() { return lepton_tree_obj.HLT_Ele17_Ele8_LeadingLeg_tag(); }
	const unsigned int &HLT_Ele17_Ele8_LeadingLeg_probe() { return lepton_tree_obj.HLT_Ele17_Ele8_LeadingLeg_probe(); }
	const unsigned int &HLT_Ele17_Ele8_LeadingLeg_version() { return lepton_tree_obj.HLT_Ele17_Ele8_LeadingLeg_version(); }
	const unsigned int &HLT_Ele17_Ele8_TrailingLeg_tag() { return lepton_tree_obj.HLT_Ele17_Ele8_TrailingLeg_tag(); }
	const unsigned int &HLT_Ele17_Ele8_TrailingLeg_probe() { return lepton_tree_obj.HLT_Ele17_Ele8_TrailingLeg_probe(); }
	const unsigned int &HLT_Ele17_Ele8_TrailingLeg_version() { return lepton_tree_obj.HLT_Ele17_Ele8_TrailingLeg_version(); }
	const unsigned int &HLT_Ele17_Ele8_tag() { return lepton_tree_obj.HLT_Ele17_Ele8_tag(); }
	const unsigned int &HLT_Ele17_Ele8_probe() { return lepton_tree_obj.HLT_Ele17_Ele8_probe(); }
	const unsigned int &HLT_Ele17_Ele8_version() { return lepton_tree_obj.HLT_Ele17_Ele8_version(); }
	const unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag() { return lepton_tree_obj.HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag(); }
	const unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe() { return lepton_tree_obj.HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe(); }
	const unsigned int &HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version() { return lepton_tree_obj.HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version(); }
	const unsigned int &HLT_Ele27_WP80_tag() { return lepton_tree_obj.HLT_Ele27_WP80_tag(); }
	const unsigned int &HLT_Ele27_WP80_probe() { return lepton_tree_obj.HLT_Ele27_WP80_probe(); }
	const unsigned int &HLT_Ele27_WP80_version() { return lepton_tree_obj.HLT_Ele27_WP80_version(); }
	const unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_tag() { return lepton_tree_obj.HLT_Ele17_Ele8_Mass50_LeadingLeg_tag(); }
	const unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_probe() { return lepton_tree_obj.HLT_Ele17_Ele8_Mass50_LeadingLeg_probe(); }
	const unsigned int &HLT_Ele17_Ele8_Mass50_LeadingLeg_version() { return lepton_tree_obj.HLT_Ele17_Ele8_Mass50_LeadingLeg_version(); }
	const unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_tag() { return lepton_tree_obj.HLT_Ele17_Ele8_Mass50_TrailingLeg_tag(); }
	const unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_probe() { return lepton_tree_obj.HLT_Ele17_Ele8_Mass50_TrailingLeg_probe(); }
	const unsigned int &HLT_Ele17_Ele8_Mass50_TrailingLeg_version() { return lepton_tree_obj.HLT_Ele17_Ele8_Mass50_TrailingLeg_version(); }
	const unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_tag() { return lepton_tree_obj.HLT_Ele20_SC4_Mass50_LeadingLeg_tag(); }
	const unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_probe() { return lepton_tree_obj.HLT_Ele20_SC4_Mass50_LeadingLeg_probe(); }
	const unsigned int &HLT_Ele20_SC4_Mass50_LeadingLeg_version() { return lepton_tree_obj.HLT_Ele20_SC4_Mass50_LeadingLeg_version(); }
	const unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_tag() { return lepton_tree_obj.HLT_Ele20_SC4_Mass50_TrailingLeg_tag(); }
	const unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_probe() { return lepton_tree_obj.HLT_Ele20_SC4_Mass50_TrailingLeg_probe(); }
	const unsigned int &HLT_Ele20_SC4_Mass50_TrailingLeg_version() { return lepton_tree_obj.HLT_Ele20_SC4_Mass50_TrailingLeg_version(); }
	const unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_tag() { return lepton_tree_obj.HLT_Ele32_SC17_Mass50_LeadingLeg_tag(); }
	const unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_probe() { return lepton_tree_obj.HLT_Ele32_SC17_Mass50_LeadingLeg_probe(); }
	const unsigned int &HLT_Ele32_SC17_Mass50_LeadingLeg_version() { return lepton_tree_obj.HLT_Ele32_SC17_Mass50_LeadingLeg_version(); }
	const unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_tag() { return lepton_tree_obj.HLT_Ele32_SC17_Mass50_TrailingLeg_tag(); }
	const unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_probe() { return lepton_tree_obj.HLT_Ele32_SC17_Mass50_TrailingLeg_probe(); }
	const unsigned int &HLT_Ele32_SC17_Mass50_TrailingLeg_version() { return lepton_tree_obj.HLT_Ele32_SC17_Mass50_TrailingLeg_version(); }
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag() { return lepton_tree_obj.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag(); }
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe() { return lepton_tree_obj.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe(); }
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version() { return lepton_tree_obj.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version(); }
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag() { return lepton_tree_obj.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag(); }
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe() { return lepton_tree_obj.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe(); }
	const unsigned int &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version() { return lepton_tree_obj.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version(); }
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag() { return lepton_tree_obj.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_tag(); }
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe() { return lepton_tree_obj.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe(); }
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version() { return lepton_tree_obj.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version(); }
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag() { return lepton_tree_obj.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_tag(); }
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe() { return lepton_tree_obj.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe(); }
	const unsigned int &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version() { return lepton_tree_obj.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version(); }
	const unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_tag() { return lepton_tree_obj.HLT_Ele8_CaloIdT_TrkIdVL_tag(); }
	const unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_probe() { return lepton_tree_obj.HLT_Ele8_CaloIdT_TrkIdVL_probe(); }
	const unsigned int &HLT_Ele8_CaloIdT_TrkIdVL_version() { return lepton_tree_obj.HLT_Ele8_CaloIdT_TrkIdVL_version(); }
}
