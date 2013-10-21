#include "TagAndProbe/Analysis/interface/Measurement.h"

#include <iostream>

using namespace tnp;
using namespace std;

// Tools
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "TagAndProbe/Analysis/interface/LeptonTree.h"
#include "TagAndProbe/Analysis/interface/LeptonSelections.h"

namespace tnp
{
    Lepton::value_type GetLeptonFromString(const std::string& lepton_name)
    {
        if (lt::string_lower(lepton_name) == "muon"    ) {return Lepton::Muon;    }
        if (lt::string_lower(lepton_name) == "electron") {return Lepton::Electron;}
        throw std::invalid_argument("[tnp::GetLeptonFromString]: ERROR - invalid value!"); 
    }

    std::string GetStringFromLepton(const Lepton::value_type lepton_type)
    {
        if (lepton_type == Lepton::Muon    ) return "muon";
        if (lepton_type == Lepton::Electron) return "electron";
        throw std::invalid_argument(Form("[tnp::GetStringFromLepton]: ERROR - invalid value! %u", lepton_type)); 
    }

    Selection::value_type GetSelectionFromString(const std::string& sel_name)
    {
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumWPDenID"  )) {return Selection::EGammaMediumWPDenID;  } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumWPDenIso" )) {return Selection::EGammaMediumWPDenIso; } 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumWPDenBoth")) {return Selection::EGammaMediumWPDenBoth;} 
        if (lt::string_lower(sel_name) == lt::string_lower("EGammaMediumWPNum"    )) {return Selection::EGammaMediumWPNum;    } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenID"       )) {return Selection::MuTightWPDenID;       } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenIso"      )) {return Selection::MuTightWPDenIso;      } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPDenBoth"     )) {return Selection::MuTightWPDenBoth;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("MuTightWPNum"         )) {return Selection::MuTightWPNum;         } 
        if (lt::string_lower(sel_name) == lt::string_lower("SameSignDenID"        )) {return Selection::SameSignDenID;      } 
        if (lt::string_lower(sel_name) == lt::string_lower("SameSignDenIso"       )) {return Selection::SameSignDenIso;     } 
        if (lt::string_lower(sel_name) == lt::string_lower("SameSignDenBoth"      )) {return Selection::SameSignDenBoth;    } 
        if (lt::string_lower(sel_name) == lt::string_lower("SameSignNum"          )) {return Selection::SameSignNum;        } 
        throw std::invalid_argument("[tnp::GetSelectionFromString]: ERROR - invalid value!"); 
    }

    std::string GetStringFromSelection(const Selection::value_type sel_type)
    {
        if (sel_type == Selection::EGammaMediumWPDenID  ) return "EGammaMediumWPDenID";
        if (sel_type == Selection::EGammaMediumWPDenIso ) return "EGammaMediumWPDenIso";
        if (sel_type == Selection::EGammaMediumWPDenBoth) return "EGammaMediumWPDenBoth";
        if (sel_type == Selection::EGammaMediumWPNum    ) return "EGammaMediumWPNum";
        if (sel_type == Selection::MuTightWPDenID       ) return "MuTightWPDenID";
        if (sel_type == Selection::MuTightWPDenIso      ) return "MuTightWPDenIso";
        if (sel_type == Selection::MuTightWPDenBoth     ) return "MuTightWPDenBoth";
        if (sel_type == Selection::MuTightWPNum         ) return "MuTightWPNum";
        if (sel_type == Selection::SameSignDenID        ) return "SameSignDenID";
        if (sel_type == Selection::SameSignDenIso       ) return "SameSignDenIso";
        if (sel_type == Selection::SameSignDenBoth      ) return "SameSignDenBoth";
        if (sel_type == Selection::SameSignNum          ) return "SameSignNum";
        throw std::invalid_argument("[tnp::GetStringFromSelection]: ERROR - invalid value!"); 
    }

    // passes selection based on above enum
    bool PassesSelection(const Lepton::value_type lepton_type, const Selection::value_type selection, const bool is_data)
    {
        using namespace lepton_tree;

        // --------------------------------------------------------------------------- //
        // electrons
        // --------------------------------------------------------------------------- //

        if (lepton_type == Lepton::Electron)
        {
            // cut values and variables
            const float el_is_barrel   = fabs(sceta()) < 1.4442;
            const float el_is_endcap   = fabs(sceta()) > 1.566;
            const float el_is_crack    = not (el_is_barrel or el_is_endcap);
            const float el_probe_pt    = probe().pt();
            const float el_tag_pt      = tag().pt();
            const float el_d0          = fabs(d0vtx()); 
            const float el_iso         = (pfchiso03() + TMath::Max(0.0f, pfemiso03() + pfnhiso03() - ea03() * TMath::Max(0.0f, rhoIsoAllCentral())))/el_probe_pt;
            const float el_iso_eg_cut  = (el_is_endcap ? (el_probe_pt < 20.0 ? 0.10 : 0.15) : 0.15);  // egammamediumwp value 
            const float el_iso_ss_cut  = 0.09;                                                        // ss2012 value
            const float el_d0_cut      = 0.010;
            const float el_tag_pt_cut  = 32.0;
            const float el_hoe_cut     = (el_is_endcap ? 0.075 : 0.1); 
            const int   el_mhits_cut   = 0; // allow maximum number of missing hits
            const bool  el_3q          = chargesAgree();

            // cut decisions 
            const bool el_passes_pt       = (el_tag_pt > el_tag_pt_cut);
            const bool el_passes_trig_tag = (is_data ? HLT_Ele27_WP80_tag() != 0 : true);
            const bool el_passes_ss_iso   = (el_iso < el_iso_ss_cut); 
            const bool el_passes_eg_iso   = (el_iso < el_iso_eg_cut); 
            const bool el_passes_no_mhits = (mhit() <= el_mhits_cut); 
            const bool el_passes_d0       = (fabs(el_d0) < el_d0_cut); 
            const bool el_passes_id       = ((mediumId() & PassNoIso) == PassNoIso);
            const bool el_passes_hoe      = (hoe() < el_hoe_cut);
            const bool el_passes_3q       = (el_3q);

            // EGamma Medium WP (2012)
            // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification 
            // --------------------------------------------------------------------------- //

            // Isolation
            if (selection == Selection::EGammaMediumWPDenIso)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
                if (not el_passes_id)       {return false;}
            }

            // ID
            if (selection == Selection::EGammaMediumWPDenID)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
                if (not el_passes_eg_iso)   {return false;}
            }

            // Both ID and isolation 
            if (selection == Selection::EGammaMediumWPDenBoth)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
            }

            // Numerator
            if (selection == Selection::EGammaMediumWPNum)
            {
                if (el_is_crack)             {return false;}
                if (not el_passes_pt)        {return false;}
                if (not el_passes_trig_tag)  {return false;}
                if (not el_passes_id)        {return false;}
                if (not el_passes_eg_iso)    {return false;}
            }

            // SUS-13-013
            // --------------------------------------------------------------------------- //

            // Denominator
            if (selection == Selection::SameSignDenIso)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
                if (not el_passes_id)       {return false;}
                if (not el_passes_hoe)      {return false;}
                if (not el_passes_no_mhits) {return false;}
                if (not el_passes_3q)       {return false;}
                if (not el_passes_d0)       {return false;}
            }

            // ID
            if (selection == Selection::SameSignDenID)
            {
                if (el_is_crack)             {return false;}
                if (not el_passes_pt)        {return false;}
                if (not el_passes_trig_tag)  {return false;}
                if (not el_passes_ss_iso)    {return false;}
            }

            // Both ID and isolation 
            if (selection == Selection::SameSignDenBoth)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
            }

            // Numerator
            if (selection == Selection::SameSignNum)
            {
                if (el_is_crack)             {return false;}
                if (not el_passes_pt)        {return false;}
                if (not el_passes_trig_tag)  {return false;}
                if (not el_passes_id)        {return false;}
                if (not el_passes_hoe)       {return false;}
                if (not el_passes_no_mhits)  {return false;}
                if (not el_passes_3q)        {return false;}
                if (not el_passes_d0)        {return false;}
                if (not el_passes_ss_iso)    {return false;}
            }
        }

        // --------------------------------------------------------------------------- //
        // muons
        // --------------------------------------------------------------------------- //

        if (lepton_type == Lepton::Muon)
        {
            // cut values and variables
            const float mu_tag_pt      = tag().pt();
            const float mu_probe_pt    = probe().pt();
            const float mu_d0          = fabs(d0vtx()); 
            const float mu_iso         = (pfchiso03() + TMath::Max(0.0f, pfemiso03() + pfnhiso03() - 0.5f * dbeta03()))/mu_probe_pt;
            const float mu_iso_ss_cut  = 0.10;  // SUS-13-013 value
            const float mu_iso_pog_cut = 0.15;  // Muon POG value (not sure about this one)
            const float mu_d0_ss_cut   = 0.005; // SUS-13-013 value
            const float mu_tag_pt_cut  = 30.0;

            // cut decisions 
            const bool mu_passes_pt       = (mu_tag_pt > mu_tag_pt_cut);
            const bool mu_passes_trig_tag = (is_data ? HLT_IsoMu24_eta2p1_tag() != 0 : true);
            const bool mu_passes_pog_iso  = (mu_iso < mu_iso_pog_cut); 
            const bool mu_passes_pog_id   = ((leptonSelection() & LeptonSelection::PassMuIsHPASS) == LeptonSelection::PassMuIsHPASS);
            const bool mu_passes_ss_iso   = (mu_iso < mu_iso_ss_cut); 
            const bool mu_passes_ss_d0    = (fabs(mu_d0) < mu_d0_ss_cut); 
            const bool mu_passes_ss_id    = ((leptonSelection() & LeptonSelection::PassMuIsHPASS) == LeptonSelection::PassMuIsHPASS) && (mu_passes_ss_d0);

            // Muon POG Selections (2012)
            // From: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
            // --------------------------------------------------------------------------- //

            // Isolation
            if (selection == Selection::MuTightWPDenIso)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
                if (not mu_passes_pog_id)   {return false;}
            }

            // ID
            if (selection == Selection::MuTightWPDenID)
            {
                if (not mu_passes_pt)        {return false;}
                if (not mu_passes_trig_tag)  {return false;}
                if (not mu_passes_pog_iso)   {return false;}
            }

            // Both ID and isolation
            if (selection == Selection::MuTightWPDenBoth)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
            }

            // Numerator
            if (selection == Selection::MuTightWPNum)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
                if (not mu_passes_pog_id)   {return false;}
                if (not mu_passes_pog_iso)  {return false;}
            }

            // SUS-13-013 Selections
            // --------------------------------------------------------------------------- //

            // Isolation
            if (selection == Selection::SameSignDenIso)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
                if (not mu_passes_ss_id)    {return false;}
            }

            // ID
            if (selection == Selection::SameSignDenID)
            {
                if (not mu_passes_pt)        {return false;}
                if (not mu_passes_trig_tag)  {return false;}
                if (not mu_passes_ss_iso)    {return false;}
            }

            // Both ID and isolation
            if (selection == Selection::SameSignDenBoth)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
            }

            // Numerator
            if (selection == Selection::SameSignNum)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
                if (not mu_passes_ss_id)    {return false;}
                if (not mu_passes_ss_iso)   {return false;}
            }
        }

        // other values are invalid
        if (lepton_type == Lepton::static_size)
        {
            std::cout << "lepton_type is invalid" << std::endl;
            return false;
        }

        // if we got here -- it's selected
        return true;
    }

} // namespace tnp
 
