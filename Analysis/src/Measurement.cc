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
        if (lt::string_lower(sel_name) == "egammadenid"      ) {return Selection::EGammaDenID;      } 
        if (lt::string_lower(sel_name) == "egammadeniso"     ) {return Selection::EGammaDenIso;     } 
        if (lt::string_lower(sel_name) == "egammadenboth"    ) {return Selection::EGammaDenBoth;    } 
        if (lt::string_lower(sel_name) == "egammanum"        ) {return Selection::EGammaNum;        } 
        if (lt::string_lower(sel_name) == "samesigneldenid"  ) {return Selection::SameSignElDenID;  } 
        if (lt::string_lower(sel_name) == "samesigneldeniso" ) {return Selection::SameSignElDenIso; } 
        if (lt::string_lower(sel_name) == "samesigneldenboth") {return Selection::SameSignElDenBoth;} 
        if (lt::string_lower(sel_name) == "samesignelnum"    ) {return Selection::SameSignElNum;    } 
        if (lt::string_lower(sel_name) == "samesignmudenid"  ) {return Selection::SameSignMuDenID;  } 
        if (lt::string_lower(sel_name) == "samesignmudeniso" ) {return Selection::SameSignMuDenIso; } 
        if (lt::string_lower(sel_name) == "samesignmudenboth") {return Selection::SameSignMuDenBoth;} 
        if (lt::string_lower(sel_name) == "samesignmunum"    ) {return Selection::SameSignMuNum;    } 
        throw std::invalid_argument("[tnp::GetSelectionFromString]: ERROR - invalid value!"); 
    }

    std::string GetStringFromSelection(const Selection::value_type sel_type)
    {
        if (sel_type == Selection::EGammaDenID       ) return "egammadenid";
        if (sel_type == Selection::EGammaDenIso      ) return "egammadeniso";
        if (sel_type == Selection::EGammaDenBoth     ) return "egammadenboth";
        if (sel_type == Selection::EGammaNum         ) return "egammanum";
        if (sel_type == Selection::SameSignElDenID   ) return "samesigneldenid";
        if (sel_type == Selection::SameSignElDenIso  ) return "samesigneldeniso";
        if (sel_type == Selection::SameSignElDenBoth ) return "samesigneldenboth";
        if (sel_type == Selection::SameSignElNum     ) return "samesignelnum";
        if (sel_type == Selection::SameSignMuDenID   ) return "samesignmudenid";
        if (sel_type == Selection::SameSignMuDenIso  ) return "samesignmudeniso";
        if (sel_type == Selection::SameSignMuDenBoth ) return "samesignmudenboth";
        if (sel_type == Selection::SameSignMuNum     ) return "samesignmunum";
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
            const float el_iso_eg_cut  = (el_is_endcap ? (el_probe_pt < 20.0 ? 0.10 : 0.15) : 0.15);  // egamma value 
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
            //const bool el_is_rand         = true;//(not ((qTag() < 0 && rnd() < 0.5) || (qTag() > 0 && rnd() >= 0.5)));

            // Denominators
            // --------------------------------------------------------------------------- //

            // Isolation

            // Same Sign Denominator Isolation selections
            if (selection == Selection::SameSignElDenIso)
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

            // EGamma Denominator Isolation selections
            if (selection == Selection::EGammaDenIso)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
                if (not el_passes_id)       {return false;}
            }

            // ID

            // Same Sign Denominator ID selections
            if (selection == Selection::SameSignElDenID)
            {
                if (el_is_crack)             {return false;}
                if (not el_passes_pt)        {return false;}
                if (not el_passes_trig_tag)  {return false;}
                if (not el_passes_ss_iso)    {return false;}
            }

            // EGamma Denominator ID selections
            if (selection == Selection::EGammaDenID)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
                if (not el_passes_eg_iso)   {return false;}
            }

            // Both

            // Same Sign Denominator Both selections
            if (selection == Selection::SameSignElDenBoth)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
            }

            // EGamma Denominator Both selections
            if (selection == Selection::EGammaDenBoth)
            {
                if (el_is_crack)            {return false;}
                if (not el_passes_pt)       {return false;}
                if (not el_passes_trig_tag) {return false;}
            }

            // Numerators
            // --------------------------------------------------------------------------- //

            // Same Sign Numerator Isolation selections
            if (selection == Selection::SameSignElNum)
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

            // EGamma Numerator Isolation selections
            if (selection == Selection::EGammaNum)
            {
                if (el_is_crack)             {return false;}
                if (not el_passes_pt)        {return false;}
                if (not el_passes_trig_tag)  {return false;}
                if (not el_passes_id)        {return false;}
                if (not el_passes_eg_iso)    {return false;}
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
            const float mu_iso_ss_cut  = 0.10;                                                        // ss2012 value
            const float mu_d0_cut      = 0.005;
            const float mu_tag_pt_cut  = 30.0;

            // cut decisions 
            const bool mu_passes_pt       = (mu_tag_pt > mu_tag_pt_cut);
            const bool mu_passes_trig_tag = (is_data ? HLT_IsoMu24_eta2p1_tag() != 0 : true);
            const bool mu_passes_ss_iso   = (mu_iso < mu_iso_ss_cut); 
            const bool mu_passes_d0       = (fabs(mu_d0) < mu_d0_cut); 
            const bool mu_passes_id       = ((leptonSelection() & LeptonSelection::PassMuIsHPASS) == LeptonSelection::PassMuIsHPASS) && (mu_passes_d0);

            // Denominators
            // --------------------------------------------------------------------------- //

            // Isolation

            // Same Sign Denominator Isolation selections
            if (selection == Selection::SameSignMuDenIso)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
                if (not mu_passes_id)       {return false;}
            }

            // ID

            // Same Sign Denominator ID selections
            if (selection == Selection::SameSignMuDenID)
            {
                if (not mu_passes_pt)        {return false;}
                if (not mu_passes_trig_tag)  {return false;}
                if (not mu_passes_ss_iso)    {return false;}
            }

            // Both

            // Same Sign Denominator Both selections
            if (selection == Selection::SameSignMuDenBoth)
            {
                if (not mu_passes_pt)       {return false;}
                if (not mu_passes_trig_tag) {return false;}
            }

            // Numerators
            // --------------------------------------------------------------------------- //

            // Same Sign Numerator Isolation selections
            if (selection == Selection::SameSignMuNum)
            {
                if (not mu_passes_pt)        {return false;}
                if (not mu_passes_trig_tag)  {return false;}
                if (not mu_passes_id)        {return false;}
                if (not mu_passes_ss_iso)    {return false;}
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
 
