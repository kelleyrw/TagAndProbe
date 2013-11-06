//#include "TagAndProbe/Analysis/interface/PerformFits.h"
//#include "TagAndProbe/Analysis/interface/Dataset.h"
//#include "TagAndProbe/Analysis/interface/LeptonSelections.h"
//#include "TagAndProbe/Analysis/interface/LeptonTree.h"
//#include "TagAndProbe/Analysis/interface/Measurement.h"
//#include "AnalysisTools/RootTools/interface/RootTools.h"
//#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
//#include "TSystem.h"
//#include <algorithm>
//
//void TestSingleFit()
{
    if (TString(gSystem->Getenv("SCRAM_ARCH")).Contains("osx"))
    {
        gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libTagAndProbeAnalysis.dylib");
    }
    else
    {
        gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libTagAndProbeAnalysis.so");
    }

    const float pt_bins[]  = {10, 15, 20, 30, 40, 50, 200};
    const size_t num_pt_bins = 6;

    const float eta_bins[] = {0, 0.8, 1.4442, 1.566, 2.0, 2.5};
    const size_t num_eta_bins = 5;

    size_t pt_bin  = 0;
    size_t eta_bin = 3;
/*     for (size_t pt_bin = 0; pt_bin != num_pt_bins; pt_bin++) */
    {
/*         for (size_t eta_bin = 0; eta_bin != num_eta_bins; eta_bin++) */
        {
/*             if (!((eta_bin==3 || eta_bin==2) && pt_bin == 4)) {continue;} */
            
            rt::TH1Container hc("plots/SameSign/electron/SameSignDenBoth_SameSignNum/data_single_el.root");
            TH1* h_pass = hc[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
            TH1* h_fail = hc[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];

            rt::TH1Container hc_mc("plots/SameSign/electron/SameSignDenBoth_SameSignNum/dy_full.root");
            TH1* h_pass_template = hc_mc[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
            TH1* h_fail_template = hc_mc[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];

            tnp::Result result = tnp::PerformSimultaneousFit
            (
                 tnp::Model::MCTemplate,
                 tnp::Model::MCTemplate,
                 tnp::Model::Exponential,
                 tnp::Model::ErfExp,
                 h_pass,
                 h_fail, 
                 60.0, 
                 120.0, 
                 2.0, 
                 /*pt_bin_label  = */Form("%1.2f < |#eta| < %1.2f", eta_bins[eta_bin], eta_bins[eta_bin+1]), 
                 /*eta_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", pt_bins[pt_bin], pt_bins[pt_bin+1]), 
                 h_pass_template, 
                 h_fail_template
            );

            result.cpass->Print(Form("plots/test/p_pass_pt%lu_vs_eta%lu.png", pt_bin, eta_bin));
            result.cfail->Print(Form("plots/test/p_fail_pt%lu_vs_eta%lu.png", pt_bin, eta_bin));
        }
    }
}
