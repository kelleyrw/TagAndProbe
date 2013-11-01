// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

// ROOT includes
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// Tools
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "TagAndProbe/Analysis/interface/Dataset.h"

using std::cout;
using std::endl;

// ------------------------------------------------------------------------------------ //
// helper functions/classes 
// ------------------------------------------------------------------------------------ //

struct contains
{
    contains(const std::string& sub_str) : m_sub_str(sub_str) {}
    bool operator() (const std::string& str) const {return lt::string_contains(str, m_sub_str);}
    std::string m_sub_str;
};
        
// ------------------------------------------------------------------------------------ //
// create plots 
// ------------------------------------------------------------------------------------ //

TCanvas* CreateEfficienyPlot1D
(
    TH1* const h_ds1,
    TH1* const h_ds2,
    TH1* const h_sf,
    const std::string& title,
    const std::string& ds1_label,
    const std::string& ds2_label,
    const std::string& option = ""
)
{
    const std::string canvas_title = lt::string_replace_all(h_ds1->GetName(), "h_", "p_");
    TCanvas* c = new TCanvas(canvas_title.c_str(), canvas_title.c_str(), 800, 600);
    c->SetGrid();

    h_ds1->SetTitle(title.c_str());
    h_ds1->GetYaxis()->SetRangeUser(0.0, 1.1);
    h_ds1->SetLineColor(kBlue);
    h_ds1->SetLineWidth(2);
    h_ds1->SetMarkerStyle(20);
    h_ds1->SetMarkerSize(1.0);
    h_ds1->SetMarkerColor(kBlue);
    h_ds1->SetStats(false);
    h_ds1->Draw(option.c_str());

    h_ds2->SetLineColor(kRed);
    h_ds2->SetLineWidth(2);
    h_ds2->SetMarkerStyle(20);
    h_ds2->SetMarkerSize(1.0);
    h_ds2->SetMarkerColor(kRed);
    h_ds2->SetStats(false);
    h_ds2->Draw(Form("%s same", option.c_str()));

    h_sf->SetLineColor(kGreen+2);
    h_sf->SetLineWidth(2);
    h_sf->SetMarkerStyle(20);
    h_sf->SetMarkerSize(1.0);
    h_sf->SetMarkerColor(kGreen+2);
    h_sf->SetStats(false);
    h_sf->Draw(Form("%s same", option.c_str()));

    TLegend* legend = new TLegend(0.2, 0.5, 0.8, 0.15);
    legend->AddEntry(h_ds1, ds1_label.c_str(), "lep");
    legend->AddEntry(h_ds2, ds2_label.c_str(), "lep");
    legend->AddEntry(h_sf , Form("%s/%s", ds1_label.c_str(), ds2_label.c_str()), "lep");
    legend->SetFillColor(0);  // 0 makes it the background clear on the pad
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.030);
    legend->Draw();


    return c;
}

TCanvas* CreateScaleFactorPlot2D
(
    TH1* const h_sf,
    const std::string& title,
    const std::string& ds1_label,
    const std::string& ds2_label,
    const std::string& option = ""
)
{
    const std::string canvas_title = lt::string_replace_all(h_sf->GetName(), "h_", "p_");
    TCanvas* c = new TCanvas(canvas_title.c_str(), canvas_title.c_str(), 800, 600);
    c->SetGrid();
    if (lt::string_contains(option, "colz"))
    {
        c->SetRightMargin(0.15);
    }

    h_sf->GetZaxis()->SetRangeUser(0.7, 1.1);
    h_sf->SetTitle(title.c_str());
    h_sf->SetTitle(title.c_str());
    h_sf->SetLineColor(kBlack);
    h_sf->SetLineWidth(2);
    h_sf->SetMarkerStyle(20);
    h_sf->SetMarkerSize(1.0);
    h_sf->SetMarkerColor(kBlack);
    h_sf->SetStats(false);
    h_sf->Draw(option.c_str());

    return c;
}

// ------------------------------------------------------------------------------------ //
// The main program 
// ------------------------------------------------------------------------------------ //

int main(int argc, char **argv)
try
{
    // FWLite libs
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();

    // global style
    gStyle->SetPaintTextFormat("1.3");

    // -------------------------------------------------------------------------------- //
    // parse the inputs
    // -------------------------------------------------------------------------------- //

    // check that the python is passed
    if (argc < 2)
    {
        throw std::invalid_argument(Form("Usage : %s [parameters.py]", argv[0]));
    }

    // check that pset contains "process" 
    const std::string pset_filename = argv[1];
    if (!edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process"))
    {
        throw std::invalid_argument(Form("[tnp_eff_plots] Error: ParametersSet 'process' is missing in your configuration file"));
    }

    // get the python configuration
    const edm::ParameterSet& process               = edm::readPSetsFrom(pset_filename)->getParameter<edm::ParameterSet>("process");
    const std::vector<edm::ParameterSet>& tnp_cfgs = process.getParameter<std::vector<edm::ParameterSet> >("tnp_compare");

    // -------------------------------------------------------------------------------- //
    // for each pset, do the comparison 
    // -------------------------------------------------------------------------------- //

    for (std::vector<edm::ParameterSet>::const_iterator tnp_cfg_iter = tnp_cfgs.begin(); tnp_cfg_iter != tnp_cfgs.end(); tnp_cfg_iter++)
    {
        // convenience
        const edm::ParameterSet& tnp_cfg = *tnp_cfg_iter;

        // parameters 
        const bool verbose               = tnp_cfg.getParameter<bool>("verbose");
        const std::string suffix         = tnp_cfg.getParameter<std::string>("suffix");

        // datasets
        const tnp::Dataset dataset1 = tnp::Dataset(tnp_cfg.getParameter<edm::ParameterSet>("dataset1"));
        const tnp::Dataset dataset2 = tnp::Dataset(tnp_cfg.getParameter<edm::ParameterSet>("dataset2"));

        // path to results 
        const std::string eff_results_path = tnp_cfg.getParameter<std::string>("eff_results_path");
        const std::string output_label     = tnp_cfg.getParameter<std::string>("output_label");

        if (verbose)
        {
            cout << "[tnp_compare] eff_results_path = " << eff_results_path << endl;
            cout << "[tnp_compare] dataset1         = " << dataset1.m_name << endl;
            cout << "[tnp_compare] dataset2         = " << dataset2.m_name << endl;
        }

        // -------------------------------------------------------------------------------- //
        // get the histograms for the efficiency 
        // -------------------------------------------------------------------------------- //

        rt::TH1Container hc1(Form("%s/%s_eff.root", eff_results_path.c_str(), dataset1.m_name.c_str()));
        rt::TH1Container hc2(Form("%s/%s_eff.root", eff_results_path.c_str(), dataset2.m_name.c_str()));

        // check that the containers have the same histogram names
        const std::vector<std::string> hc1_names = hc1.GetListOfHistograms();
        const std::vector<std::string> hc2_names = hc2.GetListOfHistograms();
        if (not std::equal(hc1_names.begin(), hc1_names.end(), hc2_names.begin()))
        {
            throw std::runtime_error("[tnp_compare] The list of histograms from the two inputs datasets do not match");
        }

        // efficiency plot names
        const std::vector<std::string> eff_names = lt::filter_container(hc1_names, contains("eff"));

        // print histogram in the file
        if (verbose)
        {
            cout << "[tnp_compare] the list of histograms in efficiency file is:" << endl;
            hc1.List();
        }

        // -------------------------------------------------------------------------------- //
        // scale factors (dataset1 / dataset2)
        // -------------------------------------------------------------------------------- //

        rt::TH1Container hc_sf;
        for (size_t i = 0; i != eff_names.size(); i++)
        {
            const std::string& eff_hist_name = eff_names.at(i);    
            const std::string& sf_hist_name  = lt::string_replace_all(eff_hist_name, "eff", "sf");
            const std::string& sf_hist_title = lt::string_replace_all(hc1[eff_hist_name]->GetTitle(), "Efficiency", "Scale Factor");
            hc_sf.Add(rt::DivideHists(hc1[eff_hist_name], hc2[eff_hist_name], sf_hist_name, sf_hist_title)); 
        }

        // -------------------------------------------------------------------------------- //
        // tables and plots of efficiency 
        // hard coded for now
        // -------------------------------------------------------------------------------- //

        // using function rt::CreateTableFromHist
        // https://github.com/kelleyrw/AnalysisTools/blob/master/RootTools/interface/TH1Tools.h#L236

        std::vector<CTable> dataset1_eff_tables;
        std::vector<CTable> dataset2_eff_tables;
        std::vector<CTable> scale_factor_tables;
        std::vector<TCanvas*> compare_plots;

        // 1D: eff(pt)
        if (std::find(eff_names.begin(), eff_names.end(), std::string("h_eff_pt")) != eff_names.end())
        {
            dataset1_eff_tables.push_back(rt::CreateTableFromHist(hc1  ["h_eff_pt"], Form("%s Efficiency vs $p_{T}$", dataset1.m_title.c_str()), "$p_{T}$", "", "GeV", "", "1.2", "1.0", "1.0", false));
            dataset2_eff_tables.push_back(rt::CreateTableFromHist(hc2  ["h_eff_pt"], Form("%s Efficiency vs $p_{T}$", dataset2.m_title.c_str()), "$p_{T}$", "", "GeV", "", "1.2", "1.0", "1.0", false));
            scale_factor_tables.push_back(rt::CreateTableFromHist(hc_sf["h_sf_pt" ], "Scale Factor vs $p_{T}$", "$p_{T}$", "", "GeV", "", "1.2", "1.0", "1.0", false));
            compare_plots.push_back(CreateEfficienyPlot1D(hc1["h_eff_pt"], hc2["h_eff_pt"], hc_sf["h_sf_pt"], "Efficiency vs p_{T};p_{T} (GeV); Efficiency", dataset1.m_title, dataset2.m_title));
        }

        // 1D: eff(eta)
        if (std::find(eff_names.begin(), eff_names.end(), std::string("h_eff_eta")) != eff_names.end())
        {
            dataset1_eff_tables.push_back(rt::CreateTableFromHist(hc1  ["h_eff_eta"], Form("%s Efficiency vs $\\eta$", dataset1.m_title.c_str()), "$\\eta$", "", "", "", "1.2", "1.2", "1.0", false));
            dataset2_eff_tables.push_back(rt::CreateTableFromHist(hc2  ["h_eff_eta"], Form("%s Efficiency vs $\\eta$", dataset2.m_title.c_str()), "$\\eta$", "", "", "", "1.2", "1.2", "1.0", false));
            scale_factor_tables.push_back(rt::CreateTableFromHist(hc_sf["h_sf_eta" ], "Scale Factor vs $\\eta$", "$\\eta$", "", "", "", "1.2", "1.2", "1.0", false));
            compare_plots.push_back(CreateEfficienyPlot1D(hc1["h_eff_eta"], hc2["h_eff_eta"], hc_sf["h_sf_eta"], "Efficiency vs #eta;#eta; Efficiency", dataset1.m_title, dataset2.m_title));
        }

        // 1D: eff(phi)
        if (std::find(eff_names.begin(), eff_names.end(), std::string("h_eff_phi")) != eff_names.end())
        {
            dataset1_eff_tables.push_back(rt::CreateTableFromHist(hc1  ["h_eff_phi"], Form("%s Efficiency vs $\\phi$", dataset1.m_title.c_str()), "$\\phi$", "", "", "", "1.2", "1.2", "1.0", false));
            dataset2_eff_tables.push_back(rt::CreateTableFromHist(hc2  ["h_eff_phi"], Form("%s Efficiency vs $\\phi$", dataset2.m_title.c_str()), "$\\phi$", "", "", "", "1.2", "1.2", "1.0", false));
            scale_factor_tables.push_back(rt::CreateTableFromHist(hc_sf["h_sf_phi" ], "Scale Factor vs $\\phi$", "$\\phi$", "", "", "", "1.2", "1.2", "1.0", false));
            compare_plots.push_back(CreateEfficienyPlot1D(hc1["h_eff_phi"], hc2["h_eff_phi"], hc_sf["h_sf_phi"], "Efficiency vs #phi;#phi; Efficiency", dataset1.m_title, dataset2.m_title));
        }

        // 1D: eff(nvtx)
        if (std::find(eff_names.begin(), eff_names.end(), std::string("h_eff_nvtx")) != eff_names.end())
        {
            dataset1_eff_tables.push_back(rt::CreateTableFromHist(hc1  ["h_eff_nvtx"], Form("%s Efficiency vs \\# vertices", dataset1.m_title.c_str()), "\\# vertices", "", "", "", "1.2", "1.0", "1.0", false));
            dataset2_eff_tables.push_back(rt::CreateTableFromHist(hc2  ["h_eff_nvtx"], Form("%s Efficiency vs \\# vertices", dataset2.m_title.c_str()), "\\# vertices", "", "", "", "1.2", "1.0", "1.0", false));
            scale_factor_tables.push_back(rt::CreateTableFromHist(hc_sf["h_sf_nvtx" ], "Scale Factor vs \\# vertices", "\\# vertices", "", "", "", "1.2", "1.0", "1.0", false));
            compare_plots.push_back(CreateEfficienyPlot1D(hc1["h_eff_nvtx"], hc2["h_eff_nvtx"], hc_sf["h_sf_nvtx"], "Efficiency vs # vertices;# vertices; Efficiency", dataset1.m_title, dataset2.m_title));
        }

        // 2D: eff(pt, eta)
        if (std::find(eff_names.begin(), eff_names.end(), std::string("h_eff_pt_vs_eta")) != eff_names.end())
        {
            dataset1_eff_tables.push_back(rt::CreateTableFromHist(hc1  ["h_eff_pt_vs_eta"], Form("%s Efficiency vs $p_{T}$ and $\\eta$", dataset1.m_title.c_str()), "$p_{T}$", "$\\eta$", "", "", "1.2", "1.0", "1.2", true));
            dataset2_eff_tables.push_back(rt::CreateTableFromHist(hc2  ["h_eff_pt_vs_eta"], Form("%s Efficiency vs $p_{T}$ and $\\eta$", dataset2.m_title.c_str()), "$p_{T}$", "$\\eta$", "", "", "1.2", "1.0", "1.2", true));
            scale_factor_tables.push_back(rt::CreateTableFromHist(hc_sf["h_sf_pt_vs_eta" ], "Scale Factor vs $p_{T}$ and $\\eta$", "$p_{T}$", "$\\eta$", "", "", "1.2", "1.0", "1.2", true));
            compare_plots.push_back(CreateScaleFactorPlot2D(hc_sf["h_sf_pt_vs_eta"], "Scale Factor vs p_{T} and #eta;#eta;p_{T} (GeV);Efficiency", dataset1.m_title, dataset2.m_title, "colz"));
        }

        // 2D: eff(eta, phi)
        if (std::find(eff_names.begin(), eff_names.end(), std::string("h_eff_eta_vs_phi")) != eff_names.end())
        {
            dataset1_eff_tables.push_back(rt::CreateTableFromHist(hc1  ["h_eff_eta_vs_phi"], Form("%s Efficiency vs $\\eta$ and $\\phi$", dataset1.m_title.c_str()), "$\\eta$", "$\\phi$", "", "", "1.2", "1.2", "1.2", true));
            dataset2_eff_tables.push_back(rt::CreateTableFromHist(hc2  ["h_eff_eta_vs_phi"], Form("%s Efficiency vs $\\eta$ and $\\phi$", dataset2.m_title.c_str()), "$\\eta$", "$\\phi$", "", "", "1.2", "1.2", "1.2", true));
            scale_factor_tables.push_back(rt::CreateTableFromHist(hc_sf["h_sf_eta_vs_phi" ], "Scale Factor vs $\\eta$ and $\\phi$", "$\\eta$", "$\\phi$", "", "", "1.2", "1.2", "1.2", true));
            compare_plots.push_back(CreateScaleFactorPlot2D(hc_sf["h_sf_eta_vs_phi"], "Scale Factor vs #eta and #phi;#phi;#eta;Efficiency", dataset1.m_title, dataset2.m_title, "colz"));
        }

        // -------------------------------------------------------------------------------- //
        // output 
        // -------------------------------------------------------------------------------- //

        // output path
        const std::string output_plot_path = Form("%s/%s", eff_results_path.c_str(), output_label.c_str());
        lt::mkdir(output_plot_path, /*force=*/true);
        rt::CopyIndexPhp(output_plot_path);

        // print the tables
        std::ofstream txt_out(Form("%s/eff_tables.txt", output_plot_path.c_str()), std::ofstream::out);
        std::ofstream tex_out(Form("%s/eff_tables.tex", output_plot_path.c_str()), std::ofstream::out);
        tex_out << "\\documentclass{article}" << endl
                << "\\begin{document}"        << endl;
        for (size_t i = 0; i != dataset1_eff_tables.size(); i++)
        {
            CTable& t_ds1 = dataset1_eff_tables.at(i);
            CTable& t_ds2 = dataset2_eff_tables.at(i);
            CTable& t_sf  = scale_factor_tables.at(i);

            // txt file
            txt_out << t_ds1 << "\n" << t_ds2 << "\n" << t_sf << endl;

            // tex file
            tex_out << t_ds1.getTex() << endl 
                    << t_ds2.getTex() << endl
                    << t_sf.getTex()  << endl;
        }
        tex_out <<"\\end{document}" << endl;

        // print the plots
        if (not suffix.empty())
        {
            for (size_t i = 0; i != compare_plots.size(); i++)
            {
                compare_plots.at(i)->Print(Form("%s/%s.%s", output_plot_path.c_str(), compare_plots.at(i)->GetName(), suffix.c_str())); 
            }
        }

        // cleanup
        lt::delete_container(compare_plots);

    } // end tnp_cfgs loop

    // done
    return 0;
}
catch (std::exception& e)
{
    std::cerr << "[tnp_compare] Error: failed..." << std::endl;
    std::cerr << e.what() << std::endl;
    return 1;
}
