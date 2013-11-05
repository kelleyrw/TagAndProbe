// C++ includes
#include <iostream>
#include <string>
#include <stdexcept>

// BOOST
#include "boost/multi_array.hpp"

// ROOT includes
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// Tools
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "TagAndProbe/Analysis/interface/PerformFits.h"
#include "TagAndProbe/Analysis/interface/Dataset.h"
#include "TagAndProbe/Analysis/interface/LeptonSelections.h"
#include "TagAndProbe/Analysis/interface/LeptonTree.h"
#include "TagAndProbe/Analysis/interface/Measurement.h"

using namespace std;

// Array of Models from the input stings 
// ------------------------------------------------------------------------------------ //

// using boost::multiarray
// http://www.boost.org/doc/libs/1_49_0/libs/multi_array/doc/user.html
typedef boost::multi_array<tnp::Model::value_type, 2> ModelArray2D;
typedef boost::multi_array<tnp::Model::value_type, 3> ModelArray3D;

// 4 models per bin: signal pass, signal fail, background pass, background fail
const unsigned int num_categories = 4;

// parse the input vstring and return an array of models (1D bins)
ModelArray2D GetModelArrayFromVString
(
    const std::vector<std::string>& model_strings, 
    const std::vector<double> bins
)
{
    // number of bins: bins array size - 1 (e.g. {1,2,3} has 2 bins)
    const unsigned int num_bins = bins.size() - 1; 
    
    // check that the # elements lines up
    if (num_bins*num_categories != model_strings.size())
    {
        throw std::invalid_argument(Form("[tnp_eff_plots] Error: # of bins (%u) do not add up with the # of models (%lu) in the configuration", num_bins*num_categories, model_strings.size()));
    }

    // loop and fill array
    ModelArray2D result(boost::extents[num_bins][num_categories]);
    for (size_t bin = 0; bin != num_bins; bin++)
    {
        for (size_t model_bin = 0; model_bin != num_categories; model_bin++)
        {

            size_t index = bin*(num_categories) + model_bin;
            result[bin][model_bin] = tnp::GetModelFromString(model_strings.at(index));
        }
    }
   
    // done
    return result;
}

// parse the input vstring and return an array of models (2D bins)

ModelArray3D GetModelArrayFromVString
(
    const std::vector<std::string>& model_strings, 
    const std::vector<double> a_bins,
    const std::vector<double> b_bins
)
{
    // number of bins: bins array size - 1 (e.g. {1,2,3} has 2 bins)
    const unsigned int num_a_bins = a_bins.size() - 1; 
    const unsigned int num_b_bins = b_bins.size() - 1; 
    
    // check that the # elements lines up
    if ((num_a_bins)*(num_b_bins)*num_categories != model_strings.size())
    {
        throw std::invalid_argument(Form("[tnp_eff_plots] Error: # of bins (%u) do not add up with the # of models (%lu) in the configuration", (num_a_bins)*(num_b_bins)*num_categories, model_strings.size()));
    }

    // loop and fill array
    ModelArray3D result(boost::extents[num_a_bins][num_b_bins][num_categories]);
    for (size_t a_bin = 0; a_bin != num_a_bins; a_bin++)
    {
        for (size_t b_bin = 0; b_bin != num_b_bins; b_bin++)
        {
            for (size_t model_bin = 0; model_bin != num_categories; model_bin++)
            {

                size_t index = a_bin*(num_b_bins * num_categories) + b_bin*(num_categories) + model_bin;
                result[a_bin][b_bin][model_bin] = tnp::GetModelFromString(model_strings.at(index));
            }
        }
    }
   
    // done
    return result;
}

// Print the fitted histogram 
// ------------------------------------------------------------------------------------ //

void PrintCanvas(TCanvas* const canvas, const std::string& filename, const std::string& suffix = "png")
{
    const std::string dirname  = lt::dirname(filename);
    const std::string basename = lt::basename(filename);
    if (suffix == "all")
    {
        lt::mkdir(dirname + "/png", /*force=*/true);
        lt::mkdir(dirname + "/eps", /*force=*/true);
        lt::mkdir(dirname + "/pdf", /*force=*/true);
        rt::CopyIndexPhp(dirname + "/png");
        rt::CopyIndexPhp(dirname + "/eps");
        rt::CopyIndexPhp(dirname + "/pdf");
        canvas->Print(Form("%s/eps/%s.eps", dirname.c_str(), basename.c_str()));
        canvas->Print(Form("%s/png/%s.png", dirname.c_str(), basename.c_str()));
        canvas->Print(Form("%s/pdf/%s.pdf", dirname.c_str(), basename.c_str()));
    }
    else
    {
        lt::mkdir(dirname + "/" + suffix, /*force=*/true);
        rt::CopyIndexPhp(dirname + "/" + suffix);
        canvas->Print(Form("%s/%s/%s.%s", dirname.c_str(), suffix.c_str(), basename.c_str(), suffix.c_str()));
    }
}

// The main program 
// ------------------------------------------------------------------------------------ //

int main(int argc, char **argv)
try
{
    // FWLite libs
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();

    // parse the inputs
    // -------------------------------------------------------------------------------------------------//

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
    const edm::ParameterSet& process = edm::readPSetsFrom(pset_filename)->getParameter<edm::ParameterSet>("process");
    const edm::ParameterSet& tnp_cfg = process.getParameter<edm::ParameterSet>("tnp_eff_plots");

    // input values 
    const tnp::Lepton::value_type lepton_type = tnp::GetLeptonFromString(tnp_cfg.getParameter<std::string>("lepton_type"));
    const double mass_low                     = tnp_cfg.getParameter<double>("mass_low" );
    const double mass_high                    = tnp_cfg.getParameter<double>("mass_high");
    const double mass_bin_width               = tnp_cfg.getParameter<double>("mass_bin_width");
    const bool verbose                        = tnp_cfg.getParameter<bool>("verbose");
    const std::string suffix                  = tnp_cfg.getParameter<std::string>("suffix");
    const std::string analysis_path           = lt::getenv("CMSSW_BASE") + "/src/TagAndProbe/Analysis";
    const std::string output_label            = tnp_cfg.getParameter<string>("output_label");
    const std::vector<tnp::Dataset> datasets  = tnp::GetDatasetsFromVPSet(tnp_cfg.getParameter<std::vector<edm::ParameterSet> >("datasets"));

    // numerator and denominator    
    const tnp::Selection::value_type num_selection = tnp::GetSelectionFromString(tnp_cfg.getParameter<std::string>("numerator"  ));
    const tnp::Selection::value_type den_selection = tnp::GetSelectionFromString(tnp_cfg.getParameter<std::string>("denominator"));

    // bins
    const std::vector<double> pt_bins   = tnp_cfg.getParameter<std::vector<double> >("pt_bins"  );
    const std::vector<double> eta_bins  = tnp_cfg.getParameter<std::vector<double> >("eta_bins" );
    const std::vector<double> phi_bins  = tnp_cfg.getParameter<std::vector<double> >("phi_bins" );
    const std::vector<double> nvtx_bins = tnp_cfg.getParameter<std::vector<double> >("nvtx_bins");
    const unsigned int num_pt_bins      = pt_bins.size()-1;
    const unsigned int num_eta_bins     = eta_bins.size()-1;
    const unsigned int num_phi_bins     = phi_bins.size()-1;
    const unsigned int num_nvtx_bins    = nvtx_bins.size()-1;

    // use |eta| ?
    const std::string phi_title = (not phi_bins.empty() and phi_bins.front() >= 0 ? "|#phi|" : "#phi");
    const std::string eta_title = (not eta_bins.empty() and eta_bins.front() >= 0 ? "|#eta|" : "#eta");

    // mc template histograms
    rt::TH1Container hc_template(tnp_cfg.getParameter<std::string>("mc_template_file"));

    // for each dataset makes the set of histograms
    // -------------------------------------------------------------------------------------------------//

    for (std::vector<tnp::Dataset>::const_iterator dataset_iter = datasets.begin(); dataset_iter != datasets.end(); dataset_iter++)
    {
        // convenience
        const tnp::Dataset& dataset = *dataset_iter;

        // mass plot ROOT file name
        // (i.e. analysis_path/plots/output_label/lepton_type/den_num/dataset.root)
        const std::string mass_plot_file_name = Form("%s/plots/%s/%s/%s_%s/%s.root",
            analysis_path.c_str(),
            output_label.c_str(),
            GetStringFromLepton(lepton_type).c_str(),
            GetStringFromSelection(den_selection).c_str(),
            GetStringFromSelection(num_selection).c_str(),
            dataset.m_name.c_str()
        );

        // output efficiency ROOT file name
        // (i.e. analysis_path/plots/output_label/lepton_type/den_num/dataset_eff.root)
        const std::string eff_plot_file_name = Form("%s/plots/%s/%s/%s_%s/%s_eff.root",
            analysis_path.c_str(),
            output_label.c_str(),
            GetStringFromLepton(lepton_type).c_str(),
            GetStringFromSelection(den_selection).c_str(),
            GetStringFromSelection(num_selection).c_str(),
            dataset.m_name.c_str()
        );

        // book efficiency hists 
        rt::TH1Container hc_mass(mass_plot_file_name);
        if (verbose)
        {
            cout << "[tnp_eff_plots] Fitting the following plots:" << endl;
            hc_mass.List();
        }

        // book efficiency hists 
        rt::TH1Container hc;

        // pt bins
        // ------------------------------------------------------ //

        if (not pt_bins.empty())
        {
            TH1* h_num = new TH1F("h_num_pt", "Numerator Counts;p_{T} (GeV)"  , num_pt_bins, pt_bins.data());
            TH1* h_den = new TH1F("h_den_pt", "Denominator Counts;p_{T} (GeV)", num_pt_bins, pt_bins.data());
            TH1* h_eff = new TH1F("h_eff_pt", "Efficiency;p_{T} (GeV)"        , num_pt_bins, pt_bins.data());
            h_eff->GetYaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray2D pt_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("pt_models"), pt_bins); 

            for (size_t pt_bin = 0; pt_bin != num_pt_bins; pt_bin++)
            {
                tnp::Result result;
                if (dataset.m_is_data)
                {
                    tnp::Model::value_type sig_pass_model = pt_models[pt_bin][0];
                    tnp::Model::value_type sig_fail_model = pt_models[pt_bin][1];
                    tnp::Model::value_type bkg_pass_model = pt_models[pt_bin][2];
                    tnp::Model::value_type bkg_fail_model = pt_models[pt_bin][3];
                    cout << Form("fitting bins: pt %lu", pt_bin) << endl; 
                    cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                    cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                    cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                    cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                    TH1* const h_pass = hc_mass[Form("h_pass_pt%lu", pt_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_pt%lu", pt_bin)];

                    TH1* const h_pass_template = hc_template[Form("h_pass_pt%lu", pt_bin)];
                    TH1* const h_fail_template = hc_template[Form("h_fail_pt%lu", pt_bin)];

                    // do the fit
                    result = PerformSimultaneousFit
                    (
                         sig_pass_model, 
                         sig_fail_model, 
                         bkg_pass_model, 
                         bkg_fail_model, 
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", pt_bins[pt_bin], pt_bins[pt_bin+1]), 
                         /*b_bin_label = */"", 
                         h_pass_template,
                         h_fail_template
                    );

                }
                else
                {
                    TH1* const h_pass = hc_mass[Form("h_pass_pt%lu", pt_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_pt%lu", pt_bin)];

                    // do the count
                    result = tnp::PerformSimpleCount
                    (
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", pt_bins[pt_bin], pt_bins[pt_bin+1]), 
                         /*b_bin_label = */"" 
                    );

                }

                // check for nan
                if (std::isnan(result.num.value)) {result.num.value = 0.0;}
                if (std::isnan(result.num.error)) {result.num.error = 0.0;}
                if (std::isnan(result.den.value)) {result.den.value = 0.0;}
                if (std::isnan(result.den.error)) {result.den.error = 0.0;}
                if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
                if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}

                // record output to histogram
                h_num->SetBinContent(pt_bin+1, std::max(result.num.value, 0.0));
                h_num->SetBinError  (pt_bin+1, std::max(result.num.error, 0.0));
                h_den->SetBinContent(pt_bin+1, std::max(result.den.value, 0.0));
                h_den->SetBinError  (pt_bin+1, std::max(result.den.error, 0.0));
                h_eff->SetBinContent(pt_bin+1, std::max(result.eff.value, 0.0));
                h_eff->SetBinError  (pt_bin+1, std::max(result.eff.error, 0.0));

                // print the fit plot
                if (not suffix.empty())
                {
                    const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_pt%lu",
                        analysis_path.c_str(),
                        output_label.c_str(),
                        GetStringFromLepton(lepton_type).c_str(),
                        GetStringFromSelection(den_selection).c_str(),
                        GetStringFromSelection(num_selection).c_str(),
                        dataset.m_name.c_str(),
                        pt_bin
                    );
                    const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 

                    PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
                    PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
                }
            }

            // add histograms to hist container
            hc.Add(h_num);
            hc.Add(h_den);
            hc.Add(h_eff);
        }

        // eta bins
        // ------------------------------------------------------ //

        if (not eta_bins.empty())
        {
            TH1* h_num = new TH1F("h_num_eta", Form("Numerator Counts;%s"  , eta_title.c_str()), num_eta_bins, eta_bins.data());
            TH1* h_den = new TH1F("h_den_eta", Form("Denominator Counts;%s", eta_title.c_str()), num_eta_bins, eta_bins.data());
            TH1* h_eff = new TH1F("h_eff_eta", Form("Efficiency;%s"        , eta_title.c_str()), num_eta_bins, eta_bins.data());
            h_eff->GetYaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray2D eta_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("eta_models"), eta_bins); 

            for (size_t eta_bin = 0; eta_bin != num_eta_bins; eta_bin++)
            {
                tnp::Result result;
                if (dataset.m_is_data)
                {
                    tnp::Model::value_type sig_pass_model = eta_models[eta_bin][0];
                    tnp::Model::value_type sig_fail_model = eta_models[eta_bin][1];
                    tnp::Model::value_type bkg_pass_model = eta_models[eta_bin][2];
                    tnp::Model::value_type bkg_fail_model = eta_models[eta_bin][3];
                    cout << Form("fitting bins: eta %lu", eta_bin) << endl; 
                    cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                    cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                    cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                    cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                    TH1* const h_pass = hc_mass[Form("h_pass_eta%lu", eta_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_eta%lu", eta_bin)];

                    TH1* const h_pass_template = hc_template[Form("h_pass_eta%lu", eta_bin)];
                    TH1* const h_fail_template = hc_template[Form("h_fail_eta%lu", eta_bin)];

                    // do the fit
                    result = PerformSimultaneousFit
                    (
                         sig_pass_model, 
                         sig_fail_model, 
                         bkg_pass_model, 
                         bkg_fail_model, 
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.2f < %s < %1.2f", eta_bins[eta_bin], eta_title.c_str(), eta_bins[eta_bin+1]), 
                         /*b_bin_label = */"", 
                         h_pass_template,
                         h_fail_template
                    );

                }
                else
                {
                    TH1* const h_pass = hc_mass[Form("h_pass_eta%lu", eta_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_eta%lu", eta_bin)];

                    // do the count
                    result = tnp::PerformSimpleCount
                    (
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.2f < %s < %1.2f", eta_bins[eta_bin], eta_title.c_str(), eta_bins[eta_bin+1]), 
                         /*b_bin_label = */"" 
                    );

                }

                // check for nan
                if (std::isnan(result.num.value)) {result.num.value = 0.0;}
                if (std::isnan(result.num.error)) {result.num.error = 0.0;}
                if (std::isnan(result.den.value)) {result.den.value = 0.0;}
                if (std::isnan(result.den.error)) {result.den.error = 0.0;}
                if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
                if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}

                // record output to histogram
                cout << "result.num = " << result.num_str() << endl;
                cout << "result.den = " << result.den_str() << endl;
                h_num->SetBinContent(eta_bin+1, std::max(result.num.value, 0.0));
                h_num->SetBinError  (eta_bin+1, std::max(result.num.error, 0.0));
                h_den->SetBinContent(eta_bin+1, std::max(result.den.value, 0.0));
                h_den->SetBinError  (eta_bin+1, std::max(result.den.error, 0.0));
                h_eff->SetBinContent(eta_bin+1, std::max(result.eff.value, 0.0));
                h_eff->SetBinError  (eta_bin+1, std::max(result.eff.error, 0.0));

                // print the fit plot
                if (not suffix.empty())
                {
                    const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_eta%lu",
                        analysis_path.c_str(),
                        output_label.c_str(),
                        GetStringFromLepton(lepton_type).c_str(),
                        GetStringFromSelection(den_selection).c_str(),
                        GetStringFromSelection(num_selection).c_str(),
                        dataset.m_name.c_str(),
                        eta_bin
                    );
                    const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 

                    PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
                    PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
                }
            }

            // add histograms to hist container
            hc.Add(h_num);
            hc.Add(h_den);
            hc.Add(h_eff);
        }

        // phi bins
        // ------------------------------------------------------ //

        if (not phi_bins.empty())
        {
            TH1* h_num = new TH1F("h_num_phi", Form("Numerator Counts;%s"  , phi_title.c_str()), num_phi_bins, phi_bins.data());
            TH1* h_den = new TH1F("h_den_phi", Form("Denominator Counts;%s", phi_title.c_str()), num_phi_bins, phi_bins.data());
            TH1* h_eff = new TH1F("h_eff_phi", Form("Efficiency;%s"        , phi_title.c_str()), num_phi_bins, phi_bins.data());
            h_eff->GetYaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray2D phi_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("phi_models"), phi_bins); 

            for (size_t phi_bin = 0; phi_bin != num_phi_bins; phi_bin++)
            {
                tnp::Result result;
                if (dataset.m_is_data)
                {
                    tnp::Model::value_type sig_pass_model = phi_models[phi_bin][0];
                    tnp::Model::value_type sig_fail_model = phi_models[phi_bin][1];
                    tnp::Model::value_type bkg_pass_model = phi_models[phi_bin][2];
                    tnp::Model::value_type bkg_fail_model = phi_models[phi_bin][3];
                    cout << Form("fitting bins: phi %lu", phi_bin) << endl; 
                    cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                    cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                    cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                    cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                    TH1* const h_pass = hc_mass[Form("h_pass_phi%lu", phi_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_phi%lu", phi_bin)];

                    TH1* const h_pass_template = hc_template[Form("h_pass_phi%lu", phi_bin)];
                    TH1* const h_fail_template = hc_template[Form("h_fail_phi%lu", phi_bin)];

                    // do the fit
                    result = PerformSimultaneousFit
                    (
                         sig_pass_model, 
                         sig_fail_model, 
                         bkg_pass_model, 
                         bkg_fail_model, 
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.2f < %s < %1.2f", phi_bins[phi_bin], phi_title.c_str(), phi_bins[phi_bin+1]), 
                         /*b_bin_label = */"", 
                         h_pass_template,
                         h_fail_template
                    );

                }
                else
                {
                    TH1* const h_pass = hc_mass[Form("h_pass_phi%lu", phi_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_phi%lu", phi_bin)];

                    // do the count
                    result = tnp::PerformSimpleCount
                    (
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.2f < %s < %1.2f", phi_bins[phi_bin], phi_title.c_str(), phi_bins[phi_bin+1]), 
                         /*b_bin_label = */"" 
                    );

                }

                // check for nan
                if (std::isnan(result.num.value)) {result.num.value = 0.0;}
                if (std::isnan(result.num.error)) {result.num.error = 0.0;}
                if (std::isnan(result.den.value)) {result.den.value = 0.0;}
                if (std::isnan(result.den.error)) {result.den.error = 0.0;}
                if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
                if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}

                // record output to histogram
                h_num->SetBinContent(phi_bin+1, result.num.value);
                h_num->SetBinError  (phi_bin+1, result.num.error);
                h_den->SetBinContent(phi_bin+1, result.den.value);
                h_den->SetBinError  (phi_bin+1, result.den.error);
                h_eff->SetBinContent(phi_bin+1, result.eff.value);
                h_eff->SetBinError  (phi_bin+1, result.eff.error);

                // print the fit plot
                if (not suffix.empty())
                {
                    const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_phi%lu",
                        analysis_path.c_str(),
                        output_label.c_str(),
                        GetStringFromLepton(lepton_type).c_str(),
                        GetStringFromSelection(den_selection).c_str(),
                        GetStringFromSelection(num_selection).c_str(),
                        dataset.m_name.c_str(),
                        phi_bin
                    );
                    const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 

                    PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
                    PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
                }
            }

            // add histograms to hist container
            hc.Add(h_num);
            hc.Add(h_den);
            hc.Add(h_eff);
        }

        // nvtx bins
        // ------------------------------------------------------ //

        if (not nvtx_bins.empty())
        {
            TH1* h_num = new TH1F("h_num_nvtx", "Numerator Counts;# vertices"  , num_nvtx_bins, nvtx_bins.data());
            TH1* h_den = new TH1F("h_den_nvtx", "Denominator Counts;# vertices", num_nvtx_bins, nvtx_bins.data());
            TH1* h_eff = new TH1F("h_eff_nvtx", "Efficiency;# vertices"        , num_nvtx_bins, nvtx_bins.data());
            h_eff->GetYaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray2D nvtx_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("nvtx_models"), nvtx_bins); 

            for (size_t nvtx_bin = 0; nvtx_bin != num_nvtx_bins; nvtx_bin++)
            {
                tnp::Result result;
                if (dataset.m_is_data)
                {
                    tnp::Model::value_type sig_pass_model = nvtx_models[nvtx_bin][0];
                    tnp::Model::value_type sig_fail_model = nvtx_models[nvtx_bin][1];
                    tnp::Model::value_type bkg_pass_model = nvtx_models[nvtx_bin][2];
                    tnp::Model::value_type bkg_fail_model = nvtx_models[nvtx_bin][3];
                    cout << Form("fitting bins: nvtx %lu", nvtx_bin) << endl; 
                    cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                    cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                    cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                    cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                    TH1* const h_pass = hc_mass[Form("h_pass_nvtx%lu", nvtx_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_nvtx%lu", nvtx_bin)];

                    TH1* const h_pass_template = hc_template[Form("h_pass_nvtx%lu", nvtx_bin)];
                    TH1* const h_fail_template = hc_template[Form("h_fail_nvtx%lu", nvtx_bin)];

                    // do the fit
                    result = PerformSimultaneousFit
                    (
                         sig_pass_model, 
                         sig_fail_model, 
                         bkg_pass_model, 
                         bkg_fail_model, 
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.0f < # vertices < %1.0f", nvtx_bins[nvtx_bin], nvtx_bins[nvtx_bin+1]), 
                         /*b_bin_label = */"", 
                         h_pass_template,
                         h_fail_template
                    );

                }
                else
                {
                    TH1* const h_pass = hc_mass[Form("h_pass_nvtx%lu", nvtx_bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_nvtx%lu", nvtx_bin)];

                    // do the count
                    result = tnp::PerformSimpleCount
                    (
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */Form("%1.0f < # vertices < %1.0f", nvtx_bins[nvtx_bin], nvtx_bins[nvtx_bin+1]), 
                         /*b_bin_label = */"" 
                    );

                }

                // check for nan
                if (std::isnan(result.num.value)) {result.num.value = 0.0;}
                if (std::isnan(result.num.error)) {result.num.error = 0.0;}
                if (std::isnan(result.den.value)) {result.den.value = 0.0;}
                if (std::isnan(result.den.error)) {result.den.error = 0.0;}
                if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
                if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}

                // record output to histogram
                h_num->SetBinContent(nvtx_bin+1, result.num.value);
                h_num->SetBinError  (nvtx_bin+1, result.num.error);
                h_den->SetBinContent(nvtx_bin+1, result.den.value);
                h_den->SetBinError  (nvtx_bin+1, result.den.error);
                h_eff->SetBinContent(nvtx_bin+1, result.eff.value);
                h_eff->SetBinError  (nvtx_bin+1, result.eff.error);

                // print the fit plot
                if (not suffix.empty())
                {
                    const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_nvtx%lu",
                        analysis_path.c_str(),
                        output_label.c_str(),
                        GetStringFromLepton(lepton_type).c_str(),
                        GetStringFromSelection(den_selection).c_str(),
                        GetStringFromSelection(num_selection).c_str(),
                        dataset.m_name.c_str(),
                        nvtx_bin
                    );
                    const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 

                    PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
                    PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
                }
            }

            // add histograms to hist container
            hc.Add(h_num);
            hc.Add(h_den);
            hc.Add(h_eff);
        }

        // pt vs eta bins
        // ------------------------------------------------------ //

        if (not (pt_bins.empty() or eta_bins.empty()))
        {
            TH1* h_num = new TH2F("h_num_pt_vs_eta", Form("Numerator Counts;%s;p_{T} (GeV)"  , eta_title.c_str()), num_eta_bins, eta_bins.data(), num_pt_bins, pt_bins.data());
            TH1* h_den = new TH2F("h_den_pt_vs_eta", Form("Denominator Counts;%s;p_{T} (GeV)", eta_title.c_str()), num_eta_bins, eta_bins.data(), num_pt_bins, pt_bins.data());
            TH1* h_eff = new TH2F("h_eff_pt_vs_eta", Form("Efficiency;%s;p_{T} (GeV)"        , eta_title.c_str()), num_eta_bins, eta_bins.data(), num_pt_bins, pt_bins.data());
            h_eff->GetZaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray3D pt_vs_eta_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("pt_vs_eta_models"), pt_bins, eta_bins); 

            for (size_t pt_bin = 0; pt_bin != num_pt_bins; pt_bin++)
            {
                for (size_t eta_bin = 0; eta_bin != num_eta_bins; eta_bin++)
                {
                    tnp::Result result;
                    if (dataset.m_is_data)
                    {
                        tnp::Model::value_type sig_pass_model = pt_vs_eta_models[pt_bin][eta_bin][0];
                        tnp::Model::value_type sig_fail_model = pt_vs_eta_models[pt_bin][eta_bin][1];
                        tnp::Model::value_type bkg_pass_model = pt_vs_eta_models[pt_bin][eta_bin][2];
                        tnp::Model::value_type bkg_fail_model = pt_vs_eta_models[pt_bin][eta_bin][3];
                        cout << Form("fitting bins: pt %lu, eta %lu ", pt_bin, eta_bin) << endl; 
                        cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                        cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                        cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                        cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                        TH1* const h_pass = hc_mass[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
                        TH1* const h_fail = hc_mass[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];

                        TH1* const h_pass_template = hc_template[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
                        TH1* const h_fail_template = hc_template[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];

                        // do the fit
                        result = PerformSimultaneousFit
                        (
                             sig_pass_model, 
                             sig_fail_model, 
                             bkg_pass_model, 
                             bkg_fail_model, 
                             h_pass, 
                             h_fail,
                             mass_low,
                             mass_high,
                             mass_bin_width,
                             /*a_bin_label = */Form("%1.2f < %s < %1.2f", eta_bins[eta_bin], eta_title.c_str(), eta_bins[eta_bin+1]), 
                             /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", pt_bins[pt_bin], pt_bins[pt_bin+1]), 
                             h_pass_template,
                             h_fail_template
                        );
                    }
                    else
                    {
                        TH1* const h_pass = hc_mass[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
                        TH1* const h_fail = hc_mass[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];

                        // do the count
                        result = tnp::PerformSimpleCount
                        (
                             h_pass, 
                             h_fail,
                             mass_low,
                             mass_high,
                             mass_bin_width,
                             /*a_bin_label = */Form("%1.2f < %s < %1.2f", eta_bins[eta_bin], eta_title.c_str(), eta_bins[eta_bin+1]), 
                             /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", pt_bins[pt_bin], pt_bins[pt_bin+1]) 
                        );

                    }

                    // check for nan
                    if (std::isnan(result.num.value)) {result.num.value = 0.0;}
                    if (std::isnan(result.num.error)) {result.num.error = 0.0;}
                    if (std::isnan(result.den.value)) {result.den.value = 0.0;}
                    if (std::isnan(result.den.error)) {result.den.error = 0.0;}
                    if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
                    if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}

                    // record output to histogram
                    h_num->SetBinContent(eta_bin+1, pt_bin+1, result.num.value);
                    h_num->SetBinError  (eta_bin+1, pt_bin+1, result.num.error);
                    h_den->SetBinContent(eta_bin+1, pt_bin+1, result.den.value);
                    h_den->SetBinError  (eta_bin+1, pt_bin+1, result.den.error);
                    h_eff->SetBinContent(eta_bin+1, pt_bin+1, result.eff.value);
                    h_eff->SetBinError  (eta_bin+1, pt_bin+1, result.eff.error);

                    // print the fit plot
                    if (not suffix.empty())
                    {
                        const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_pt%lu_vs_eta%lu",
                            analysis_path.c_str(),
                            output_label.c_str(),
                            GetStringFromLepton(lepton_type).c_str(),
                            GetStringFromSelection(den_selection).c_str(),
                            GetStringFromSelection(num_selection).c_str(),
                            dataset.m_name.c_str(),
                            pt_bin,
                            eta_bin
                        );
                        const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 

                        PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
                        PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
                    }
                }
            }

            // add histograms to hist container
            hc.Add(h_num);
            hc.Add(h_den);
            hc.Add(h_eff);
        }

        // eta vs phi bins
        // ------------------------------------------------------ //

        if (not (eta_bins.empty() or phi_bins.empty()))
        {
            TH1* h_num = new TH2F("h_num_eta_vs_phi", Form("Numerator Counts;%s;p_{T} (GeV)"  , phi_title.c_str()), num_phi_bins, phi_bins.data(), num_eta_bins, eta_bins.data());
            TH1* h_den = new TH2F("h_den_eta_vs_phi", Form("Denominator Counts;%s;p_{T} (GeV)", phi_title.c_str()), num_phi_bins, phi_bins.data(), num_eta_bins, eta_bins.data());
            TH1* h_eff = new TH2F("h_eff_eta_vs_phi", Form("Efficiency;%s;p_{T} (GeV)"        , phi_title.c_str()), num_phi_bins, phi_bins.data(), num_eta_bins, eta_bins.data());
            h_eff->GetZaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray3D eta_vs_phi_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("eta_vs_phi_models"), eta_bins, phi_bins); 

            for (size_t eta_bin = 0; eta_bin != num_eta_bins; eta_bin++)
            {
                for (size_t phi_bin = 0; phi_bin != num_phi_bins; phi_bin++)
                {
                    tnp::Result result;
                    if (dataset.m_is_data)
                    {
                        tnp::Model::value_type sig_pass_model = eta_vs_phi_models[eta_bin][phi_bin][0];
                        tnp::Model::value_type sig_fail_model = eta_vs_phi_models[eta_bin][phi_bin][1];
                        tnp::Model::value_type bkg_pass_model = eta_vs_phi_models[eta_bin][phi_bin][2];
                        tnp::Model::value_type bkg_fail_model = eta_vs_phi_models[eta_bin][phi_bin][3];
                        cout << Form("fitting bins: eta %lu, phi %lu ", eta_bin, phi_bin) << endl; 
                        cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                        cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                        cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                        cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                        TH1* const h_pass = hc_mass[Form("h_pass_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
                        TH1* const h_fail = hc_mass[Form("h_fail_eta%lu_vs_phi%lu", eta_bin, phi_bin)];

                        TH1* const h_pass_template = hc_template[Form("h_pass_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
                        TH1* const h_fail_template = hc_template[Form("h_fail_eta%lu_vs_phi%lu", eta_bin, phi_bin)];

                        // do the fit
                        result = PerformSimultaneousFit
                        (
                             sig_pass_model, 
                             sig_fail_model, 
                             bkg_pass_model, 
                             bkg_fail_model, 
                             h_pass, 
                             h_fail,
                             mass_low,
                             mass_high,
                             mass_bin_width,
                             /*a_bin_label = */Form("%1.2f < %s < %1.2f", phi_bins[phi_bin], phi_title.c_str(), phi_bins[phi_bin+1]), 
                             /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", eta_bins[eta_bin], eta_bins[eta_bin+1]), 
                             h_pass_template,
                             h_fail_template
                        );

                    }
                    else
                    {
                        TH1* const h_pass = hc_mass[Form("h_pass_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
                        TH1* const h_fail = hc_mass[Form("h_fail_eta%lu_vs_phi%lu", eta_bin, phi_bin)];

                        // do the count
                        result = tnp::PerformSimpleCount
                        (
                             h_pass, 
                             h_fail,
                             mass_low,
                             mass_high,
                             mass_bin_width,
                             /*a_bin_label = */Form("%1.2f < %s < %1.2f", phi_bins[phi_bin], phi_title.c_str(), phi_bins[phi_bin+1]), 
                             /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", eta_bins[eta_bin], eta_bins[eta_bin+1]) 
                        );

                    }

                    // check for nan
                    if (std::isnan(result.num.value)) {result.num.value = 0.0;}
                    if (std::isnan(result.num.error)) {result.num.error = 0.0;}
                    if (std::isnan(result.den.value)) {result.den.value = 0.0;}
                    if (std::isnan(result.den.error)) {result.den.error = 0.0;}
                    if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
                    if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}

                    // record output to histogram
                    h_num->SetBinContent(phi_bin+1, eta_bin+1, result.num.value);
                    h_num->SetBinError  (phi_bin+1, eta_bin+1, result.num.error);
                    h_den->SetBinContent(phi_bin+1, eta_bin+1, result.den.value);
                    h_den->SetBinError  (phi_bin+1, eta_bin+1, result.den.error);
                    h_eff->SetBinContent(phi_bin+1, eta_bin+1, result.eff.value);
                    h_eff->SetBinError  (phi_bin+1, eta_bin+1, result.eff.error);

                    // print the fit plot
                    if (not suffix.empty())
                    {
                        const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_eta%lu_vs_phi%lu",
                            analysis_path.c_str(),
                            output_label.c_str(),
                            GetStringFromLepton(lepton_type).c_str(),
                            GetStringFromSelection(den_selection).c_str(),
                            GetStringFromSelection(num_selection).c_str(),
                            dataset.m_name.c_str(),
                            eta_bin,
                            phi_bin
                        );
                        const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 

                        PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
                        PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
                    }
                }
            }

            // add histograms to hist container
            hc.Add(h_num);
            hc.Add(h_den);
            hc.Add(h_eff);
        }

        // write and print the output
        // ------------------------------------------------------ //

        rt::SetStyle();
        hc.Write(eff_plot_file_name);
        if (not suffix.empty())
        {
            std::string output_print_path = lt::string_replace_all(eff_plot_file_name, ".root", "") + "/" + suffix;
            cout << "Printing histograms to: " << output_print_path  << endl;
            hc.Print(output_print_path, suffix);
        }

    } // end dataset loop

    // done
    return 0;
}
catch (std::exception& e)
{
    cerr << "[tnp_eff_plots] Error: failed..." << endl;
    cerr << e.what() << endl;
    return 1;
}
