// C++ includes
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>

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
    const std::vector<double>& bins
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
    const std::vector<double>& a_bins,
    const std::vector<double>& b_bins
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
//     const unsigned int num_pt_bins      = pt_bins.size()-1;
//     const unsigned int num_eta_bins     = eta_bins.size()-1;
//     const unsigned int num_phi_bins     = phi_bins.size()-1;
//     const unsigned int num_nvtx_bins    = nvtx_bins.size()-1;

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

        // efficiency hists container
        rt::TH1Container hc;

        // simple struct to hold the 1D binning info
        // ------------------------------------------------------ //

        struct bin_info_t
        {
            // construct:
            bin_info_t() {}

            // construct:
            bin_info_t
            (
                const std::string& n,
                const std::string& t,
                const std::string& u,
                const size_t nd,
                const std::vector<double>& b
            )
                : name(n)
                , title(t)
                , unit(u)
                , num_digits(nd)
                , bins(b)
            {}

            // methods
            const size_t nbins()           const {return bins.size()-1;}
            const double* data()           const {return bins.data();}
            const std::string model_name() const {return name + "_models";}
            const std::string bin_label(const size_t bin) const
            {
                using namespace std;

                if (unit.empty())
                {
                    ostringstream os;
                    os << setprecision(num_digits) << bins.at(bin) << " < " << title << " < " << bins.at(bin+1);
                    return os.str(); 
                }
                else
                {
                    ostringstream os;
                    os << setprecision(num_digits) << bins.at(bin) << " " << unit << " < " << title << " < " << bins.at(bin+1) << " " << unit;
                    return os.str(); 
                }
            }

            // members:
            std::string name;
            std::string title;
            std::string unit;
            size_t num_digits;
            std::vector<double> bins;
        };

        const bin_info_t pt_bin_info  ("pt"  , "p_{T}"                                                            , "(GeV)", 2, pt_bins  );
        const bin_info_t eta_bin_info ("eta" , (not eta_bins.empty() and eta_bins.front() < 0 ? "#eta" : "|#eta|"), ""     , 3, eta_bins );
        const bin_info_t phi_bin_info ("phi" , (not phi_bins.empty() and phi_bins.front() < 0 ? "#phi" : "|#phi|"), ""     , 3, phi_bins );
        const bin_info_t nvtx_bin_info("nvtx", "# of vertices"                                                    , ""     , 2, nvtx_bins);

        std::vector<bin_info_t> bin_infos;
        if (not pt_bins.empty()  ) {bin_infos.push_back(pt_bin_info  );}
        if (not eta_bins.empty() ) {bin_infos.push_back(eta_bin_info );}
        if (not phi_bins.empty() ) {bin_infos.push_back(phi_bin_info );}
        if (not nvtx_bins.empty()) {bin_infos.push_back(nvtx_bin_info);}

        // create 1D bin efficiency plots
        // ------------------------------------------------------ //

        for (size_t i = 0; i != bin_infos.size(); i++)
        {
            const bin_info_t& bi = bin_infos.at(i);

            // book hists
            TH1* h_num = new TH1F(Form("h_num_%s", bi.name.c_str()), Form("Numerator Counts;%s"  , bi.title.c_str()), bi.nbins(), bi.data());
            TH1* h_den = new TH1F(Form("h_den_%s", bi.name.c_str()), Form("Denominator Counts;%s", bi.title.c_str()), bi.nbins(), bi.data());
            TH1* h_eff = new TH1F(Form("h_eff_%s", bi.name.c_str()), Form("Efficiency;%s"        , bi.title.c_str()), bi.nbins(), bi.data());
            h_eff->GetYaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray2D models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >(bi.model_name()), bi.bins); 

            for (size_t bin = 0; bin != bi.nbins(); bin++)
            {
                tnp::Result result;
                if (dataset.m_is_data)
                {
                    tnp::Model::value_type sig_pass_model = models[bin][0];
                    tnp::Model::value_type sig_fail_model = models[bin][1];
                    tnp::Model::value_type bkg_pass_model = models[bin][2];
                    tnp::Model::value_type bkg_fail_model = models[bin][3];
                    cout << Form("[tnp_eff_plots] fitting bins: %s %lu", bi.name.c_str(), bin) << endl; 
                    cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                    cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                    cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                    cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                    TH1* const h_pass = hc_mass[Form("h_pass_%s%lu", bi.name.c_str(), bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_%s%lu", bi.name.c_str(), bin)];

                    TH1* const h_pass_template = hc_template[Form("h_pass_%s%lu", bi.name.c_str(), bin)];
                    TH1* const h_fail_template = hc_template[Form("h_fail_%s%lu", bi.name.c_str(), bin)];

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
                         /*a_bin_label = */bi.bin_label(bin),
                         /*b_bin_label = */"", 
                         h_pass_template,
                         h_fail_template
                    );

                }
                else
                {
                    TH1* const h_pass = hc_mass[Form("h_pass_%s%lu", bi.name.c_str(), bin)];
                    TH1* const h_fail = hc_mass[Form("h_fail_%s%lu", bi.name.c_str(), bin)];

                    // do the count
                    result = tnp::PerformSimpleCount
                    (
                         h_pass, 
                         h_fail,
                         mass_low,
                         mass_high,
                         mass_bin_width,
                         /*a_bin_label = */bi.bin_label(bin),
                         /*b_bin_label = */"" 
                    );
                }

                // check results for nan
                if (std::isnan(result.num.value)) {result.num.value = 0.0;}
                if (std::isnan(result.num.error)) {result.num.error = 0.0;}
                if (std::isnan(result.den.value)) {result.den.value = 0.0;}
                if (std::isnan(result.den.error)) {result.den.error = 0.0;}
                if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
                if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}

                // record results to histogram
                h_num->SetBinContent(bin+1, std::max(result.num.value, 0.0));
                h_num->SetBinError  (bin+1, std::max(result.num.error, 0.0));
                h_den->SetBinContent(bin+1, std::max(result.den.value, 0.0));
                h_den->SetBinError  (bin+1, std::max(result.den.error, 0.0));
                h_eff->SetBinContent(bin+1, std::max(result.eff.value, 0.0));
                h_eff->SetBinError  (bin+1, std::max(result.eff.error, 0.0));

                // print the fit plot
                if (not suffix.empty())
                {
                    const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_%s%lu",
                        analysis_path.c_str(),
                        output_label.c_str(),
                        GetStringFromLepton(lepton_type).c_str(),
                        GetStringFromSelection(den_selection).c_str(),
                        GetStringFromSelection(num_selection).c_str(),
                        dataset.m_name.c_str(),
                        bi.name.c_str(),
                        bin
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

        // simple struct to hold the 2D binning info
        // ------------------------------------------------------ //

        struct bin_info_pair_t
        {
            bin_info_pair_t() {}
            bin_info_pair_t(const bin_info_t& xbi, const bin_info_t& ybi) : x(xbi), y(ybi) {}
            bin_info_t x;
            bin_info_t y;
        };

        std::vector<bin_info_pair_t> bin_info_pairs;
        if (not pt_bins.empty()  and not eta_bins.empty()) {bin_info_pairs.push_back(bin_info_pair_t(eta_bin_info, pt_bin_info ));}
        if (not eta_bins.empty() and not phi_bins.empty()) {bin_info_pairs.push_back(bin_info_pair_t(phi_bin_info, eta_bin_info));}

        // create 2D bin efficiency plots
        // ------------------------------------------------------ //

        for (size_t i = 0; i != bin_info_pairs.size(); i++)
        {
            const bin_info_t& xbi = bin_info_pairs.at(i).x;
            const bin_info_t& ybi = bin_info_pairs.at(i).y;

            TH1* h_num = new TH2F(Form("h_num_%s_vs_%s", ybi.name.c_str(), xbi.name.c_str()), Form("Numerator Counts;%s;%s"  , xbi.title.c_str(), ybi.title.c_str()), xbi.nbins(), xbi.data(), ybi.nbins(), ybi.data());
            TH1* h_den = new TH2F(Form("h_den_%s_vs_%s", ybi.name.c_str(), xbi.name.c_str()), Form("Denominator Counts;%s;%s", xbi.title.c_str(), ybi.title.c_str()), xbi.nbins(), xbi.data(), ybi.nbins(), ybi.data());
            TH1* h_eff = new TH2F(Form("h_eff_%s_vs_%s", ybi.name.c_str(), xbi.name.c_str()), Form("Efficiency;%s;%s"        , xbi.title.c_str(), ybi.title.c_str()), xbi.nbins(), xbi.data(), ybi.nbins(), ybi.data());
            h_eff->GetZaxis()->SetRangeUser(0, 1.1);

            // models
            const ModelArray3D models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >(Form("%s_vs_%s_models", ybi.name.c_str(), xbi.name.c_str())), xbi.bins, ybi.bins); 

            for (size_t xbin = 0; xbin != xbi.nbins(); xbin++)
            {
                for (size_t ybin = 0; ybin != ybi.nbins(); ybin++)
                {
                    tnp::Result result;
                    if (dataset.m_is_data)
                    {
                        tnp::Model::value_type sig_pass_model = models[xbin][ybin][0];
                        tnp::Model::value_type sig_fail_model = models[xbin][ybin][1];
                        tnp::Model::value_type bkg_pass_model = models[xbin][ybin][2];
                        tnp::Model::value_type bkg_fail_model = models[xbin][ybin][3];
                        cout << Form("[tnp_eff_plots] fitting bins: %s %lu, %s %lu ", xbi.name.c_str(), xbin, ybi.name.c_str(), ybin) << endl; 
                        cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
                        cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
                        cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
                        cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

                        TH1* const h_pass = hc_mass[Form("h_pass_%s%lu_vs_%s%lu", ybi.name.c_str(), ybin, xbi.name.c_str(), xbin)];
                        TH1* const h_fail = hc_mass[Form("h_fail_%s%lu_vs_%s%lu", ybi.name.c_str(), ybin, xbi.name.c_str(), xbin)];

                        TH1* const h_pass_template = hc_template[Form("h_pass_%s%lu_vs_%s%lu", ybi.name.c_str(), ybin, xbi.name.c_str(), xbin)];
                        TH1* const h_fail_template = hc_template[Form("h_fail_%s%lu_vs_%s%lu", ybi.name.c_str(), ybin, xbi.name.c_str(), xbin)];

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
                             xbi.bin_label(xbin),
                             ybi.bin_label(ybin),
                             h_pass_template,
                             h_fail_template
                        );
                    }
                    else
                    {
                        TH1* const h_pass = hc_mass[Form("h_pass_%s%lu_vs_%s%lu", ybi.name.c_str(), ybin, xbi.name.c_str(), xbin)];
                        TH1* const h_fail = hc_mass[Form("h_fail_%s%lu_vs_%s%lu", ybi.name.c_str(), ybin, xbi.name.c_str(), xbin)];

                        // do the count
                        result = tnp::PerformSimpleCount
                        (
                             h_pass, 
                             h_fail,
                             mass_low,
                             mass_high,
                             mass_bin_width,
                             xbi.bin_label(xbin),
                             ybi.bin_label(ybin)
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
                    h_num->SetBinContent(xbin+1, ybin+1, result.num.value);
                    h_num->SetBinError  (xbin+1, ybin+1, result.num.error);
                    h_den->SetBinContent(xbin+1, ybin+1, result.den.value);
                    h_den->SetBinError  (xbin+1, ybin+1, result.den.error);
                    h_eff->SetBinContent(xbin+1, ybin+1, result.eff.value);
                    h_eff->SetBinError  (xbin+1, ybin+1, result.eff.error);

                    // print the fit plot
                    if (not suffix.empty())
                    {
                        const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_%s%lu_vs_%s%lu",
                            analysis_path.c_str(),
                            output_label.c_str(),
                            GetStringFromLepton(lepton_type).c_str(),
                            GetStringFromSelection(den_selection).c_str(),
                            GetStringFromSelection(num_selection).c_str(),
                            dataset.m_name.c_str(),
                            ybi.name.c_str(),
                            ybin,
                            xbi.name.c_str(),
                            xbin
                        );
                        const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 

                        PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
                        PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
                    }
                }
            }
        }


        // pt vs eta bins
        // ------------------------------------------------------ //

//         if (not (pt_bins.empty() or eta_bins.empty()))
//         {
//             TH1* h_num = new TH2F("h_num_pt_vs_eta", Form("Numerator Counts;%s;p_{T} (GeV)"  , eta_title.c_str()), num_eta_bins, eta_bins.data(), num_pt_bins, pt_bins.data());
//             TH1* h_den = new TH2F("h_den_pt_vs_eta", Form("Denominator Counts;%s;p_{T} (GeV)", eta_title.c_str()), num_eta_bins, eta_bins.data(), num_pt_bins, pt_bins.data());
//             TH1* h_eff = new TH2F("h_eff_pt_vs_eta", Form("Efficiency;%s;p_{T} (GeV)"        , eta_title.c_str()), num_eta_bins, eta_bins.data(), num_pt_bins, pt_bins.data());
//             h_eff->GetZaxis()->SetRangeUser(0, 1.1);
// 
//             // models
//             const ModelArray3D pt_vs_eta_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("pt_vs_eta_models"), pt_bins, eta_bins); 
// 
//             for (size_t pt_bin = 0; pt_bin != num_pt_bins; pt_bin++)
//             {
//                 for (size_t eta_bin = 0; eta_bin != num_eta_bins; eta_bin++)
//                 {
//                     tnp::Result result;
//                     if (dataset.m_is_data)
//                     {
//                         tnp::Model::value_type sig_pass_model = pt_vs_eta_models[pt_bin][eta_bin][0];
//                         tnp::Model::value_type sig_fail_model = pt_vs_eta_models[pt_bin][eta_bin][1];
//                         tnp::Model::value_type bkg_pass_model = pt_vs_eta_models[pt_bin][eta_bin][2];
//                         tnp::Model::value_type bkg_fail_model = pt_vs_eta_models[pt_bin][eta_bin][3];
//                         cout << Form("fitting bins: pt %lu, eta %lu ", pt_bin, eta_bin) << endl; 
//                         cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
//                         cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
//                         cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
//                         cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 
// 
//                         TH1* const h_pass = hc_mass[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
//                         TH1* const h_fail = hc_mass[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
// 
//                         TH1* const h_pass_template = hc_template[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
//                         TH1* const h_fail_template = hc_template[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
// 
//                         // do the fit
//                         result = PerformSimultaneousFit
//                         (
//                              sig_pass_model, 
//                              sig_fail_model, 
//                              bkg_pass_model, 
//                              bkg_fail_model, 
//                              h_pass, 
//                              h_fail,
//                              mass_low,
//                              mass_high,
//                              mass_bin_width,
//                              /*a_bin_label = */Form("%1.2f < %s < %1.2f", eta_bins[eta_bin], eta_title.c_str(), eta_bins[eta_bin+1]), 
//                              /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", pt_bins[pt_bin], pt_bins[pt_bin+1]), 
//                              h_pass_template,
//                              h_fail_template
//                         );
//                     }
//                     else
//                     {
//                         TH1* const h_pass = hc_mass[Form("h_pass_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
//                         TH1* const h_fail = hc_mass[Form("h_fail_pt%lu_vs_eta%lu", pt_bin, eta_bin)];
// 
//                         // do the count
//                         result = tnp::PerformSimpleCount
//                         (
//                              h_pass, 
//                              h_fail,
//                              mass_low,
//                              mass_high,
//                              mass_bin_width,
//                              /*a_bin_label = */Form("%1.2f < %s < %1.2f", eta_bins[eta_bin], eta_title.c_str(), eta_bins[eta_bin+1]), 
//                              /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", pt_bins[pt_bin], pt_bins[pt_bin+1]) 
//                         );
// 
//                     }
// 
//                     // check for nan
//                     if (std::isnan(result.num.value)) {result.num.value = 0.0;}
//                     if (std::isnan(result.num.error)) {result.num.error = 0.0;}
//                     if (std::isnan(result.den.value)) {result.den.value = 0.0;}
//                     if (std::isnan(result.den.error)) {result.den.error = 0.0;}
//                     if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
//                     if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}
// 
//                     // record output to histogram
//                     h_num->SetBinContent(eta_bin+1, pt_bin+1, result.num.value);
//                     h_num->SetBinError  (eta_bin+1, pt_bin+1, result.num.error);
//                     h_den->SetBinContent(eta_bin+1, pt_bin+1, result.den.value);
//                     h_den->SetBinError  (eta_bin+1, pt_bin+1, result.den.error);
//                     h_eff->SetBinContent(eta_bin+1, pt_bin+1, result.eff.value);
//                     h_eff->SetBinError  (eta_bin+1, pt_bin+1, result.eff.error);
// 
//                     // print the fit plot
//                     if (not suffix.empty())
//                     {
//                         const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_pt%lu_vs_eta%lu",
//                             analysis_path.c_str(),
//                             output_label.c_str(),
//                             GetStringFromLepton(lepton_type).c_str(),
//                             GetStringFromSelection(den_selection).c_str(),
//                             GetStringFromSelection(num_selection).c_str(),
//                             dataset.m_name.c_str(),
//                             pt_bin,
//                             eta_bin
//                         );
//                         const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 
// 
//                         PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
//                         PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
//                     }
//                 }
//             }
// 
//             // add histograms to hist container
//             hc.Add(h_num);
//             hc.Add(h_den);
//             hc.Add(h_eff);
//         }
// 
//         // eta vs phi bins
//         // ------------------------------------------------------ //
// 
//         if (not (eta_bins.empty() or phi_bins.empty()))
//         {
//             TH1* h_num = new TH2F("h_num_eta_vs_phi", Form("Numerator Counts;%s;p_{T} (GeV)"  , phi_title.c_str()), num_phi_bins, phi_bins.data(), num_eta_bins, eta_bins.data());
//             TH1* h_den = new TH2F("h_den_eta_vs_phi", Form("Denominator Counts;%s;p_{T} (GeV)", phi_title.c_str()), num_phi_bins, phi_bins.data(), num_eta_bins, eta_bins.data());
//             TH1* h_eff = new TH2F("h_eff_eta_vs_phi", Form("Efficiency;%s;p_{T} (GeV)"        , phi_title.c_str()), num_phi_bins, phi_bins.data(), num_eta_bins, eta_bins.data());
//             h_eff->GetZaxis()->SetRangeUser(0, 1.1);
// 
//             // models
//             const ModelArray3D eta_vs_phi_models = GetModelArrayFromVString(tnp_cfg.getParameter<std::vector<std::string> >("eta_vs_phi_models"), eta_bins, phi_bins); 
// 
//             for (size_t eta_bin = 0; eta_bin != num_eta_bins; eta_bin++)
//             {
//                 for (size_t phi_bin = 0; phi_bin != num_phi_bins; phi_bin++)
//                 {
//                     tnp::Result result;
//                     if (dataset.m_is_data)
//                     {
//                         tnp::Model::value_type sig_pass_model = eta_vs_phi_models[eta_bin][phi_bin][0];
//                         tnp::Model::value_type sig_fail_model = eta_vs_phi_models[eta_bin][phi_bin][1];
//                         tnp::Model::value_type bkg_pass_model = eta_vs_phi_models[eta_bin][phi_bin][2];
//                         tnp::Model::value_type bkg_fail_model = eta_vs_phi_models[eta_bin][phi_bin][3];
//                         cout << Form("fitting bins: eta %lu, phi %lu ", eta_bin, phi_bin) << endl; 
//                         cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
//                         cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
//                         cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
//                         cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 
// 
//                         TH1* const h_pass = hc_mass[Form("h_pass_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
//                         TH1* const h_fail = hc_mass[Form("h_fail_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
// 
//                         TH1* const h_pass_template = hc_template[Form("h_pass_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
//                         TH1* const h_fail_template = hc_template[Form("h_fail_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
// 
//                         // do the fit
//                         result = PerformSimultaneousFit
//                         (
//                              sig_pass_model, 
//                              sig_fail_model, 
//                              bkg_pass_model, 
//                              bkg_fail_model, 
//                              h_pass, 
//                              h_fail,
//                              mass_low,
//                              mass_high,
//                              mass_bin_width,
//                              /*a_bin_label = */Form("%1.2f < %s < %1.2f", phi_bins[phi_bin], phi_title.c_str(), phi_bins[phi_bin+1]), 
//                              /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", eta_bins[eta_bin], eta_bins[eta_bin+1]), 
//                              h_pass_template,
//                              h_fail_template
//                         );
// 
//                     }
//                     else
//                     {
//                         TH1* const h_pass = hc_mass[Form("h_pass_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
//                         TH1* const h_fail = hc_mass[Form("h_fail_eta%lu_vs_phi%lu", eta_bin, phi_bin)];
// 
//                         // do the count
//                         result = tnp::PerformSimpleCount
//                         (
//                              h_pass, 
//                              h_fail,
//                              mass_low,
//                              mass_high,
//                              mass_bin_width,
//                              /*a_bin_label = */Form("%1.2f < %s < %1.2f", phi_bins[phi_bin], phi_title.c_str(), phi_bins[phi_bin+1]), 
//                              /*b_bin_label = */Form("%1.0f GeV < p_{T} < %1.0f GeV", eta_bins[eta_bin], eta_bins[eta_bin+1]) 
//                         );
// 
//                     }
// 
//                     // check for nan
//                     if (std::isnan(result.num.value)) {result.num.value = 0.0;}
//                     if (std::isnan(result.num.error)) {result.num.error = 0.0;}
//                     if (std::isnan(result.den.value)) {result.den.value = 0.0;}
//                     if (std::isnan(result.den.error)) {result.den.error = 0.0;}
//                     if (std::isnan(result.eff.value)) {result.eff.value = 0.0;}
//                     if (std::isnan(result.eff.error)) {result.eff.error = 0.0;}
// 
//                     // record output to histogram
//                     h_num->SetBinContent(phi_bin+1, eta_bin+1, result.num.value);
//                     h_num->SetBinError  (phi_bin+1, eta_bin+1, result.num.error);
//                     h_den->SetBinContent(phi_bin+1, eta_bin+1, result.den.value);
//                     h_den->SetBinError  (phi_bin+1, eta_bin+1, result.den.error);
//                     h_eff->SetBinContent(phi_bin+1, eta_bin+1, result.eff.value);
//                     h_eff->SetBinError  (phi_bin+1, eta_bin+1, result.eff.error);
// 
//                     // print the fit plot
//                     if (not suffix.empty())
//                     {
//                         const std::string fit_plot_pass_name = Form("%s/plots/%s/%s/%s_%s/%s_eff/p_pass_eta%lu_vs_phi%lu",
//                             analysis_path.c_str(),
//                             output_label.c_str(),
//                             GetStringFromLepton(lepton_type).c_str(),
//                             GetStringFromSelection(den_selection).c_str(),
//                             GetStringFromSelection(num_selection).c_str(),
//                             dataset.m_name.c_str(),
//                             eta_bin,
//                             phi_bin
//                         );
//                         const std::string fit_plot_fail_name = lt::string_replace_last(fit_plot_pass_name, "pass", "fail"); 
// 
//                         PrintCanvas(result.cpass, fit_plot_pass_name, suffix);
//                         PrintCanvas(result.cfail, fit_plot_fail_name, suffix);
//                     }
//                 }
//             }
// 
//             // add histograms to hist container
//             hc.Add(h_num);
//             hc.Add(h_den);
//             hc.Add(h_eff);
//         }

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
