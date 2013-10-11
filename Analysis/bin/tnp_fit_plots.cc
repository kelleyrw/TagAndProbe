// C++ includes
#include <iostream>
#include <string>
#include <stdexcept>

// BOOST
#include "boost/multi_array.hpp"

// ROOT includes
#include "TH1F.h"
// #include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
// #include "TSystem.h"
// #include "TChain.h"
// #include "TBenchmark.h"
// #include "TTreeCache.h"

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// Tools
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "TagAndProbe/Analysis/interface/PerformFits.h"
// #include "TagAndProbe/Analysis/interface/LeptonTree.h"
// #include "TagAndProbe/Analysis/interface/DorkyEventIdentifier.h"
// #include "TagAndProbe/Analysis/interface/goodrun.h"
// #include "TagAndProbe/Analysis/interface/Measurement.h"
// #include "TagAndProbe/Analysis/interface/LeptonSelections.h"

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
        throw std::invalid_argument(Form("[tnp_fit_plots] Error: # of bins (%u) do not add up with the # of models (%lu) in the configuration", num_bins*num_categories, model_strings.size()));
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
        throw std::invalid_argument(Form("[tnp_fit_plots] Error: # of bins (%u) do not add up with the # of models (%lu) in the configuration", (num_a_bins)*(num_b_bins)*num_categories, model_strings.size()));
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



// The main program 
// ------------------------------------------------------------------------------------ //

int main(int argc, char **argv)
try
{
    // parse the inputs
    // -------------------------------------------------------------------------------------------------//

    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();

    // check that the python is passed
    if (argc < 2)
    {
        throw std::invalid_argument(Form("Usage : %s [parameters.py]", argv[0]));
    }

    // check that 
    const std::string pset_filename = argv[1];
    if (!edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process"))
    {
        throw std::invalid_argument(Form("[tnp_fit_plots] Error: ParametersSet 'process' is missing in your configuration file"));
    }

    // get the python configuration
    const edm::ParameterSet& process = edm::readPSetsFrom(pset_filename)->getParameter<edm::ParameterSet>("process");
    const edm::ParameterSet& tnp_cfg = process.getParameter<edm::ParameterSet>("tnp_fit_plots");

    const std::vector<double> pt_bins   = tnp_cfg.getParameter<std::vector<double> >("pt_bins");
    const std::vector<double> eta_bins  = tnp_cfg.getParameter<std::vector<double> >("eta_bins");
    const std::vector<double> phi_bins  = tnp_cfg.getParameter<std::vector<double> >("phi_bins");
    const std::vector<double> nvtx_bins = tnp_cfg.getParameter<std::vector<double> >("nvtx_bins");
    const unsigned int npt_bins         = pt_bins.size()-1;
    const unsigned int neta_bins        = eta_bins.size()-1;
    const unsigned int nphi_bins        = phi_bins.size()-1;
    const unsigned int nnvtx_bins       = nvtx_bins.size()-1;

    // get the bin models pt vs eta
    const std::vector<std::string> pt_vs_eta_model_strings = tnp_cfg.getParameter<std::vector<std::string> >("pt_vs_eta_models");
    const std::vector<std::string> pt_model_strings        = tnp_cfg.getParameter<std::vector<std::string> >("pt_models");
    if ((npt_bins)*(neta_bins)*num_categories != pt_vs_eta_model_strings.size())
    {
        throw invalid_argument("[tnp_fit_plots] Error: pt vs eta bins do not line up with pt and eta models in configuation");
    }

    cout << "num strings = " << pt_vs_eta_model_strings.size() << endl;
//     cout << "index\t" << "pt_bin" << "\t" << "eta_bin" << "\t" << "model_bin" << " = " << "pt_vs_eta_model_strings.at(i)" << endl;
    for (size_t pt_bin = 0; pt_bin != npt_bins; pt_bin++)
    {
        for (size_t eta_bin = 0; eta_bin != neta_bins; eta_bin++)
        {
            for (size_t model_bin = 0; model_bin != num_categories; model_bin++)
            {
                int index = pt_bin*(neta_bins * num_categories) + eta_bin*(num_categories) + model_bin;
//                 cout << index << "\t" <<  pt_bin << "\t" << eta_bin << "\t" << model_bin << " = " << pt_vs_eta_model_strings.at(index) << endl;
            }
        }
    }

//     typedef Array3D::index index;
    ModelArray3D A(boost::extents[npt_bins][neta_bins][num_categories]);

    for (size_t pt_bin = 0; pt_bin != npt_bins; pt_bin++)
    {
        for (size_t eta_bin = 0; eta_bin != neta_bins; eta_bin++)
        {
            for (size_t model_bin = 0; model_bin != num_categories; model_bin++)
            {

                int index = pt_bin*(neta_bins * num_categories) + eta_bin*(num_categories) + model_bin;
                A[pt_bin][eta_bin][model_bin] = tnp::GetModelFromString(pt_vs_eta_model_strings.at(index));
//                 cout << index << "\t" <<  pt_bin << "\t" << eta_bin << "\t" << model_bin << " = " << pt_vs_eta_model_strings.at(index) << endl;
            }
        }
    }

    ModelArray3D pt_vs_eta_models = GetModelArrayFromVString(pt_vs_eta_model_strings, pt_bins, eta_bins);
    for (size_t pt_bin = 0; pt_bin != npt_bins; pt_bin++)
    {
        for (size_t eta_bin = 0; eta_bin != neta_bins; eta_bin++)
        {
            for (size_t model_bin = 0; model_bin != num_categories; model_bin++)
            {

                int index = pt_bin*(neta_bins * num_categories) + eta_bin*(num_categories) + model_bin;
                tnp::Model::value_type model = pt_vs_eta_models[pt_bin][eta_bin][model_bin];
                cout << index << "\t" <<  pt_bin << "\t" << eta_bin << "\t" << model_bin << " = " << pt_vs_eta_model_strings.at(index) << "\t" << GetStringFromModel(model) << endl;
            }
        }
    }
    ModelArray2D pt_models = GetModelArrayFromVString(pt_model_strings, pt_bins);
    for (size_t pt_bin = 0; pt_bin != npt_bins; pt_bin++)
    {
        for (size_t model_bin = 0; model_bin != num_categories; model_bin++)
        {

            int index = pt_bin*(num_categories) + model_bin;
            tnp::Model::value_type model = pt_models[pt_bin][model_bin];
            cout << index << "\t" <<  pt_bin << "\t" << model_bin << " = " << pt_model_strings.at(index) << "\t" << GetStringFromModel(model) << endl;
        }
    }

    // done
    return 0;
}
catch (std::exception& e)
{
    cerr << "[tnp_fit_plots] Error: failed..." << endl;
    cerr << e.what() << endl;
    return 1;
}
