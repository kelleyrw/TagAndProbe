#include <iostream>
#include <string>
#include "boost/program_options.hpp"

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

// #include "DataFormats/FWLite/interface/Event.h"
// #include "DataFormats/Common/interface/Handle.h"
// #include "FWCore/FWLite/interface/AutoLibraryLoader.h"

// #include "DataFormats/FWLite/interface/InputSource.h"
// #include "DataFormats/FWLite/interface/OutputFiles.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// #include "DataFormats/MuonReco/interface/Muon.h"
// #include "DataFormats/PatCandidates/interface/Muon.h"
// #include "PhysicsTools/FWLite/interface/TFileService.h"

using namespace std;

int main(int argc, char* argv[])
try
{
    // inputs
    // -------------------------------------------------------------------------------------------------//

    std::string test_string = "[tnp] testing";

    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help"          , "print this menu")
        ("ts"            , po::value<std::string>(&test_string)->required(), "REQUIRED: test_string")
        ;

    // parse it
    try
    {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) 
        {
            cout << desc << "\n";
            return 1;
        }

        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        cerr << e.what() << "\nexiting" << endl;
        cout << desc << "\n";
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown error!" << "\n";
        return false;
    }

    // print inputs
    cout << "inputs:" << endl;
    cout << "test_string :\t" << test_string << endl;

    // test python
    // -------------------------------------------------------------------------------------------------//
    const std::string pset_filename = "config/tnp_test.py";

    // test if it exists
    if( !edm::readPSetsFrom(pset_filename)->existsAs<edm::ParameterSet>("process") ){
        std::cout << " ERROR: ParametersSet 'process' is missing in your configuration file" << std::endl; exit(0);
    }

    // get the python configuration
    const edm::ParameterSet& process = edm::readPSetsFrom(pset_filename)->getParameter<edm::ParameterSet>("process");
    const edm::ParameterSet& cfg     = process.getParameter<edm::ParameterSet>("TagAndProbeAnalysis");
    const std::string baby_file      = cfg.getParameter<std::string>("baby_file");
    const std::string plot_file      = cfg.getParameter<std::string>("plot_file");
    const int max_events             = cfg.getParameter<int>("maxEvents");

    // get the python configuration
    cout << "baby ntuple file name = " << baby_file << endl;
    cout << "output plot file name = " << plot_file << endl;
    cout << "max events            = " << max_events << endl;

    // done
    return 0;
}
catch (std::exception& e)
{
    cerr << "[ss2012_anlaysis] Error: failed..." << endl;
    cerr << e.what() << endl;
    return 1;
}

