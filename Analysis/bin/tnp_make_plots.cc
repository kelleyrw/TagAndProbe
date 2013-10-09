// C++ includes
#include <iostream>
#include <string>
#include <stdexcept>

// ROOT includes
#include "TH1F.h"
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TSystem.h"
#include "TChain.h"
#include "TBenchmark.h"
#include "TTreeCache.h"

// CMSSW includes
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// Tools
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"
#include "TagAndProbe/Analysis/interface/LeptonTree.h"
#include "TagAndProbe/Analysis/interface/DorkyEventIdentifier.h"
#include "TagAndProbe/Analysis/interface/goodrun.h"
#include "TagAndProbe/Analysis/interface/Measurement.h"

using namespace std;


// hold the information relevant to the a dataset
// ------------------------------------------------------------------------------------ //

struct Dataset
{
    // construct:
    Dataset(const edm::ParameterSet& pset);

    Dataset
    (
        const std::string& name,
        const std::vector<std::string>& input_file_names,
        const std::string& run_list = "",
        const bool is_data = false
    );

    // members:
    std::string m_name;
    std::vector<std::string> m_input_file_names;
    std::string m_run_list;
    bool m_is_data;
};

Dataset::Dataset
(
    const std::string& name,
    const std::vector<std::string>& input_file_names,
    const std::string& run_list,
    const bool is_data
)
    : m_name(name)
    , m_input_file_names(input_file_names)
    , m_run_list(run_list)
    , m_is_data(is_data)
{
}

Dataset::Dataset(const edm::ParameterSet& pset)
    : m_name(pset.getParameter<std::string>("name"))
    , m_input_file_names(pset.getParameter<std::vector<std::string> >("files"))
    , m_run_list(pset.getParameter<std::string>("run_list"))
    , m_is_data(pset.getParameter<bool>("is_data"))
{
}

// non member functions
std::vector<Dataset> GetDatasetsFromVPSet(const std::vector<edm::ParameterSet>& psets)
{
    std::vector<Dataset> results;
    results.assign(psets.begin(), psets.end());
    return results;
}


// Looper class to hold all the variables and make the histograms 
// ------------------------------------------------------------------------------------ //

class LeptonTreeLooper
{
    public:
        // consstruct:
        LeptonTreeLooper
        (
            const std::string& output_file_name, 
            const tnp::Lepton::value_type lepton_type,
            const tnp::Selection::value_type numerator,
            const tnp::Selection::value_type denominator,
            const float fit_mass_low,
            const float fit_mass_high,
            const bool verbose
        ); 
        ~LeptonTreeLooper() {EndJob();}
        void BookHists();
//         int Analyze(long event);
        void EndJob();

        // operator wrapper to act as a function
        int operator()(long event);

    private:
        // analysis parameters
        std::string m_output_file_name;
        tnp::Lepton::value_type m_lepton_type;
        tnp::Selection::value_type m_num;
        tnp::Selection::value_type m_den;
        float m_mass_low;
        float m_mass_high;
        TH1D* h_pu;
        bool m_verbose;

        // members
        rt::TH1Container m_hist_container;
};

LeptonTreeLooper::LeptonTreeLooper
(
    const std::string& output_file_name, 
    const tnp::Lepton::value_type lepton_type,
    const tnp::Selection::value_type numerator,
    const tnp::Selection::value_type denominator,
    const float fit_mass_low,
    const float fit_mass_high,
    const bool verbose
)
    : m_output_file_name(output_file_name)
    , m_num(numerator)
    , m_den(denominator)
    , m_mass_low(fit_mass_low)
    , m_mass_high(fit_mass_high)
//     , h_pileup(rt::GetHistFromRootFile<TH1D>("data/puWeights_Summer12_53x_Observed.root", "puWeights"))
    , m_verbose(verbose)
{
	BookHists();
}

void LeptonTreeLooper::BookHists()
{
	rt::TH1Container& hc = m_hist_container;

    hc.Add(new TH1F("h_test", "h_test", 30, 60, 120));

//    // mass bins
//    const int nmass_bins = static_cast<int>(fabs(m_mhigh - m_mlow)/tp::MassBinWidth);
//
//    // book pt hists 
//    for (size_t ptbin = 0; ptbin != npt_bins; ptbin++)
//    {
//        for (size_t etabin = 0; etabin != neta_bins; etabin++)
//        {
//            const std::string bin_title = Form("%1.0f GeV < p_{T} < %1.0f GeV, %1.2f < |#eta| < %1.2f", pt_bins[ptbin], pt_bins[ptbin+1], eta_bins[etabin], eta_bins[etabin+1]);
//
//            hc.Add(new TH1F(Form("h_pass_pt%lu_eta%lu", ptbin, etabin), Form("Passing probes (%s); tag & probe mass (GeV); Events / %1.1f (GeV)", bin_title.c_str(), tp::MassBinWidth), nmass_bins, m_mlow, m_mhigh));
//            hc.Add(new TH1F(Form("h_fail_pt%lu_eta%lu", ptbin, etabin), Form("Failing probes (%s); tag & probe mass (GeV); Events / %1.1f (GeV)", bin_title.c_str(), tp::MassBinWidth), nmass_bins, m_mlow, m_mhigh));
//        }
//    }

    // sumw2
    hc.SetMarkerStyle(20);
    hc.SetMarkerSize(1.0);
    hc.Sumw2();
}

int LeptonTreeLooper::operator()(long event)
{
    using namespace lepton_tree;
    return 0;
}

void LeptonTreeLooper::EndJob()
{
    // write output
    cout << "Writing histogram root file to: " << m_output_file_name << endl;
    m_hist_container.Write(m_output_file_name);
//     if (m_print)
//     {
//         std::string output_print_path = rt::dirname(m_root_file->GetName()) + "/" + m_suffix;
//         m_hist_container.Print(output_print_path, m_suffix);
//     }
}

// wrapper to call multiple loopers on each event 
// ------------------------------------------------------------------------------------ //

struct MultiLeptonTreeLooper
{
    // construct:
    MultiLeptonTreeLooper(const std::vector<LeptonTreeLooper*>& loopers)
        : m_loopers(loopers)
    {}

    // call each loopers analysis function
    int operator()(long event)
    {
        for (size_t i = 0; i != m_loopers.size(); i++)
        {
            m_loopers.at(i)->operator()(event);
        }
        return 0;
    }

    // member
    std::vector<LeptonTreeLooper*> m_loopers;
};


// Peform a analysis given by Function on a chain
// ------------------------------------------------------------------------------------ //

template <typename Function>
int ScanChain
(
    Function analyze, 
    const Dataset& dataset,
    const int num_events = -1,
    const bool verbose = false,
    const int evt_run = -1,
    const int evt_lumi = -1,
    const int evt_event = -1
)
{
    using namespace std;

    TChain* chain = new TChain("leptons");
    for (size_t i = 0; i != dataset.m_input_file_names.size(); i++)
    {
        chain = rt::MakeTChain(dataset.m_input_file_names.at(i), "leptonTree", chain, verbose);
        rt::PrintFilesFromTChain(chain);
    }

    // test chain
    if (!chain)
    {
        throw std::invalid_argument("[ScanChain] Error: chain is NULL!");
    }
    if (chain->GetListOfFiles()->GetEntries()<1)
    {
        throw std::invalid_argument("[ScanChain] Error: chain has no files!");
    }
    if (not chain->GetFile())
    {
        throw std::invalid_argument("[ScanChain] Error: chain has no files or file path is invalid!");
    }

    // tree name
    string tree_name = chain->GetName();

    // set the "good run" list 
    const std::string& run_list = dataset.m_run_list;
    if (!run_list.empty())
    {
        set_goodrun_file(run_list.c_str());
    }

    // reset duplicate counter
    reset();

    // benchmark
    TBenchmark bmark;
    bmark.Start("benchmark");

    // events counts and max events
    int i_permilleOld = 0;
    long num_events_total = 0;
    long num_events_chain = (num_events >= 0 && num_events < chain->GetEntries()) ? num_events : chain->GetEntries();
    TObjArray* list_of_files = chain->GetListOfFiles();
    TIter file_iter(list_of_files);
    TFile* current_file = NULL;

    // count the duplicates and bad events
    unsigned long duplicates = 0;
    unsigned long bad_events = 0;

    // loop over files in the chain
    while ((current_file = (TFile*)file_iter.Next()))
    {
        TFile *file = TFile::Open(current_file->GetTitle());
        if (!file || file->IsZombie())
        {
            throw std::runtime_error(Form("[ScanChain] Error: File from TChain is invalid or corrupt: %s", current_file->GetTitle()));
        }

        // get the trees in each file
        TTree *tree = dynamic_cast<TTree*>(file->Get(tree_name.c_str()));
        if (!tree || tree->IsZombie())
        {
            throw std::runtime_error(Form("[ScanChain] Error: File from TChain has an invalid TTree or is corrupt: %s", current_file->GetTitle()));
        }
        TTreeCache::SetLearnEntries(10);
        tree->SetCacheSize(128*1024*1024);
        lepton_tree_obj.Init(tree);

        // Loop over Events in current file
        if (num_events_total >= num_events_chain) continue;
        long num_events_tree = tree->GetEntriesFast();

        // loop over events to Analyze
        for (long event = 0; event != num_events_tree; ++event)
        {
            // quit if the total is > the number in the chain
            if (num_events_total >= num_events_chain) continue;

            // load the entry
            tree->LoadTree(event);
            lepton_tree_obj.GetEntry(event);
            ++num_events_total;

            // pogress
            int i_permille = (int)floor(1000 * num_events_total / float(num_events_chain));
            if (i_permille != i_permilleOld)
            {
                printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                fflush(stdout);
                i_permilleOld = i_permille;
            }

            // check run/ls/evt
            unsigned int run = lepton_tree_obj.run();
            unsigned int ls  = lepton_tree_obj.lumi();
            unsigned int evt = lepton_tree_obj.event();

            if (evt_event >= 0)
            {
                if (evt==static_cast<unsigned int>(evt_event))
                {
                    if (verbose) {cout << "[ScanChain] selected event:\t" << evt << endl;}
                }
                else
                {
                    continue;
                }
            }
            if (evt_lumi >= 0)
            {
                if (ls==static_cast<unsigned int>(evt_lumi))
                {
                    if (verbose) {cout << "[ScanChain] selected lumi:\t" << ls << endl;}
                }
                else
                {
                    continue;
                }
            }
            if (evt_run >= 0)
            {
                if (ls==static_cast<unsigned int>(evt_run))
                {
                    if (verbose) {cout << "[ScanChain] selected run:\t" << run << endl;}
                }
                else
                {
                    continue;
                }
            }

            // filter out events
            if (dataset.m_is_data)
            {
                if (!run_list.empty())
                {
                    // check for good run and events
                    if(!goodrun(run, ls)) 
                    {
                        if (verbose) {cout << "[ScanChain] Bad run and lumi:\t" << run << ", " << ls << endl;}
                        bad_events++;
                        continue;
                    }
                }

                // check for dupiclate run and events
                DorkyEventIdentifier id = {run, evt, ls};
                if (is_duplicate(id))
                {
                    if (verbose) {cout << "[ScanChain] good run file = " << run_list << endl;}
                    duplicates++;
                    continue;
                }
            }

            // print run/ls/event
            if (verbose)
            {
                cout << Form("[ScanChain] run %d, ls %d, evt %d", run, ls, evt) << endl;
            }

            // analysis
            analyze(event);

        } // end event loop

        // close current file
        file->Close();
        delete file;

    } // end file loop

    // print warning if the totals don't line up
    if (num_events_chain != num_events_total) 
    {
        cout << "[ScanChain] Error: number of events from the files " 
            << "(" << num_events_chain << ") " 
            << "is not equal to the total number of events "
            << "(" << num_events_total << ")." 
            << endl;
    }

    // the benchmark results 
    // -------------------------------------------------------------------------------------------------//
    bmark.Stop("benchmark");
    cout << endl;
    cout << num_events_total << " Events Processed" << endl;
    cout << "# of bad events filtered = " << bad_events << endl; 
    cout << "# of duplicates filtered = " << duplicates << endl; 
    cout << "------------------------------" << endl;
    cout << "CPU  Time: " << Form("%.01f", bmark.GetCpuTime("benchmark" )) << endl;
    cout << "Real Time: " << Form("%.01f", bmark.GetRealTime("benchmark")) << endl;
    cout << endl;

    // done
    return 0;
}

int main(int argc, char **argv)
try
{
    // parse the inputs
    // -------------------------------------------------------------------------------------------------//

    gSystem->Load( "libFWCoreFWLite" );
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
        throw std::invalid_argument(Form("[tnp_make_plots] Error: ParametersSet 'process' is missing in your configuration file"));
    }

    // get the python configuration
    const edm::ParameterSet& process               = edm::readPSetsFrom(pset_filename)->getParameter<edm::ParameterSet>("process");
    const std::vector<edm::ParameterSet>& tnp_cfgs = process.getParameter<std::vector<edm::ParameterSet>>("tnp_make_plots");

    // Looper of parameter sets
    // Makes a set of histogram for each element of tnp_cfgs vector for each dataset
    // -------------------------------------------------------------------------------------------------//

    for (std::vector<edm::ParameterSet>::const_iterator tnp_cfg_iter = tnp_cfgs.begin(); tnp_cfg_iter != tnp_cfgs.end(); tnp_cfg_iter++)
    {
        // convenience
        const edm::ParameterSet& tnp_cfg = *tnp_cfg_iter;

        // get the inputs 
        const long long max_events                = tnp_cfg.getParameter<long long>("max_events");
        const tnp::Lepton::value_type lepton_type = tnp::GetLeptonFromString(tnp_cfg.getParameter<std::string>("lepton_type"));
        const float mass_low                      = tnp_cfg.getParameter<double>("mass_low" );
        const float mass_high                     = tnp_cfg.getParameter<double>("mass_high");
        const bool verbose                        = tnp_cfg.getParameter<bool>("verbose");
        const std::string analysis_path           = lt::getenv("TNP");
        const std::string output_label            = tnp_cfg.getParameter<string>("output_label");


        const std::vector<Dataset> datasets = GetDatasetsFromVPSet(tnp_cfg.getParameter<std::vector<edm::ParameterSet> >("datasets"));

        // for each dataset makes the set of histograms
        // -------------------------------------------------------------------------------------------------//

        for (std::vector<Dataset>::const_iterator dataset_iter = datasets.begin(); dataset_iter != datasets.end(); dataset_iter++)
        {
            // convenience
            const Dataset& dataset = *dataset_iter;

            // for each dataset makes the set of histograms
            // -------------------------------------------------------------------------------------------------//

            std::vector<std::string> num_sel_strings = tnp_cfg.getParameter<std::vector<std::string> >("numerator"  );
            std::vector<std::string> den_sel_strings = tnp_cfg.getParameter<std::vector<std::string> >("denominator");

            // numerator and denominator selections should match up one-to-one
            // and thus be the same size vectors
            if (num_sel_strings.size() != den_sel_strings.size())
            {
                throw std::invalid_argument(Form("[tnp_make_plots] Error: ParametersSet value of 'numerator' and 'denominator' must be the same size"));
            }

            // build a vector of loopers
            std::vector<LeptonTreeLooper*> tnp_loopers;
            for (size_t i = 0; i != num_sel_strings.size(); i++)
            {
                // numerator and denominator    
                const tnp::Selection::value_type num_selection = tnp::GetSelectionFromString(num_sel_strings.at(i));
                const tnp::Selection::value_type den_selection = tnp::GetSelectionFromString(den_sel_strings.at(i));

                // output ROOT file name
                // (i.e. analysis_path/plots/output_label/lepton_type/den_num/dataset.root)
                const std::string output_file_name = Form("%s/plots/%s/%s/%s_%s/%s.root",
                    analysis_path.c_str(),
                    output_label.c_str(),
                    GetStringFromLepton(lepton_type).c_str(),
                    den_sel_strings.at(i).c_str(),
                    num_sel_strings.at(i).c_str(),
                    dataset.m_name.c_str()
                );

                // make the plots
                LeptonTreeLooper* tnp_looper = new LeptonTreeLooper 
                (
                     output_file_name,
                     lepton_type,
                     num_selection,
                     den_selection,
                     mass_low,
                     mass_high,
                     verbose
                ); 
                tnp_loopers.push_back(tnp_looper);
            }

            // print out the parameters for each run
            cout << "\n[tnp_make_plots] running with the following inputs:" << endl;
            printf("%-15s = %lld\n", "max_events"   , max_events                              );
            printf("%-15s = %s\n"  , "name"         , dataset.m_name.c_str()                  );
            printf("%-15s = %s\n"  , "run_list"     , dataset.m_run_list.c_str()              );
            printf("%-15s = %d\n"  , "is_data"      , dataset.m_is_data                       );
            printf("%-15s = {%s}\n", "num_selection", lt::string_join(num_sel_strings).c_str());
            printf("%-15s = {%s}\n", "den_selection", lt::string_join(den_sel_strings).c_str());

            // scan the chain and actually make the plots
            ScanChain(MultiLeptonTreeLooper(tnp_loopers), dataset, max_events);

            // cleanup
            lt::delete_container(tnp_loopers);
        }

    } // end loop over tnp_cfgs

    // done
    return 0;
}
catch (std::exception& e)
{
    cerr << "[tnp_make_plots] Error: failed..." << endl;
    cerr << e.what() << endl;
    return 1;
}

