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
#include "TagAndProbe/Analysis/interface/LeptonSelections.h"
#include "TagAndProbe/Analysis/interface/Dataset.h"

using namespace std;


// hold the information relevant to the a dataset
// ------------------------------------------------------------------------------------ //

//struct Dataset
//{
//    // construct:
//    Dataset(const edm::ParameterSet& pset);
//
//    Dataset
//    (
//        const std::string& name,
//        const std::vector<std::string>& input_file_names,
//        const std::string& run_list = "",
//        const bool is_data = false
//    );
//
//    // members:
//    std::string m_name;
//    std::vector<std::string> m_input_file_names;
//    std::string m_run_list;
//    bool m_is_data;
//};
//
//Dataset::Dataset
//(
//    const std::string& name,
//    const std::vector<std::string>& input_file_names,
//    const std::string& run_list,
//    const bool is_data
//)
//    : m_name(name)
//    , m_input_file_names(input_file_names)
//    , m_run_list(run_list)
//    , m_is_data(is_data)
//{
//}
//
//Dataset::Dataset(const edm::ParameterSet& pset)
//    : m_name(pset.getParameter<std::string>("name"))
//    , m_input_file_names(pset.getParameter<std::vector<std::string> >("files"))
//    , m_run_list(pset.getParameter<std::string>("run_list"))
//    , m_is_data(pset.getParameter<bool>("is_data"))
//{
//}
//
//// non member functions
//std::vector<Dataset> GetDatasetsFromVPSet(const std::vector<edm::ParameterSet>& psets)
//{
//    std::vector<Dataset> results;
//    results.assign(psets.begin(), psets.end());
//    return results;
//}


// Looper class to hold all the variables and make the histograms 
// ------------------------------------------------------------------------------------ //

class MassPlotLooper
{
    public:
        // consstruct:
        MassPlotLooper
        (
            const std::string& output_file_name, 
            const tnp::Lepton::value_type lepton_type,
            const tnp::Selection::value_type numerator,
            const tnp::Selection::value_type denominator,
            const double mass_low,
            const double mass_high,
            const double mass_bin_width,
            const std::vector<double> pt_bins,
            const std::vector<double> eta_bins,
            const std::vector<double> phi_bins,
            const std::vector<double> nvtx_bins,
            const bool is_data,
            const std::string& pileup_hist_file, 
            const std::string& pileup_hist_name, 
            const std::string suffix,
            const bool verbose
        ); 

        // destroy:
        ~MassPlotLooper() {EndJob();}

        // book all the histograms
        void BookHists();

        // analyze the event
        int Analyze(long long entry);

        int operator()(long long entry);

        // end of the job tasks
        void EndJob();

    private:
        // analysis parameters
        std::string m_output_file_name;
        tnp::Lepton::value_type m_lepton_type;
        tnp::Selection::value_type m_num;
        tnp::Selection::value_type m_den;
        double m_mass_low;
        double m_mass_high;
        double m_mass_bin_width;
        std::vector<double> m_pt_bins;
        std::vector<double> m_eta_bins;
        std::vector<double> m_phi_bins;
        std::vector<double> m_nvtx_bins;
        bool m_is_data;
        std::string m_suffix;
        bool m_verbose;
        TH1D* h_pileup;

        // members
        rt::TH1Container m_hist_container;
};

MassPlotLooper::MassPlotLooper
(
    const std::string& output_file_name, 
    const tnp::Lepton::value_type lepton_type,
    const tnp::Selection::value_type numerator,
    const tnp::Selection::value_type denominator,
    const double mass_low,
    const double mass_high,
    const double mass_bin_width,
    const std::vector<double> pt_bins,
    const std::vector<double> eta_bins,
    const std::vector<double> phi_bins,
    const std::vector<double> nvtx_bins,
    const bool is_data,
    const std::string& pileup_hist_file, 
    const std::string& pileup_hist_name, 
    const std::string suffix,
    const bool verbose
)
    : m_output_file_name(output_file_name)
    , m_lepton_type(lepton_type)
    , m_num(numerator)
    , m_den(denominator)
    , m_mass_low(mass_low)
    , m_mass_high(mass_high)
    , m_mass_bin_width(mass_bin_width)
    , m_pt_bins(pt_bins)
    , m_eta_bins(eta_bins)
    , m_phi_bins(phi_bins)
    , m_nvtx_bins(nvtx_bins)
    , m_is_data(is_data)
    , m_suffix(suffix)
    , m_verbose(verbose)
    , h_pileup(rt::GetHistFromRootFile<TH1D>(pileup_hist_file, pileup_hist_name))
{
	BookHists();
}

void MassPlotLooper::BookHists()
{
    rt::TH1Container& hc = m_hist_container;

    // mass bins
    const int num_mass_bins = static_cast<int>(fabs(m_mass_high - m_mass_low)/m_mass_bin_width);

    // use |eta| ?
    const std::string phi_title = (not m_phi_bins.empty() and m_phi_bins.front() >= 0 ? "|#phi|" : "#phi");
    const std::string eta_title = (not m_eta_bins.empty() and m_eta_bins.front() >= 0 ? "|#eta|" : "#eta");

    // book pt histograms 
    if (m_pt_bins.size() >= 2)
    {
        for (size_t pt_bin = 0; pt_bin != m_pt_bins.size()-1; pt_bin++)
        {
            const std::string bin_title = Form("%1.0f GeV < p_{T} < %1.0f GeV", m_pt_bins[pt_bin], m_pt_bins[pt_bin+1]);
            hc.Add(new TH1F(Form("h_pass_pt%lu", pt_bin), Form("Passing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
            hc.Add(new TH1F(Form("h_fail_pt%lu", pt_bin), Form("Failing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
        }
    }

    // book eta histograms
    if (m_eta_bins.size() >= 2)
    {
        for (size_t eta_bin = 0; eta_bin != m_eta_bins.size()-1; eta_bin++)
        {
            const std::string bin_title = Form("%1.2f < %s < %1.2f", m_eta_bins[eta_bin], eta_title.c_str(), m_eta_bins[eta_bin+1]);
            hc.Add(new TH1F(Form("h_pass_eta%lu", eta_bin), Form("Passing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
            hc.Add(new TH1F(Form("h_fail_eta%lu", eta_bin), Form("Failing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
        }
    }

    // book phi histogram
    if (m_phi_bins.size() >= 2)
    {
        for (size_t phi_bin = 0; phi_bin != m_phi_bins.size()-1; phi_bin++)
        {
            const std::string bin_title = Form("%1.2f < %s < %1.2f", m_phi_bins[phi_bin], phi_title.c_str(), m_phi_bins[phi_bin+1]);
            hc.Add(new TH1F(Form("h_pass_phi%lu", phi_bin), Form("Passing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
            hc.Add(new TH1F(Form("h_fail_phi%lu", phi_bin), Form("Failing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
        }
    }

    // book # vertices histogram
    if (m_nvtx_bins.size() >= 2)
    {
        for (size_t nvtx_bin = 0; nvtx_bin != m_nvtx_bins.size()-1; nvtx_bin++)
        {
            const std::string bin_title = Form("%1.0f < # Vertices < %1.0f", m_nvtx_bins[nvtx_bin], m_nvtx_bins[nvtx_bin+1]);
            hc.Add(new TH1F(Form("h_pass_nvtx%lu", nvtx_bin), Form("Passing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
            hc.Add(new TH1F(Form("h_fail_nvtx%lu", nvtx_bin), Form("Failing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
        }
    }

    // book pt vs eta histograms 
    if (m_pt_bins.size() >= 2 and m_eta_bins.size() >= 2)
    {
        for (size_t pt_bin = 0; pt_bin != m_pt_bins.size()-1; pt_bin++)
        {
            for (size_t eta_bin = 0; eta_bin != m_eta_bins.size()-1; eta_bin++)
            {
                const std::string bin_title = Form("%1.0f GeV < p_{T} < %1.0f GeV, %1.2f < %s < %1.2f", m_pt_bins[pt_bin], m_pt_bins[pt_bin+1], m_eta_bins[eta_bin], eta_title.c_str(), m_eta_bins[eta_bin+1]);
                hc.Add(new TH1F(Form("h_pass_pt%lu_eta%lu", pt_bin, eta_bin), Form("Passing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
                hc.Add(new TH1F(Form("h_fail_pt%lu_eta%lu", pt_bin, eta_bin), Form("Failing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
            }
        }
    }

    // book eta vs phi histograms 
    if (m_eta_bins.size() >= 2 and m_phi_bins.size() >= 2)
    {
        for (size_t eta_bin = 0; eta_bin != m_eta_bins.size()-1; eta_bin++)
        {
            for (size_t phi_bin = 0; phi_bin != m_phi_bins.size()-1; phi_bin++)
            {
                const std::string bin_title = Form("%1.2f < %s < %1.2f, %1.2f < %s < %1.2f", m_eta_bins[eta_bin], eta_title.c_str(), m_eta_bins[eta_bin+1], m_phi_bins[phi_bin], phi_title.c_str(), m_phi_bins[phi_bin+1]);
                hc.Add(new TH1F(Form("h_pass_eta%lu_phi%lu", eta_bin, phi_bin), Form("Passing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
                hc.Add(new TH1F(Form("h_fail_eta%lu_phi%lu", eta_bin, phi_bin), Form("Failing probes (%s);tag & probe mass (GeV);Events / %1.1f (GeV)", bin_title.c_str(), m_mass_bin_width), num_mass_bins, m_mass_low, m_mass_high));
            }
        }
    }

    // TDRStyle
    rt::SetStyle();

    // sumw2
    hc.SetMarkerStyle(20);
    hc.SetMarkerSize(1.0);
    hc.SetLineWidth(1.0);
    hc.Sumw2();
}

int MassPlotLooper::Analyze(long long entry)
{
    // convenience
	rt::TH1Container& hc = m_hist_container;

    using namespace lepton_tree;

	try
    {
        // tag and probe
        // ------------------------------------------------------------------------------------ //

        if (m_verbose) {cout << Form("\nrun %u, ls %u, evt %u", run(), lumi(), event()) << endl;}

        // mode
        const bool is_mu   = (m_lepton_type == tnp::Lepton::Muon);
        const bool is_el   = (m_lepton_type == tnp::Lepton::Electron);
        const bool is_data = m_is_data;
        const bool is_mc   = (not is_data);

        // passes the probe denominator 
        if (not tnp::PassesSelection(m_lepton_type, m_den, is_data))
        {
            if (m_verbose) {cout << "Did not pass the denominator selection" << endl;}
            return 0;
        }

        unsigned int evt_sel = eventSelection();
        // Z/onia --> ee
        if (is_el)
        {
            if (((evt_sel & tnp::EventSelection::ZeeTagAndProbe   ) != tnp::EventSelection::ZeeTagAndProbe   ) && 
                ((evt_sel & tnp::EventSelection::OniaEETagAndProbe) != tnp::EventSelection::OniaEETagAndProbe))
            {
                if (m_verbose) {cout << "Did not pass Z/onia --> ee" << endl;}
                return 0;
            }
        }
        // Z/onia --> mm
        if (is_mu)
        {
            if (((evt_sel & tnp::EventSelection::ZmmTagAndProbe     ) != tnp::EventSelection::ZmmTagAndProbe     ) && 
                ((evt_sel & tnp::EventSelection::OniaMuMuTagAndProbe) != tnp::EventSelection::OniaMuMuTagAndProbe))
            {
                if (m_verbose) {cout << "Did not pass Z/onia --> mm" << endl;}
                return 0;
            }
        }

        // require OS leptons
        if((qProbe() * qTag()) > 0)
        {
            if (m_verbose) {cout << "Did not pass OS requirement" << endl;}
            return 0;
        }

        // require mass window around the Z
        const float mass = tagAndProbeMass();
        if (not (m_mass_low < mass && mass < m_mass_high))
        {
            if (m_verbose) {cout << "Did not pass Z mass requirement" << endl;}
            return 0;
        }

        // MC reqruied DeltaR(reco lepton, status 1 gen level lepton) < 0.2
        if (is_mc && gen_drs1() > 0.2)
        {
            if (m_verbose) {cout << "Did not DeltaR(lep, s1) < 0.2 requirement" << endl;}
            return 0;
        }

        // check pt boundaries
        const bool has_pt_bins = m_pt_bins.size() >= 2;
        const float probe_pt   = probe().pt();
        const float pt_min     = (m_pt_bins.empty() ? 999999.0 : m_pt_bins.front());
        const float pt_max     = (m_pt_bins.empty() ? 999999.0 : m_pt_bins.back() );
        if (has_pt_bins and not (pt_min < probe_pt && probe_pt < pt_max))
        {
            if (m_verbose) {cout << "outside pt bins" << endl;}
            return 0;
        }

        // check eta boundaries
        const bool has_eta_bins = m_eta_bins.size() >= 2;
        const float probe_eta   = fabs(is_el ? sceta() : probe().eta());
        const float eta_min     = (m_eta_bins.empty() ? 999999.0 : m_eta_bins.front());
        const float eta_max     = (m_eta_bins.empty() ? 999999.0 : m_eta_bins.back() );
        if (has_eta_bins and not (eta_min < probe_eta && probe_eta < eta_max))
        {
            if (m_verbose) {cout << "outside eta bins" << endl;}
            return 0;
        }

        // check phi boundaries
        const bool has_phi_bins = m_phi_bins.size() >= 2;
        const float probe_phi   = probe().phi();
        const float phi_min     = (m_phi_bins.empty() ? 999999.0 : m_phi_bins.front());
        const float phi_max     = (m_phi_bins.empty() ? 999999.0 : m_phi_bins.back() );
        if (has_phi_bins and not (phi_min < probe_phi && probe_phi < phi_max))
        {
            if (m_verbose) {cout << "outside phi bins" << endl;}
            return 0;
        }

        // check the nvtxs boundaries
        const bool has_nvtx_bins = m_nvtx_bins.size() >= 2;
        const float nvtxs        = nvtx();
        const float nvtx_min     = (m_nvtx_bins.empty() ? 999999.0 : m_nvtx_bins.front());
        const float nvtx_max     = (m_nvtx_bins.empty() ? 999999.0 : m_nvtx_bins.back() );
        if (has_nvtx_bins and not (nvtx_min < nvtxs && nvtxs < nvtx_max))
        {
            if (m_verbose) {cout << "outside nvtx bins" << endl;}
            return 0;
        }

        // find pt,eta,phi,nvtx bin
        unsigned int pt_bin   = (has_pt_bins   ? rt::find_bin(probe_pt , m_pt_bins  ) : -999999);
        unsigned int eta_bin  = (has_eta_bins  ? rt::find_bin(probe_eta, m_eta_bins ) : -999999);
        unsigned int phi_bin  = (has_phi_bins  ? rt::find_bin(probe_phi, m_phi_bins ) : -999999);
        unsigned int nvtx_bin = (has_nvtx_bins ? rt::find_bin(nvtxs    , m_nvtx_bins) : -999999);

        if (m_verbose and has_pt_bins  ) {cout << Form("pt %f  , pt_bin %u"  , probe_pt , pt_bin  ) << endl;}
        if (m_verbose and has_eta_bins ) {cout << Form("eta %f , eta_bin %u" , probe_eta, eta_bin ) << endl;}
        if (m_verbose and has_phi_bins ) {cout << Form("phi %f , phi %u"     , probe_phi, phi_bin ) << endl;}
        if (m_verbose and has_nvtx_bins) {cout << Form("nvtx %f, nvtx_bin %u", nvtxs    , nvtx_bin) << endl;}

        // PU re-weight
        //const float weight = 1.0;
        const float weight = (is_mc ? (scale1fb() * h_pileup->GetBinContent(nvtxs+1)) : 1.0);

        // Fill the histograms
        // ------------------------------------------------------------------------------------ //

        // passes the probe numerator 
        if (tnp::PassesSelection(m_lepton_type, m_num, is_data))
        {
            if (m_verbose) {cout << "passes the numerator selection" << endl;}

            // fill hists
            if (has_pt_bins                  ) {hc[Form("h_pass_pt%u", pt_bin)                 ]->Fill(mass, weight);}
            if (has_eta_bins                 ) {hc[Form("h_pass_eta%u", eta_bin)               ]->Fill(mass, weight);}
            if (has_phi_bins                 ) {hc[Form("h_pass_phi%u", phi_bin)               ]->Fill(mass, weight);}
            if (has_nvtx_bins                ) {hc[Form("h_pass_nvtx%u", nvtx_bin)             ]->Fill(mass, weight);}
            if (has_pt_bins and has_eta_bins ) {hc[Form("h_pass_pt%u_eta%u", pt_bin, eta_bin)  ]->Fill(mass, weight);}
            if (has_eta_bins and has_phi_bins) {hc[Form("h_pass_eta%u_phi%u", eta_bin, phi_bin)]->Fill(mass, weight);}
        }
        // fails the probe numerator 
        else
        {
            if (m_verbose) {cout << "fails the numerator selection" << endl;}

            // fill hists
            if (has_pt_bins                  ) {hc[Form("h_fail_pt%u", pt_bin)                 ]->Fill(mass, weight);}
            if (has_eta_bins                 ) {hc[Form("h_fail_eta%u", eta_bin)               ]->Fill(mass, weight);}
            if (has_phi_bins                 ) {hc[Form("h_fail_phi%u", phi_bin)               ]->Fill(mass, weight);}
            if (has_nvtx_bins                ) {hc[Form("h_fail_nvtx%u", nvtx_bin)             ]->Fill(mass, weight);}
            if (has_pt_bins and has_eta_bins ) {hc[Form("h_fail_pt%u_eta%u", pt_bin, eta_bin)  ]->Fill(mass, weight);}
            if (has_eta_bins and has_phi_bins) {hc[Form("h_fail_eta%u_phi%u", eta_bin, phi_bin)]->Fill(mass, weight);}
        }

        // done
        return 0;
    }
    catch (std::exception& e)
    {
        cout << Form("Fatal error on run %d, ls %d, event %d", run(), lumi(), event()) << endl;
        cout << e.what() << endl;
        cout << "Exiting..." << endl;
        exit(1);
    }

    // done
    return 0;
}

void MassPlotLooper::EndJob()
{
    // write output
    cout << "Writing histogram root file to: " << m_output_file_name << endl;
    m_hist_container.Write(m_output_file_name);
    if (not m_suffix.empty())
    {
        std::string output_print_path = lt::string_replace_all(m_output_file_name, ".root", "") + "/" + m_suffix;
        cout << "Printing histograms to: " << output_print_path  << endl;
        m_hist_container.Print(output_print_path, m_suffix);
    }
}

// call each loopers analysis function
int MassPlotLooper::operator()(long long entry)
{
    return Analyze(entry);
}

// wrapper to call multiple loopers on each event 
// ------------------------------------------------------------------------------------ //

struct MultiMassPlotLooper
{
    // construct:
    MultiMassPlotLooper(const std::vector<MassPlotLooper*>& loopers)
        : m_loopers(loopers)
    {}

    // call each loopers analysis function
    int operator()(long long entry)
    {
        for (size_t i = 0; i != m_loopers.size(); i++)
        {
            m_loopers.at(i)->Analyze(entry);
        }
        return 0;
    }

    // member
    std::vector<MassPlotLooper*> m_loopers;
};


// Peform a analysis given by Function on a chain
// ------------------------------------------------------------------------------------ //

template <typename Function>
int ScanChain
(
    Function analyze, 
    const tnp::Dataset& dataset,
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
                // Turned off since this doesn't work since there are multple lepton pairs per event...
                //DorkyEventIdentifier id = {run, evt, ls};
                //if (is_duplicate(id))
                //{
                //    if (verbose) {cout << "[ScanChain] good run file = " << run_list << endl;}
                //    duplicates++;
                //    continue;
                //}
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
        throw std::invalid_argument(Form("[tnp_make_plots] Error: ParametersSet 'process' is missing in your configuration file"));
    }

    // get the python configuration
    const edm::ParameterSet& process = edm::readPSetsFrom(pset_filename)->getParameter<edm::ParameterSet>("process");
    const edm::ParameterSet& tnp_cfg = process.getParameter<edm::ParameterSet>("tnp_make_plots");

    // Looper of parameter sets
    // Makes a set of histogram for each element of tnp_cfgs vector for each dataset
    // -------------------------------------------------------------------------------------------------//

    // get the inputs 
    const long long max_events                = tnp_cfg.getParameter<long long>("max_events");
    const tnp::Lepton::value_type lepton_type = tnp::GetLeptonFromString(tnp_cfg.getParameter<std::string>("lepton_type"));
    const double mass_low                     = tnp_cfg.getParameter<double>("mass_low" );
    const double mass_high                    = tnp_cfg.getParameter<double>("mass_high");
    const double mass_bin_width               = tnp_cfg.getParameter<double>("mass_bin_width");
    const bool verbose                        = tnp_cfg.getParameter<bool>("verbose");
    const std::string suffix                  = tnp_cfg.getParameter<std::string>("suffix");
    const std::string analysis_path           = lt::getenv("CMSSW_BASE") + "/src/TagAndProbe/Analysis";
    const std::string output_label            = tnp_cfg.getParameter<std::string>("output_label");
    const std::string pileup_hist_file        = tnp_cfg.getParameter<std::string>("pileup_hist_file");
    const std::string pileup_hist_name        = tnp_cfg.getParameter<std::string>("pileup_hist_name");
    const std::vector<double> pt_bins         = tnp_cfg.getParameter<std::vector<double> >("pt_bins");
    const std::vector<double> eta_bins        = tnp_cfg.getParameter<std::vector<double> >("eta_bins");
    const std::vector<double> phi_bins        = tnp_cfg.getParameter<std::vector<double> >("phi_bins");
    const std::vector<double> nvtx_bins       = tnp_cfg.getParameter<std::vector<double> >("nvtx_bins");
    const std::vector<tnp::Dataset> datasets  = tnp::GetDatasetsFromVPSet(tnp_cfg.getParameter<std::vector<edm::ParameterSet> >("datasets"));

    // numerator and denominator    
    const tnp::Selection::value_type num_selection = tnp::GetSelectionFromString(tnp_cfg.getParameter<std::string>("numerator"  ));
    const tnp::Selection::value_type den_selection = tnp::GetSelectionFromString(tnp_cfg.getParameter<std::string>("denominator"));

    // for each dataset makes the set of histograms
    // -------------------------------------------------------------------------------------------------//

    for (std::vector<tnp::Dataset>::const_iterator dataset_iter = datasets.begin(); dataset_iter != datasets.end(); dataset_iter++)
    {
        // convenience
        const tnp::Dataset& dataset = *dataset_iter;

        // for each dataset makes the set of histograms
        // -------------------------------------------------------------------------------------------------//


        // output ROOT file name
        // (i.e. analysis_path/plots/output_label/lepton_type/den_num/dataset.root)
        const std::string output_file_name = Form("%s/plots/%s/%s/%s_%s/%s.root",
            analysis_path.c_str(),
            output_label.c_str(),
            GetStringFromLepton(lepton_type).c_str(),
            GetStringFromSelection(den_selection).c_str(),
            GetStringFromSelection(num_selection).c_str(),
            dataset.m_name.c_str()
        );

        // print out the parameters for each run
        cout << "\n[tnp_make_plots] running with the following inputs:" << endl;
        printf("%-15s = %lld\n", "max_events"   , max_events                                   );
        printf("%-15s = %s\n"  , "name"         , dataset.m_name.c_str()                       );
        printf("%-15s = %s\n"  , "run_list"     , dataset.m_run_list.c_str()                   );
        printf("%-15s = %d\n"  , "is_data"      , dataset.m_is_data                            );
        printf("%-15s = %s\n"  , "num_selection", GetStringFromSelection(den_selection).c_str());
        printf("%-15s = %s\n"  , "den_selection", GetStringFromSelection(num_selection).c_str());

        // analysis looper
        ScanChain
        (
            MassPlotLooper
            (
                output_file_name,
                lepton_type,
                num_selection,
                den_selection,
                mass_low,
                mass_high,
                mass_bin_width,
                pt_bins,
                eta_bins,
                phi_bins,
                nvtx_bins,
                dataset.m_is_data,
                pileup_hist_file,
                pileup_hist_name,
                suffix,
                verbose
            ), 
            dataset,
            max_events
        );

    } // end loop over datasets 

    // done
    return 0;
}
catch (std::exception& e)
{
    cerr << "[tnp_make_plots] Error: failed..." << endl;
    cerr << e.what() << endl;
    return 1;
}

