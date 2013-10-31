#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"

void CreateMCPileupHist(const std::string& input_chain_path, const std::string& output_file_name, const long num_events = 1000000)
{
    // get the chain
    TChain e("Events");
    e.Add(input_chain_path.c_str());

    // fill the hist
    TH1D* h_pileup = new TH1D("hNPUTrue", "# True PU interactions;# interactions", 2000, 0, 200);
    h_pileup->Sumw2();
    h_pileup->SetLineColor(kBlue);
    h_pileup->SetLineWidth(2);
    h_pileup->SetMarkerStyle(20);
    h_pileup->SetMarkerColor(kBlack);
    h_pileup->SetMarkerSize(0.2);
    if (num_events >= 0)
    {
        e.Draw("puInfo_trueNumInteractions[1]>>hNPUTrue", "", "goff", num_events);
    }
    else
    {
        e.Draw("puInfo_trueNumInteractions[1]>>hNPUTrue", "", "goff");
    }

    // output
    TFile output_file(output_file_name.c_str(), "RECREATE");
    h_pileup->Write();
    output_file.Close();
}
