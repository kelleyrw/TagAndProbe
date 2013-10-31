#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <stdexcept>
#include <iostream>

void ComputePUWeightsTrue
(
    const std::string& source_file_name = "data/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_NPUTrue.root",
	const std::string& target_file_name = "data/MyDataPileupHistogramTrue_8TeV.root",
	const std::string& output_file_name = "data/puWeights_Summer12_53x_True.root"
) 
{
    try
    {
        // PU distribution from data
        const int re_bin = 10;
        TFile target_file(target_file_name.c_str());
        if (target_file.IsZombie())
        {
            throw std::runtime_error("[ComputePUWeightsTrue] Error: target file not found!");
        }
        TH1D* const h_target = dynamic_cast<TH1D*>(target_file.Get("pileup"));
        if (not h_target)
        {
            throw std::runtime_error("[ComputePUWeightsTrue] Error: target histogram not found!");
        }
        //h_target->Sumw2();
        h_target->Rebin(re_bin);
        h_target->Scale(1.0/h_target->Integral());

        // PU distribution from data MC
        TFile source_file(source_file_name.c_str());
        if (source_file.IsZombie())
        {
            throw std::runtime_error("[ComputePUWeightsTrue] Error: source file not found!");
        }
        TH1D* const h_den  = dynamic_cast<TH1D*>(source_file.Get("hNPUTrue"));
        if (not h_den)
        {
            throw std::runtime_error("[ComputePUWeightsTrue] Error: denominator histogram not found!");
        }
        //h_den->Sumw2();
        h_den->Rebin(re_bin);
        h_den->Scale(1.0/h_den->Integral());

        // create the PU distribtuions
        TH1D h_pu_weights("puWeights", "puWeights", 2000/re_bin, 0.0, 200.0);
        h_pu_weights.Sumw2();
        h_pu_weights.Divide(h_target, h_den);

        // save the output
        TFile output_file(output_file_name.c_str(), "RECREATE");
        h_pu_weights.Write();

        // cleanup
        output_file.Close();
        target_file.Close();
        source_file.Close();

        // done
        return;
    }
    catch (std::exception& e)
    {
        std::cerr << "[ComputePUWeightsTrue] failed..." << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }

}
