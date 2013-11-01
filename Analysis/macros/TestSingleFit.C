{
    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libTagAndProbeAnalysis.dylib");

    int pt  = 2;
    int eta = 3;
    rt::TH1Container hc("plots/SameSign/electron/SameSignDenBoth_SameSignNum/data_single_el.root");
    TH1* h_pass = hc[Form("h_pass_pt%d_vs_eta%d", pt, eta)];
    TH1* h_fail = hc[Form("h_fail_pt%d_vs_eta%d", pt, eta)];

    rt::TH1Container hc_mc("plots/SameSign/electron/SameSignDenBoth_SameSignNum/dy_full.root");
    TH1* h_pass_template = hc_mc[Form("h_pass_pt%d_vs_eta%d", pt, eta)];
    TH1* h_fail_template = hc_mc[Form("h_fail_pt%d_vs_eta%d", pt, eta)];

    result = tnp::PerformSimultaneousFit
    (
        tnp::Model::BreitWignerCB,
        tnp::Model::MCTemplate,
        tnp::Model::Exponential,
        tnp::Model::Exponential,
        h_pass,
        h_fail, 
        60.0, 
        120.0, 
        2.0, 
        "p_{T}", 
        "#eta", 
        h_pass_template, 
        h_fail_template
    );
}
