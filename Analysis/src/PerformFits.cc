#include "TagAndProbe/Analysis/interface/PerformFits.h"

// c++
#include <string>
#include <stdexcept>
#include <algorithm>

// ROOT
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPaveText.h"

// RooFit
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooExtendPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooWorkspace.h"

// CMSSW
#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"

// Tools
#include "AnalysisTools/RootTools/interface/RootTools.h"
#include "AnalysisTools/LanguageTools/interface/LanguageTools.h"

using namespace std;

typedef RooAbsPdf* RooAbsPdfPtr;

// physical constants
// from: http://pdg.lbl.gov/
const double Mz     = 91.1876;
const double GammaZ = 2.4952;

namespace tnp
{
    // ----------------------------------------------------------------- //
    // PDF defintions
    // ----------------------------------------------------------------- //

    // wrapper to get the PDFs
    // ----------------------------------------------------------------- //

    void AddModelToWorkspace(Model::value_type model, RooWorkspace &w, const std::string& model_name)
    {
        // convenience
        const std::string unique_label = lt::string_replace_first(model_name, "model_", "_").c_str(); 
        char const * const ul = unique_label.c_str(); 

        switch (model)
        {
            case Model::BreitWignerCB: 
            {
                // breitwigner
                w.factory(Form("BreitWigner::bw%s(mass,mz%s[%f,80,100],gammaz%s[%f,0.1,10])", ul, ul, Mz, ul, GammaZ));

                // crystal ball
                w.factory(Form("RooCBShape::cb%s(mass,mean%s[0,-10,10],sigma%s[1,0.1,10],alpha%s[5,0,20],n%s[1,0,10])", ul, ul, ul, ul, ul));

                // convolution
                w.factory(Form("FCONV::%s(mass,bw%s,cb%s)", model_name.c_str(), ul, ul));
                break;
            }
            case Model::MCTemplate:
            {
                // gaussian
                w.factory(Form("Gaussian::gaus%s(mass,mean%s[0,-10,10],sigma%s[2,0,10])", ul, ul, ul));

                // template histogram
                TH1* const h_template = dynamic_cast<TH1*>(w.obj(Form("h_template%s", ul)));
                if (h_template == NULL)
                {
                    throw std::invalid_argument("[tnp::AddModelToWorkspace] Error: template histogram is NULL!");
                }

                // mc template
                RooDataHist data_hist(Form("data_hist%s", ul), Form("data_hist%s", ul), RooArgSet(*w.var("mass")), h_template);
                w.import(data_hist);
                w.factory(Form("HistPdf::hist_pdf%s(mass, data_hist%s, 1)", ul, ul));

                // convolution
                w.factory(Form("FCONV::%s(mass,hist_pdf%s,gaus%s)", model_name.c_str(), ul, ul));
                break;
            }
            case Model::Exponential: 
            {
                w.factory(Form("Exponential::%s(mass,t%s[-0.1,-1.0,0.0])", model_name.c_str(), ul)); 
                break;
            }
            case Model::Argus:
            {
                w.factory(Form("ArgusBG::%s(mass,m%s[-20],c%s[-100],p%s[-1])",  model_name.c_str(), ul, ul, ul));
                break;
            }
            case Model::ErfExp:
            {
                w.factory(Form("CMSShape::%s(mass,alfa%s[50,5,200],beta%s[0.01,0,10],gamma%s[0.1,0,1.0],peak%s[%f,85,97])", model_name.c_str(), ul, ul, ul, ul, Mz));
                w.var(Form("peak%s",ul))->setConstant(kTRUE);  
                break;
            }
            case Model::Chebychev:
            {
                w.factory(Form("Chebychev::%s(mass,{a0%s[0,-10,10],a1%s[0,-10,10]})", model_name.c_str(), ul, ul));
                break;
            }
            case Model::ChebyExp:
            {
                // chebychev polynomial
                w.factory(Form("Chebychev::chebychev%s(mass,{a0%s[0,-10,10],a1%s[0,-10,10]})", ul, ul, ul));

                // exponential
                w.factory(Form("Exponential::exp%s(mass,t%s[-0.1,-1.0,0.0])", ul, ul)); 

                // convolution
                w.factory(Form("FCONV::%s(mass,chebychev%s,exp%s)", model_name.c_str(), ul, ul));
                break;
            }
            case Model::Linear:
            {
                w.factory(Form("Polynomial::%s(mass,{a0%s[0],a1%s[0]})", model_name.c_str(), ul, ul));
                break;
            }
            case Model::Poly2:
            {
                w.factory(Form("Polynomial::%s(mass,{a0%s[-0.02,0,-0.05],a1%s[0.000, -0.01, 0.01]},a2%s[0]})", model_name.c_str(), ul, ul, ul));
                break;
            }
            case Model::Poly3:
            {
                w.factory(Form("Polynomial::%s(mass,{a0%s[0],a1%s[0],a2%s[0],a3%s[0]})", model_name.c_str(), ul, ul, ul, ul));
                break;
            }
            case Model::Poly6:
            {
                w.factory(Form("Polynomial::%s(mass,{a0%s[0],a1%s[0],a2%s[0],a3%s[0],a4%s[0],a5%s[0],a6%s[0]})", model_name.c_str(), ul, ul, ul, ul, ul, ul, ul));
                break;
            }
            case Model::Poly8:
            {
                w.factory(Form("Polynomial::%s(mass,{a0%s[0],a1%s[0],a2%s[0],a3%s[0],a4%s[0],a5%s[0],a6%s[0],a7%s[0],a8%s[0]})", model_name.c_str(), ul, ul, ul, ul, ul, ul, ul, ul, ul));
                break;
            }
            default:
                throw std::invalid_argument("[tnp::AddModelToWorkspace] Error: model not supported");
        }

        // done 
        return;
    }

    // get Model::value_type from a string
    // ----------------------------------------------------------------- //

    Model::value_type GetModelFromString(const std::string& model_name)
    {
        if(lt::string_lower(model_name) == "breitwignercb") {return Model::BreitWignerCB;}
        if(lt::string_lower(model_name) == "mctemplate"   ) {return Model::MCTemplate;   }
        if(lt::string_lower(model_name) == "exponential"  ) {return Model::Exponential;  }
        if(lt::string_lower(model_name) == "argus"        ) {return Model::Argus;        }
        if(lt::string_lower(model_name) == "erfexp"       ) {return Model::ErfExp;       }
        if(lt::string_lower(model_name) == "chebychev"    ) {return Model::Chebychev;    }
        if(lt::string_lower(model_name) == "chebyexp"     ) {return Model::ChebyExp;     }
        if(lt::string_lower(model_name) == "linear"       ) {return Model::Linear;       }
        if(lt::string_lower(model_name) == "poly2"        ) {return Model::Poly2;        }
        if(lt::string_lower(model_name) == "poly3"        ) {return Model::Poly3;        }
        if(lt::string_lower(model_name) == "poly6"        ) {return Model::Poly6;        }
        if(lt::string_lower(model_name) == "poly8"        ) {return Model::Poly8;        }

        // if here, didn't find a match
        throw std::invalid_argument(Form("[tnp::GetModelFromString] Error: model %s not found", model_name.c_str()));
    }

    // get Model::value_type from a string
    // ----------------------------------------------------------------- //

    std::string GetStringFromModel(const Model::value_type model)
    {
            switch(model)
            {
                case Model::BreitWignerCB : return "BreitWignerCB";
                case Model::MCTemplate    : return "MCTemplate";
                case Model::Exponential   : return "Exponential";
                case Model::Argus         : return "Argus";
                case Model::ErfExp        : return "ErfExp";
                case Model::Chebychev     : return "Chebychev";
                case Model::ChebyExp      : return "ChebyExp";
                case Model::Linear        : return "Linear";
                case Model::Poly2         : return "Poly2";
                case Model::Poly3         : return "Poly3";
                case Model::Poly6         : return "Poly6";
                case Model::Poly8         : return "Poly8";
                default:
                    throw std::invalid_argument("[tnp::GetStringFromModel] Error: model not found");
            }
    }

    // simple class to hold the results
    // ----------------------------------------------------------------- //

    // construct:
    Result::Result
    (
        const double efficiency, 
        const double efficiency_err, 
        const double numerator, 
        const double numerator_err, 
        const double denominator, 
        const double denominator_err, 
        const double background, 
        const double background_err
    )
        : cpass(new TCanvas("cpass", "cpass", 800, 600))
        , cfail(new TCanvas("cfail", "cfail", 800, 600))
    {
        eff.value = efficiency; 
        eff.error = efficiency_err; 
        num.value = numerator; 
        num.error = numerator_err; 
        den.value = denominator; 
        den.error = denominator_err; 
        bkg.value = background; 
        bkg.error = background_err;
    }

    Result::Result()
        : cpass(new TCanvas("cpass", "cpass", 800, 600))
        , cfail(new TCanvas("cfail", "cfail", 800, 600))
    {
        eff.value = -999999.0;
        eff.error = -999999.0;
        num.value = -999999.0;
        num.error = -999999.0;
        den.value = -999999.0;
        den.error = -999999.0;
        bkg.value = -999999.0;
        bkg.error = -999999.0;
    }

    Result::~Result()
    {
        // fixme: not working quite right
//         delete cpass;
//         delete cfail;
    }

    std::string Result::eff_str() const
    {
        return lt::pm(eff.value, eff.value, "1.3");
    }

    std::string Result::num_str() const
    {
        return lt::pm(num.value, num.value, "1.3");
    }

    std::string Result::den_str() const
    {
        return lt::pm(den.value, den.value, "1.3");
    }

    std::string Result::bkg_str() const
    {
        return lt::pm(bkg.value, bkg.value, "1.3");
    }


    // helper function to create a text box
    TPaveText* CreateTextBox(double x1, double y1, double x2, double y2, const std::string& text, const Color_t color = kBlack)
    {
        TPaveText *text_box = new TPaveText(x1, y1, x2, y2, "NDC");
        text_box->SetTextColor(color);
        text_box->SetFillStyle(0);
        text_box->SetBorderSize(0);
        text_box->SetTextAlign(12);
        text_box->AddText(text.c_str());
        return text_box;
    }

    // Perform the fit (with factor)
    Result PerformSimultaneousFit
    (
        const Model::value_type sig_pass_model, 
        const Model::value_type sig_fail_model, 
        const Model::value_type bkg_pass_model, 
        const Model::value_type bkg_fail_model, 
        const TH1* const h_pass, 
        const TH1* const h_fail,
        const float mass_low,
        const float mass_high,
        const float mass_bin_width,
        const std::string a_bin_label, 
        const std::string b_bin_label, 
        TH1* const h_pass_template,
        TH1* const h_fail_template
    )
    {
        // print the models used:
        cout << "[tnp::PerformSumultaneousFit] Fitting using the following PDF models:" << endl;
        cout << "sig pass model = " << tnp::GetStringFromModel(sig_pass_model) << endl; 
        cout << "sig fail model = " << tnp::GetStringFromModel(sig_fail_model) << endl; 
        cout << "bkg pass model = " << tnp::GetStringFromModel(bkg_pass_model) << endl; 
        cout << "bkg fail model = " << tnp::GetStringFromModel(bkg_fail_model) << endl; 

        // test template histogram's existence
        if (sig_pass_model == Model::MCTemplate && h_pass_template == NULL)
        {
            throw std::invalid_argument("[tnp::PerformSumultaneousFit] Error: pass template histogram is NULL!");
        }
        if (sig_fail_model == Model::MCTemplate && h_fail_template == NULL)
        {
            throw std::invalid_argument("[tnp::PerformSumultaneousFit] Error: fail template histogram is NULL!");
        }

        // build the RooWorkspace
        // ----------------------------------------- // 

        // workspace
        const std::string w_name = lt::string_replace_all(h_pass->GetName(), "h_", "w_name_");
        RooWorkspace w(w_name.c_str(), /*doCINTexport=*/false);

        // add independent variable (mass)
        RooRealVar mass("mass", "mass", mass_low, mass_high, "GeV");
        w.import(mass);

        // signal model (need ptr for polymorphism to work -- references for convenience)
//         const std::string sig_pass_model_name = lt::string_replace_all(h_pass->GetName(), "h_", "model_sig_"); 
//         const std::string sig_fail_model_name = lt::string_replace_all(h_fail->GetName(), "h_", "model_sig_"); 
//         const std::string bkg_pass_model_name = lt::string_replace_all(h_pass->GetName(), "h_", "model_bkg_"); 
//         const std::string bkg_fail_model_name = lt::string_replace_all(h_fail->GetName(), "h_", "model_bkg_"); 
        const std::string sig_pass_model_name = "model_sig_pass"; 
        const std::string sig_fail_model_name = "model_sig_fail"; 
        const std::string bkg_pass_model_name = "model_bkg_pass"; 
        const std::string bkg_fail_model_name = "model_bkg_fail"; 

        // add template hist
        if (sig_pass_model == Model::MCTemplate)
        {
            const std::string h_pass_newname = lt::string_replace_first(sig_pass_model_name, "model_", "h_template_"); 
            TH1* const h_pass_template_temp = dynamic_cast<TH1*>(h_pass_template->Clone(h_pass_newname.c_str()));
            w.import(*h_pass_template_temp);
            delete h_pass_template_temp;
        }
        if (sig_fail_model == Model::MCTemplate)
        {
            const std::string h_fail_newname = lt::string_replace_first(sig_fail_model_name, "model_", "h_template_"); 
            TH1* const h_fail_template_temp = dynamic_cast<TH1*>(h_fail_template->Clone(h_fail_newname.c_str()));
            w.import(*h_fail_template_temp);
            delete h_fail_template_temp;
        }

        AddModelToWorkspace(sig_pass_model, w, sig_pass_model_name); 
        AddModelToWorkspace(sig_fail_model, w, sig_fail_model_name); 
        AddModelToWorkspace(bkg_pass_model, w, bkg_pass_model_name); 
        AddModelToWorkspace(bkg_fail_model, w, bkg_fail_model_name); 
        w.Print();
        
        // extract the PDF's from the workspace 
        // ----------------------------------------- // 

        RooAbsPdf& spass_model = *w.pdf(sig_pass_model_name.c_str());
        RooAbsPdf& sfail_model = *w.pdf(sig_fail_model_name.c_str());
        RooAbsPdf& bpass_model = *w.pdf(bkg_pass_model_name.c_str());
        RooAbsPdf& bfail_model = *w.pdf(bkg_fail_model_name.c_str());

        // do the simultaneous fit 
        // ----------------------------------------- // 

        // Define categories
        RooCategory sample("sample","");
        sample.defineType("pass",1);
        sample.defineType("fail",2);

        // count maximums
        double nsig_max      = h_pass->Integral() + h_fail->Integral();
        double nbkg_fail_max = h_fail->Integral();
        double nbkg_pass_max = h_pass->Integral();

        RooRealVar eff("eff", "Efficiency", 0.8 ,0 ,1.0);
        RooRealVar nsig("nsig", "Signal Yield", 0.80*nsig_max, 0, nsig_max);
        RooRealVar nbkg_pass("nbkg_pass","Background count in PASS sample", 50, 0, nbkg_pass_max);
        RooRealVar nbkg_fail("nbkg_fail","Background count in FAIL sample", 0.1*nbkg_fail_max, 0.01, nbkg_fail_max);
        RooFormulaVar nsig_pass("nsig_pass" ,"eff*nsig"      , RooArgList(eff, nsig));
        RooFormulaVar nsig_fail("nsig_fail" ,"(1.0-eff)*nsig", RooArgList(eff, nsig));
        RooFormulaVar ntot_pass("ntot_pass" ,"nsig_pass+nbkg_pass", RooArgList(nsig_pass, nbkg_pass));
        RooFormulaVar ntot_fail("ntot_fail" ,"nsig_fail+nbkg_fail", RooArgList(nsig_fail, nbkg_fail));

        // add signal and background PDF for the total shape
        RooExtendPdf esig_pass("esig_pass","esig_pass", spass_model, nsig_pass);
        RooExtendPdf ebkg_pass("ebkg_pass","ebkg_pass", bpass_model, nbkg_pass);
        RooAddPdf model_pass("model_pass","Model for PASS sample", RooArgList(esig_pass, ebkg_pass));

        RooExtendPdf esig_fail("esig_fail","esig_fail", sfail_model, nsig_fail);
        RooExtendPdf ebkg_fail("ebkg_fail","ebkg_fail", bfail_model, nbkg_fail);
        RooAddPdf model_fail("model_fail","Model for Fail sample", RooArgList(esig_fail, ebkg_fail));

        // import data histograms
        RooDataHist data_pass("data_pass", "data_pass", RooArgSet(mass), h_pass);
        RooDataHist data_fail("data_fail", "data_fail", RooArgSet(mass), h_fail);
        RooDataHist data_comb("data_comb", "data_comb", RooArgList(mass), RooFit::Index(sample), RooFit::Import("pass", data_pass), RooFit::Import("fail", data_fail));  

        // for the total PDF
        RooSimultaneous total_pdf("total_pdf","total_pdf",sample);
        total_pdf.addPdf(model_pass, "pass");  
        total_pdf.addPdf(model_fail, "fail");

        // do the fit
        RooFitResult& roo_fit_result = *total_pdf.fitTo(data_comb, RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

        // integrate on a subrange
        mass.setRange("zwindow", mass_low, mass_high);
        RooAbsReal& int_nsig_pass = *esig_pass.createIntegral (mass, RooFit::NormSet(mass), RooFit::Range("zwindow"));
        RooAbsReal& int_nbkg_pass = *ebkg_pass.createIntegral (mass, RooFit::NormSet(mass), RooFit::Range("zwindow"));
        RooAbsReal& int_ntot_pass = *model_pass.createIntegral(mass, RooFit::NormSet(mass), RooFit::Range("zwindow"));
        RooAbsReal& int_nsig_fail = *esig_fail.createIntegral (mass, RooFit::NormSet(mass), RooFit::Range("zwindow"));
        RooAbsReal& int_nbkg_fail = *ebkg_fail.createIntegral (mass, RooFit::NormSet(mass), RooFit::Range("zwindow"));
        RooAbsReal& int_ntot_fail = *model_fail.createIntegral(mass, RooFit::NormSet(mass), RooFit::Range("zwindow"));

        // conventience varariables for the full fit counts
        Result::value_t npass_sig = {nsig_pass.getVal(), nsig_pass.getPropagatedError(roo_fit_result)};
        Result::value_t npass_bkg = {nbkg_pass.getVal(), nbkg_pass.getPropagatedError(roo_fit_result)};
        Result::value_t npass_tot = {ntot_pass.getVal(), ntot_pass.getPropagatedError(roo_fit_result)};
        Result::value_t nfail_sig = {nsig_fail.getVal(), nsig_fail.getPropagatedError(roo_fit_result)};
        Result::value_t nfail_bkg = {nbkg_fail.getVal(), nbkg_fail.getPropagatedError(roo_fit_result)};
        Result::value_t nfail_tot = {ntot_fail.getVal(), ntot_fail.getPropagatedError(roo_fit_result)};
        Result::value_t data_eff  = {eff.getVal()      , eff.getPropagatedError(roo_fit_result)      };

        // the value integrated on the subrange
        Result::value_t int_npass_sig = {int_nsig_pass.getVal() * npass_sig.value, int_nsig_pass.getVal() * npass_sig.error};
        Result::value_t int_npass_bkg = {int_nbkg_pass.getVal() * npass_bkg.value, int_nbkg_pass.getVal() * npass_bkg.error};
        Result::value_t int_npass_tot = {int_ntot_pass.getVal() * npass_tot.value, int_ntot_pass.getVal() * npass_tot.error};
        Result::value_t int_nfail_sig = {int_nsig_fail.getVal() * nfail_sig.value, int_nsig_fail.getVal() * nfail_sig.error};
        Result::value_t int_nfail_bkg = {int_nbkg_fail.getVal() * nfail_bkg.value, int_nbkg_fail.getVal() * nfail_bkg.error};
        Result::value_t int_nfail_tot = {int_ntot_fail.getVal() * nfail_tot.value, int_ntot_fail.getVal() * nfail_tot.error};

        // the new efficiency
        Result::value_t num = int_npass_sig;
        Result::value_t den = {int_npass_sig.value + int_nfail_sig.value, sqrt(pow(int_npass_sig.error,2) + pow(int_nfail_sig.error,2))};
        Result::value_t int_data_eff = {num.value/den.value, data_eff.error};
 
        // results of eff */
        Result simple_result
        (
             data_eff.value, 
             data_eff.error,
             num.value,
             num.error,
             den.value,
             den.error,
             int_npass_bkg.value,
             int_npass_bkg.error
        );
 
        // passing plot
        simple_result.cpass->cd(); 
        simple_result.cpass->SetName(lt::string_replace_all(h_pass->GetName(), "h_", "canvas_").c_str());
        simple_result.cpass->SetTitle(simple_result.cpass->GetName());
        RooPlot *mframe_pass = mass.frame(RooFit::Bins(static_cast<int>(fabs(mass_high-mass_low)/2.0)));
        data_pass.plotOn(mframe_pass, RooFit::MarkerStyle(kFullCircle), RooFit::MarkerSize(0.8), RooFit::DrawOption("ZP"));
        model_pass.plotOn(mframe_pass);
        model_pass.plotOn(mframe_pass,RooFit::Components("ebkg_pass"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
        mframe_pass->SetTitle("Passing Probes");
        mframe_pass->Draw();
        TPaveText *a_box         = CreateTextBox(0.15, 0.80, 0.41, 0.85, a_bin_label); a_box->Draw();
        TPaveText *b_box         = CreateTextBox(0.15, 0.75, 0.41, 0.80, b_bin_label); b_box->Draw();
        TPaveText *npass_box     = CreateTextBox(0.15, 0.70, 0.33, 0.75, Form("%1.0f Events", nbkg_pass_max)); npass_box->Draw();
        TPaveText *eff_box       = CreateTextBox(0.65, 0.80, 0.85, 0.84, Form("#varepsilon = %1.3f #pm %1.3f"   , int_data_eff.value , int_data_eff.error )); eff_box->Draw();
        TPaveText *nsig_pass_box = CreateTextBox(0.65, 0.75, 0.85, 0.79, Form("N^{pass}_{sig} = %1.0f #pm %1.0f", npass_sig.value    , npass_sig.error    )); nsig_pass_box->Draw();
        TPaveText *nbkg_pass_box = CreateTextBox(0.65, 0.70, 0.85, 0.74, Form("N^{pass}_{bkg} = %1.0f #pm %1.0f", npass_bkg.value    , npass_bkg.error    )); nbkg_pass_box->Draw();
        TPaveText *ntot_pass_box = CreateTextBox(0.65, 0.65, 0.85, 0.69, Form("N^{pass}_{tot} = %1.0f #pm %1.0f", npass_tot.value    , npass_tot.error    )); ntot_pass_box->Draw();
        TPaveText *nsig_fail_box = CreateTextBox(0.65, 0.60, 0.85, 0.64, Form("N^{fail}_{sig} = %1.0f #pm %1.0f", nfail_sig.value    , nfail_sig.error    )); nsig_fail_box->Draw();
        TPaveText *nbkg_fail_box = CreateTextBox(0.65, 0.55, 0.85, 0.59, Form("N^{fail}_{bkg} = %1.0f #pm %1.0f", nfail_bkg.value    , nfail_bkg.error    )); nbkg_fail_box->Draw();
        TPaveText *ntot_fail_box = CreateTextBox(0.65, 0.50, 0.85, 0.54, Form("N^{fail}_{tot} = %1.0f #pm %1.0f", nfail_tot.value    , nfail_tot.error    )); ntot_fail_box->Draw();

        // failing plot
        simple_result.cfail->cd(); 
        simple_result.cfail->SetName(lt::string_replace_all(h_fail->GetName(), "h_", "canvas_").c_str());
        simple_result.cfail->SetTitle(simple_result.cfail->GetName());
        RooPlot *mframe_fail = mass.frame(RooFit::Bins(static_cast<int>(fabs(mass_high-mass_low)/mass_bin_width)));
        data_fail.plotOn(mframe_fail, RooFit::MarkerStyle(kFullCircle), RooFit::MarkerSize(0.8), RooFit::DrawOption("ZP"));
        model_fail.plotOn(mframe_fail);
        model_fail.plotOn(mframe_fail,RooFit::Components("ebkg_fail"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
        mframe_fail->SetTitle("Failing Probes");
        mframe_fail->Draw();
        a_box->Draw();
        b_box->Draw();
        TPaveText *nfail_box = CreateTextBox(0.15, 0.70, 0.33, 0.75, Form("%1.0f Events", nbkg_fail_max)); nfail_box->Draw();
        eff_box->Draw();
        nsig_pass_box->Draw();
        nbkg_pass_box->Draw();
        ntot_pass_box->Draw();
        nsig_fail_box->Draw();
        nbkg_fail_box->Draw();
        ntot_fail_box->Draw();

        // print the results
        CTable t1;
        t1.useTitle();
        t1.setTitle("data fit results full range");
        t1.setTable() (                         "value",         "error")
                      ( "signal_pass" , npass_sig.value, npass_sig.error)
                      ( "signal_fail" , nfail_sig.value, nfail_sig.error)
                      ( "bkg_pass"    , npass_bkg.value, npass_bkg.error)
                      ( "bkg_fail"    , nfail_bkg.value, nfail_bkg.error)
                      ( "total_pass"  , npass_tot.value, npass_tot.error)
                      ( "total_fail"  , nfail_tot.value, nfail_tot.error)
                      ( "eff"         , data_eff.value , data_eff.error );
        t1.print();

        CTable t2;
        t2.useTitle();
        t2.setTitle("data fit results subrange");
        t2.setTable() (                                 "value",         "error")
                      ( "signal_pass" , int_npass_sig.value, int_npass_sig.error)
                      ( "signal_fail" , int_nfail_sig.value, int_nfail_sig.error)
                      ( "bkg_pass"    , int_npass_bkg.value, int_npass_bkg.error)
                      ( "bkg_fail"    , int_nfail_bkg.value, int_nfail_bkg.error)
                      ( "total_pass"  , int_npass_tot.value, int_npass_tot.error)
                      ( "total_fail"  , int_nfail_tot.value, int_nfail_tot.error)
                      ( "eff"         , int_data_eff.value , int_data_eff.error );
        t2.print();
 
        // done 
        return simple_result;
    }

    // Perform a simple count 
    Result PerformSimpleCount
    (
        TH1* const h_pass, 
        TH1* const h_fail,
        const float mass_low,
        const float mass_high,
        const float mass_bin_width,
        const std::string a_bin_label, 
        const std::string b_bin_label
    )
    {
        const float zwin_low  = mass_low;
        const float zwin_high = mass_high;
        const unsigned int num_mass_bins = static_cast<unsigned int>(fabs(mass_high-mass_low)/mass_bin_width);

        // counts
        const double num_pass      = rt::IntegralAndError(h_pass, zwin_low, zwin_high).first;
        const double num_pass_err  = rt::IntegralAndError(h_pass, zwin_low, zwin_high).second;
        const double num_fail      = rt::IntegralAndError(h_fail, zwin_low, zwin_high).first;
        const double num_fail_err  = rt::IntegralAndError(h_fail, zwin_low, zwin_high).second;
        const double num_total     = num_pass + num_fail;
        const double num_total_err = sqrt(num_pass_err*num_pass_err + num_fail_err*num_fail_err);
        const double eff           = num_pass/num_total;
        const double eff_err       = sqrt(eff * (1.0 - eff)/num_total); // binomial
        Result simple_result(eff, eff_err, num_pass, num_pass_err, num_total, num_total_err, 0, 0); // taking symmetric for now

        // independent mass var
        RooRealVar mass("mass", "mass", mass_low, mass_high, "GeV");
        RooDataHist hist_pass("data_pass", "data_pass", RooArgSet(mass), h_pass);
        RooDataHist hist_fail("data_fail", "data_fail", RooArgSet(mass), h_fail);

        // passing plot
        simple_result.cpass->cd(); 
        RooPlot *mframe_pass = mass.frame(RooFit::Bins(num_mass_bins));
        hist_pass.plotOn(mframe_pass);
        mframe_pass->SetTitle("MC Passing Probes");
        mframe_pass->Draw();
        TPaveText *a_box     = CreateTextBox(0.15, 0.80, 0.41, 0.85, a_bin_label); a_box->Draw();
        TPaveText *b_box     = CreateTextBox(0.15, 0.75, 0.41, 0.80, b_bin_label); b_box->Draw();
        TPaveText *npass_box = CreateTextBox(0.15, 0.70, 0.33, 0.75, Form("%1.0f Events", h_pass->Integral())); npass_box->Draw();
        TPaveText *eff_box   = CreateTextBox(0.65, 0.80, 0.85, 0.84, "#varepsilon = " + simple_result.eff_str()); eff_box->Draw();
        TPaveText *num_box   = CreateTextBox(0.65, 0.75, 0.85, 0.79, Form("N_{pass} = %1.0f" , simple_result.num.value)); num_box->Draw();
        TPaveText *fail_box  = CreateTextBox(0.65, 0.70, 0.85, 0.74, Form("N_{fail} = %1.0f" , num_fail)); num_box->Draw();
        TPaveText *den_box   = CreateTextBox(0.65, 0.65, 0.85, 0.69, Form("N_{total} = %1.0f", simple_result.den.value)); den_box->Draw();

        // failing plot
        simple_result.cfail->cd(); 
        RooPlot *mframe_fail = mass.frame(RooFit::Bins(num_mass_bins));
        hist_fail.plotOn(mframe_fail);
        mframe_fail->SetTitle("MC Failing Probes");
        mframe_fail->Draw();
        a_box->Draw();
        b_box->Draw();
        TPaveText *nfail_box = CreateTextBox(0.15, 0.70, 0.33, 0.75, Form("%1.0f Events", h_fail->Integral())); nfail_box->Draw();
        eff_box->Draw();
        num_box->Draw();
        fail_box->Draw();
        den_box->Draw();

        return simple_result;
    }

} // namespace tnp
