#ifndef TNP_PERFORM_FITS_H
#define TNP_PERFORM_FITS_H

#include "TCanvas.h"
#include "TH1.h"
#include <string>

namespace tnp
{
    // models available
    struct Model
    {
        enum value_type
        {
            BreitWignerCB, // Breit-Wigner convolved with Crystal Ball function
            MCTemplate,    // MC Template, need to provide histogram 
            Exponential,   // Exponential 
            Argus,         // Argus function (http://en.wikipedia.org/wiki/ARGUS_distribution)
            ErfExp,        // Error function convolved with Exponential
            Linear,        // Linear function
            Poly2,         // 2nd order polynomial
            Poly3,         // 3rd order polynomial
            Poly6,         // 6th order polynomial
            Poly8,         // 8th order polynomial
            LinearExp,     // Linear convolved with exponential
            Poly2Exp,      // 2nd order polynomial convolved with exponential
            Poly3Exp,      // 3rd order polynomial convolved with exponential
            Poly4Exp,      // 4th order polynomial convolved with exponential
            Poly8Exp,      // 8th order polynomial convolved with exponential
            Chebychev,     // Chebychev polynomial 
            ChebyExp,      // Chebychev polynomial convolved with exponential
            static_size
        };
    };

    // get Model::value_type from a string
    Model::value_type GetModelFromString(const std::string& model_name);

    // get string from Model::value_type
    std::string GetStringFromModel(const Model::value_type model);

    // a simple struct to hold the results
    struct Result
    {
        // simple value/error pair
        struct value_t
        {
            double value;
            double error;
        };

        // construct:
        Result();
        Result
        (
            const double efficiency, 
            const double efficiency_err, 
            const double numerator, 
            const double numerator_err, 
            const double denominator, 
            const double denominator_err, 
            const double background, 
            const double background_err
        );

        // destroy:
        ~Result();

        // members:
        value_t eff;
        value_t num;
        value_t den;
        value_t bkg;
        TCanvas* cpass;
        TCanvas* cfail;

        // methods:
        std::string eff_str() const;
        std::string num_str() const;
        std::string den_str() const;
        std::string bkg_str() const;
    };

    // Get the model
    
    // Peform simultaneous fit
    Result PerformSimultaneousFit
    (
        const Model::value_type sig_pass_model, 
        const Model::value_type sig_fail_model, 
        const Model::value_type bkg_pass_model, 
        const Model::value_type bkg_fail_model, 
        const TH1* const h_pass, 
        const TH1* const h_fail,
        const float mass_low = 60.0,
        const float mass_high = 120.0,
        const float mass_bin_width = 2.0,
        const std::string a_bin_label = "", 
        const std::string b_bin_label = "", 
        TH1* const h_pass_template = NULL,
        TH1* const h_fail_template = NULL
    );

    // Peform smple count (no background fitting/subtraction) 
    Result PerformSimpleCount
    (
        TH1* const h_pass, 
        TH1* const h_fail,
        const float mlow = 60.0,
        const float mhigh = 120.0,
        const float mass_bin_width = 2.0,
        const std::string a_bin_label = "", 
        const std::string b_bin_label = "" 
    );

} // namespace tnp

#endif // TNP_PERFORM_FITS_H
