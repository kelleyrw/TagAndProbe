#ifndef TNP_MEASUREMENT_H
#define TNP_MEASUREMENT_H

#include <string>

namespace tnp
{
    // simple lepton type
    struct Lepton
    {
        enum value_type 
        {
            Muon,
            Electron,
            Both,
            static_size
        };
    };

    // simple selection type
    struct Selection
    {
        enum value_type 
        {
            // EGamma
            EGammaDenID,
            EGammaDenIso,
            EGammaDenBoth,
            EGammaNum,

            // same sign electrons
            SameSignElDenID,
            SameSignElDenIso,
            SameSignElDenBoth,
            SameSignElNum,

            // same sign muons
            SameSignMuDenID,
            SameSignMuDenIso,
            SameSignMuDenBoth,
            SameSignMuNum,

            static_size
        };
    };

    // EGamma cut based
    struct EG
    {
        enum value_type
        {
            DETAIN        = (1<<0),
            DPHIIN        = (1<<1),
            SIGMAIETAIETA = (1<<2),
            HOE           = (1<<3),
            OOEMOOP       = (1<<4),
            D0VTX         = (1<<5),
            DZVTX         = (1<<6),
            ISO           = (1<<7),
            VTXFIT        = (1<<8),
            MHITS         = (1<<9),
            static_size
        };
    };

    // all possible cuts pass
    static const unsigned int PassAll = 
    (
        EG::DETAIN        | 
        EG::DPHIIN        | 
        EG::SIGMAIETAIETA | 
        EG::HOE           |
        EG::OOEMOOP       | 
        EG::D0VTX         | 
        EG::DZVTX         | 
        EG::ISO           | 
        EG::VTXFIT        | 
        EG::MHITS
    );

    // all possible cuts pass without isolation
    static const unsigned int PassNoIso = 
    (
        EG::DETAIN |
        EG::DPHIIN |
        EG::SIGMAIETAIETA |
        EG::HOE |
        EG::OOEMOOP |
        EG::D0VTX |
        EG::DZVTX |
        EG::VTXFIT |
        EG::MHITS
    );

    // get the lepton lepton from a string
    Lepton::value_type GetLeptonFromString(const std::string& lepton_name);

    // get the string from the lepton type 
    std::string GetStringFromLepton(const Lepton::value_type lepton_type);

    // get the selection from a string
    Selection::value_type GetSelectionFromString(const std::string& sel_name);

    // passes selection based on above enum
    bool PassesSelection(const Lepton::value_type lepton, const Selection::value_type selection, const bool is_data);

} // namespace tnp

#endif //TNP_MEASUREMENT_H
