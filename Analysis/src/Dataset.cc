#include "TagAndProbe/Analysis/interface/Dataset.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace tnp
{
    Dataset::Dataset
    (
        const std::string& name,
        const std::string& title,
        const std::vector<std::string>& input_file_names,
        const std::string& run_list,
        const bool is_data
    )
        : m_name(name)
        , m_title(title)
        , m_input_file_names(input_file_names)
        , m_run_list(run_list)
        , m_is_data(is_data)
    {
    }

    Dataset::Dataset(const edm::ParameterSet& pset)
        : m_name(pset.getParameter<std::string>("name"))
        , m_title(pset.getParameter<std::string>("title"))
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

} // namespace tnp
