#ifndef TNP_DATASET_H
#define TNP_DATASET_H

#include <vector>
#include <string>

// foward declare
namespace edm {class ParameterSet;}

namespace tnp
{
    // simple class to hold the information relevant to the a dataset
    struct Dataset
    {
        // construct:
        Dataset(const edm::ParameterSet& pset);

        Dataset
        (
             const std::string& name,
             const std::string& title,
             const std::vector<std::string>& input_file_names,
             const std::string& run_list = "",
             const bool is_data = false
        );

        // members:
        std::string m_name;
        std::string m_title;
        std::vector<std::string> m_input_file_names;
        std::string m_run_list;
        bool m_is_data;
    };

    // get the datasets from the VPset
    std::vector<Dataset> GetDatasetsFromVPSet(const std::vector<edm::ParameterSet>& psets);

} // namespace tnp

#endif // TNP_DATASET_H
