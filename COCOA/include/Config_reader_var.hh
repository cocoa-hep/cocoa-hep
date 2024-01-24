#ifndef __CONFIG_READER_VAR_H__
#define __CONFIG_READER_VAR_H__

#include <string>
#include "G4NistManager.hh"

struct Geometry_definition
{
    int kMaxPixel;
    int kNLayers;
    // int Low_kNLayers;
    std::vector<std::vector<long double>> layer_noise_ECAL;
    std::vector<std::vector<long double>> layer_noise_HCAL;
    std::vector<long double> layer_noise;

    std::vector<std::vector<long double>> layer_deta_ECAL;
    std::vector<std::vector<long double>> layer_deta_HCAL;
    std::vector<long double> layer_deta_flatten;

    std::vector<std::vector<long double>> layer_dphi_ECAL;
    std::vector<std::vector<long double>> layer_dphi_HCAL;
    std::vector<long double> layer_dphi_flatten;

    std::vector<std::vector<long double>> layer_inn_radius_ECAL;
    std::vector<std::vector<long double>> layer_inn_radius_HCAL;
    std::vector<long double> layer_inn_radius_flatten;

    std::vector<std::vector<long double>> layer_mid_radius_ECAL;
    std::vector<std::vector<long double>> layer_mid_radius_HCAL;
    std::vector<long double> layer_mid_radius_flatten;

    std::vector<std::vector<long double>> layer_out_radius_ECAL;
    std::vector<std::vector<long double>> layer_out_radius_HCAL;
    std::vector<long double> layer_out_radius_flatten;

    std::vector<std::vector<int>> number_of_pixels_ECAL;
    std::vector<std::vector<int>> number_of_pixels_HCAL;
    std::vector<int> number_of_pixels_flatten;

    std::vector<std::vector<long double>> resolution_width_of_ECAL_layers_in_X0;
    std::vector<std::vector<long double>> resolution_width_of_HCAL_layers_in_Lambda_int;
};
struct Graph_construction
{
    std::vector<int>   max_samelayer_edges;  // how many nearest neighbors to create for cells in each layer
    std::vector<float> max_dr;               // dr range wherein we consider neighbors (same layer)
    std::vector<int>   max_interlayer_edges; // how many edges a cell can have to cells in other layers
    int                max_layer_sep;        // max number of layers that interlayer edges can span
    std::vector<int>   max_celltrack_edges;  // how many nearest cells in each layer to connect to a track
    std::vector<float> max_celltrack_dr;     // dr range wherein we consider track-cell connections, per layer
};
struct Topological_clustering_var
{
    float sigma_threshold_for_seed_cells;
    float sigma_threshold_for_neighboring_cells;
    float sigma_threshold_for_last_cells;
    float cluster_negative_energy_cells;
    float local_max_seed_energy;
};
struct Particle_flow_alg_var
{
    float S_discriminant_threshold;
    float E_div_p_threshold;
    std::vector<float> delta_Rprime_threshold;
    std::vector<float> momentum_delta_Rprime_threshold;
    float factor_sigma_E_div_p_template;
    float Moliere_radius;
};
struct Jet_parameters
{
    std::string algorithm;
    std::string recombination_scheme;
    long double radius;
    long double ptmin;
};
struct Fiducial_cuts
{
    float min_pT;  // remove particles at truth level with pT small than min_pT. Unit is GeV.
    float dR_cut;  // remove particles at truth level separated farther than dR_cut from the primary particle (used to generate single-jet data)
};

class Config_reader_var
{
private:
    Config_reader_var();

public:
    static Config_reader_var &GetInstance()
    {
        static Config_reader_var config;
        return config;
    };
    std::string Output_file_path;
    std::string Macro_file_path;
    std::string Type_of_running;
    bool Save_truth_particle_graph;
    bool Use_high_granularity;
    bool Skip_unuseable_tracks;
    bool check_geometry_overlap;
    bool check_geometry_overlap_only;
    bool use_inner_detector;
    bool use_ID_support;

    long double r_inn_calo;
    long double Layer_gap;
    long double fieldValue;
    long double max_eta_barrel;
    long double max_eta_endcap;
    long double max_phi;
    float       samplingFraction_ECAL;
    float       samplingFraction_HCAL;
    G4Material *Material_ECAL;
    G4Material *Material_HCAL;
    Geometry_definition low_resolution;
    Geometry_definition high_resolution;
    Graph_construction graph_construction;
    Topological_clustering_var topological_clustering;
    bool doPFlow;
    Particle_flow_alg_var particle_flow;
    Jet_parameters jet_parameter;
    Fiducial_cuts fiducial_cuts;

    bool doSuperclustering;

    bool run_hadron_test;
    bool run_piZero_test;
    bool run_jets_test;
    
};

#endif // __CONFIG_READER_VAR_H__
