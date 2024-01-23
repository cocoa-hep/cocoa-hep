#include "Config_reader_func.hh"

#include <iostream>
#include <fstream>

Config_reader_func::Config_reader_func(std::string path, Config_reader_var &config_var)
{
    
    std::ifstream config_doc(path, std::ifstream::binary);

    config_doc >> configs;
    
    config_var.Output_file_path = configs["Output_file_path"].asString();
    config_var.Type_of_running = configs["Type_of_running"].asString();
    config_var.Macro_file_path = configs["Macro_file_path"].asString();
    config_var.Save_truth_particle_graph = configs["Save_truth_particle_graph"].asBool();
    config_var.Use_high_granularity = configs["Use_high_granularity"].asBool();
    config_var.Skip_unuseable_tracks = configs["Skip_unuseable_tracks"].asBool();
    config_var.doSuperclustering = configs.get( "Do_superclustering", false ).asBool();

    config_var.r_inn_calo = configs["Geometry_definition"]["Inner_calorimeter_layer"].asDouble();
    config_var.Layer_gap = configs["Geometry_definition"]["Layer_gap"].asDouble();
    config_var.fieldValue = configs["Geometry_definition"]["Magnetic_field"].asDouble() * tesla;
    config_var.max_eta_barrel = configs["Geometry_definition"]["Max_eta_of_barrel_region"].asDouble();
    config_var.max_eta_endcap = configs["Geometry_definition"]["Max_eta_of_endcap_region"].asDouble();
    config_var.max_phi = 2. * M_PI * rad;
    config_var.check_geometry_overlap = configs["Geometry_definition"]["Check_Geometry_overlap"].asBool();
    config_var.check_geometry_overlap_only = configs["Geometry_definition"].get( "Check_Geometry_overlap_only", false ).asBool();
    if ( config_var.check_geometry_overlap_only )
	config_var.check_geometry_overlap = true;
    config_var.use_inner_detector = configs["Geometry_definition"].get( "Use_inner_detector", true ).asBool();
    config_var.use_ID_support = configs["Geometry_definition"].get( "Use_ID_support", true ).asBool();

    Json::Value &layervals = configs["Graph_construction"]["max_samelayer_edges"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_samelayer_edges);
                 layervals = configs["Graph_construction"]["max_dr"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_dr);
                 layervals = configs["Graph_construction"]["max_interlayer_edges"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_interlayer_edges);
    config_var.graph_construction.max_layer_sep = configs["Graph_construction"]["max_layer_sep"].asInt();
                 layervals = configs["Graph_construction"]["max_celltrack_edges"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_celltrack_edges);
                 layervals = configs["Graph_construction"]["max_celltrack_dr"];
    Fill_1D_vector(layervals, config_var.graph_construction.max_celltrack_dr);

    config_var.topological_clustering.sigma_threshold_for_seed_cells = configs["TopoClustering"]["sigma_threshold_for_seed_cells"].asFloat();
    config_var.topological_clustering.sigma_threshold_for_neighboring_cells = configs["TopoClustering"]["sigma_threshold_for_neighboring_cells"].asFloat();
    config_var.topological_clustering.sigma_threshold_for_last_cells = configs["TopoClustering"]["sigma_threshold_for_last_cells"].asFloat();
    config_var.topological_clustering.cluster_negative_energy_cells = configs["TopoClustering"]["cluster_negative_energy_cells"].asBool();
    config_var.topological_clustering.local_max_seed_energy = configs["TopoClustering"]["local_max_seed_energy"].asFloat();
    
    config_var.particle_flow.S_discriminant_threshold = configs["Particle_flow"].get( "S_discriminant_threshold", -1.0 ).asFloat();
    config_var.particle_flow.E_div_p_threshold = configs["Particle_flow"].get( "E_div_p_threshold", 0.1 ).asFloat();
    config_var.particle_flow.factor_sigma_E_div_p_template = configs["Particle_flow"].get( "factor_sigma_E_div_p_template", 1.25 ).asFloat();
    config_var.particle_flow.Moliere_radius = configs["Particle_flow"].get( "Moliere_radius", 0.035 ).asFloat();
    
    config_var.jet_parameter.algorithm = configs["Jet_parameters"]["algorithm"].asString();
    config_var.jet_parameter.recombination_scheme = configs["Jet_parameters"]["recombination_scheme"].asString();
    config_var.jet_parameter.ptmin = configs["Jet_parameters"]["ptmin"].asDouble();
    config_var.jet_parameter.radius = configs["Jet_parameters"]["radius"].asDouble();

    config_var.fiducial_cuts.min_pT = configs["Fiducial_cuts"].get( "min_pT", -1.0 ).asFloat();
    config_var.fiducial_cuts.dR_cut = configs["Fiducial_cuts"].get( "dR_cut", -1.0 ).asFloat();

    if (config_var.Type_of_running != "Standard")
    {
        std::cout<<config_var.Type_of_running<<std::endl;
        config_var.Use_high_granularity = false;
    }

    // ==============================================================
    // *Materials
    // ==============================================================
    config_var.Material_ECAL = Material_build("ECAL");
    config_var.Material_HCAL = Material_build("HCAL");

    //* Fill low resolution
    Json::Value &characters = configs["Geometry_definition"]["Detector_granularity"]["Number_of_pixels_ECAL"];
    Fill_1D_vector(characters, config_var.low_resolution.number_of_pixels_ECAL); // Low_number_of_pixels_ECAL
    characters = configs["Geometry_definition"]["Detector_granularity"]["Number_of_pixels_HCAL"];
    Fill_1D_vector(characters, config_var.low_resolution.number_of_pixels_HCAL); //Low_resolution_number_of_pixels_HCAL
    characters = configs["Geometry_definition"]["Detector_granularity"]["Width_of_ECAL_layers_in_X0"];
    Fill_1D_vector(characters, config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0); //Low_resolution_width_of_ECAL_layers_in_X0
    characters = configs["Geometry_definition"]["Detector_granularity"]["Width_of_HCAL_layers_in_Lambda_int"];
    config_var.samplingFraction_ECAL = configs["Geometry_definition"].get("SamplingFraction_ECAL", 0.02 ).asFloat();
    config_var.samplingFraction_HCAL = configs["Geometry_definition"].get("SamplingFraction_HCAL", 0.02 ).asFloat();
    Fill_1D_vector(characters, config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int); //Low_resolution_width_of_HCAL_layers_in_Lambda_int
    characters = configs["Geometry_definition"]["Noise_in_ECAL"];
    Fill_1D_vector(characters, config_var.low_resolution.layer_noise_ECAL); //Low_layer_noise_ECAL
    characters = configs["Geometry_definition"]["Noise_in_HCAL"];
    Fill_1D_vector(characters, config_var.low_resolution.layer_noise_HCAL); //Low_layer_noise_HCAL
    
    if ( configs["Particle_flow"].isObject() ) {
	config_var.doPFlow = true;
	characters = configs["Particle_flow"]["delta_Rprime_threshold"];
	Fill_1D_vector(characters, config_var.particle_flow.delta_Rprime_threshold);
	characters = configs["Particle_flow"]["momentum_delta_Rprime_threshold"];
	Fill_1D_vector(characters, config_var.particle_flow.momentum_delta_Rprime_threshold);
    } else {
	config_var.doPFlow = false;
	config_var.particle_flow.delta_Rprime_threshold = { 2.4, 1.25, 0.8 };
	config_var.particle_flow.momentum_delta_Rprime_threshold = { 2000, 5000 };
    }
    
    config_var.low_resolution.layer_inn_radius_ECAL = config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0;
    config_var.low_resolution.layer_mid_radius_ECAL = config_var.low_resolution.layer_inn_radius_ECAL;
    config_var.low_resolution.layer_out_radius_ECAL = config_var.low_resolution.layer_mid_radius_ECAL;
    config_var.low_resolution.layer_deta_ECAL = config_var.low_resolution.layer_out_radius_ECAL;
    config_var.low_resolution.layer_dphi_ECAL = config_var.low_resolution.layer_deta_ECAL;
    config_var.low_resolution.layer_inn_radius_HCAL = config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int;
    config_var.low_resolution.layer_mid_radius_HCAL = config_var.low_resolution.layer_inn_radius_HCAL;
    config_var.low_resolution.layer_out_radius_HCAL = config_var.low_resolution.layer_mid_radius_HCAL;
    config_var.low_resolution.layer_deta_HCAL = config_var.low_resolution.layer_out_radius_HCAL;
    config_var.low_resolution.layer_dphi_HCAL = config_var.low_resolution.layer_deta_HCAL;
    
    int kNLayers = 0;
    int nLow_Layers = config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0.size();
    kNLayers += nLow_Layers;
    long double r_inn = config_var.r_inn_calo;
    for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
    {
        config_var.low_resolution.layer_noise.push_back(config_var.low_resolution.layer_noise_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_ECAL.at(ilow_layer).at(0) = 2 * config_var.max_eta_endcap / config_var.low_resolution.number_of_pixels_ECAL.at(ilow_layer).at(0); // Low_number_of_pixels_ECAL.at(ilow_layer);
        config_var.low_resolution.layer_dphi_ECAL.at(ilow_layer).at(0) = config_var.max_phi / config_var.low_resolution.number_of_pixels_ECAL.at(ilow_layer).at(0);
        config_var.low_resolution.number_of_pixels_flatten.push_back(config_var.low_resolution.number_of_pixels_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_flatten.push_back(config_var.low_resolution.layer_deta_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_dphi_flatten.push_back(config_var.low_resolution.layer_dphi_ECAL.at(ilow_layer).at(0));
        long double r_out = r_inn + config_var.low_resolution.resolution_width_of_ECAL_layers_in_X0.at(ilow_layer).at(0) * config_var.Material_ECAL->GetRadlen();
        config_var.low_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(0) = r_inn;
        config_var.low_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(0) = r_out - 0.5 * (r_out - r_inn);
        config_var.low_resolution.layer_out_radius_ECAL.at(ilow_layer).at(0) = r_out;
        config_var.low_resolution.layer_inn_radius_flatten.push_back(config_var.low_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_mid_radius_flatten.push_back(config_var.low_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_out_radius_flatten.push_back(config_var.low_resolution.layer_out_radius_ECAL.at(ilow_layer).at(0));
        r_inn = r_out;
    }
    r_inn = r_inn + config_var.Layer_gap;
    nLow_Layers = config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int.size();
    kNLayers += nLow_Layers;
    for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
    {
        config_var.low_resolution.layer_noise.push_back(config_var.low_resolution.layer_noise_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_HCAL.at(ilow_layer).at(0) = 2 * config_var.max_eta_endcap / config_var.low_resolution.number_of_pixels_HCAL.at(ilow_layer).at(0);
        config_var.low_resolution.layer_dphi_HCAL.at(ilow_layer).at(0) = config_var.max_phi / config_var.low_resolution.number_of_pixels_HCAL.at(ilow_layer).at(0);
        config_var.low_resolution.number_of_pixels_flatten.push_back(config_var.low_resolution.number_of_pixels_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_dphi_flatten.push_back(config_var.low_resolution.layer_dphi_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_deta_flatten.push_back(config_var.low_resolution.layer_deta_HCAL.at(ilow_layer).at(0));
        long double r_out = r_inn + config_var.low_resolution.resolution_width_of_HCAL_layers_in_Lambda_int.at(ilow_layer).at(0) * config_var.Material_HCAL->GetNuclearInterLength();
        config_var.low_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(0) = r_inn;
        config_var.low_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(0) = r_out - 0.5 * (r_out - r_inn);
        config_var.low_resolution.layer_out_radius_HCAL.at(ilow_layer).at(0) = r_out;
        config_var.low_resolution.layer_inn_radius_flatten.push_back(config_var.low_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_mid_radius_flatten.push_back(config_var.low_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(0));
        config_var.low_resolution.layer_out_radius_flatten.push_back(config_var.low_resolution.layer_out_radius_HCAL.at(ilow_layer).at(0));
        r_inn = r_out;
    }
    //* Fill low resolution end
    config_var.low_resolution.kNLayers = kNLayers;
    if (config_var.Use_high_granularity)
    {
        kNLayers = 0;
        //* Fill high resolution
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Number_of_pixels_ECAL"];
        Fill_2D_vector(characters, config_var.high_resolution.number_of_pixels_ECAL);
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Number_of_pixels_HCAL"];
        Fill_2D_vector(characters, config_var.high_resolution.number_of_pixels_HCAL);
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Width_of_ECAL_layers_in_X0"];
        Fill_2D_vector(characters, config_var.high_resolution.resolution_width_of_ECAL_layers_in_X0);
        characters = configs["Geometry_definition"]["High_Granularity_detector"]["Width_of_HCAL_layers_in_Lambda_int"];
        Fill_2D_vector(characters, config_var.high_resolution.resolution_width_of_HCAL_layers_in_Lambda_int);
        //* Fill high resolution end
        nLow_Layers = config_var.high_resolution.number_of_pixels_ECAL.size();
        config_var.high_resolution.layer_inn_radius_ECAL = config_var.high_resolution.resolution_width_of_ECAL_layers_in_X0;
        config_var.high_resolution.layer_mid_radius_ECAL = config_var.high_resolution.layer_inn_radius_ECAL;
        config_var.high_resolution.layer_out_radius_ECAL = config_var.high_resolution.layer_mid_radius_ECAL;
        config_var.high_resolution.layer_deta_ECAL = config_var.high_resolution.layer_out_radius_ECAL;
        config_var.high_resolution.layer_dphi_ECAL = config_var.high_resolution.layer_deta_ECAL;
        config_var.high_resolution.layer_inn_radius_HCAL = config_var.high_resolution.resolution_width_of_HCAL_layers_in_Lambda_int;
        config_var.high_resolution.layer_mid_radius_HCAL = config_var.high_resolution.layer_inn_radius_HCAL;
        config_var.high_resolution.layer_out_radius_HCAL = config_var.high_resolution.layer_mid_radius_HCAL;
        config_var.high_resolution.layer_deta_HCAL = config_var.high_resolution.layer_out_radius_HCAL;
        config_var.high_resolution.layer_dphi_HCAL = config_var.high_resolution.layer_deta_HCAL;
        r_inn = config_var.r_inn_calo;
        for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
        {

            int nHigh_Layers = config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).size();
            kNLayers += nHigh_Layers;
            for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
            {
                config_var.high_resolution.layer_deta_ECAL.at(ilow_layer).at(ihigh_layer) = 2 * config_var.max_eta_endcap / config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).at(ihigh_layer);
                config_var.high_resolution.layer_dphi_ECAL.at(ilow_layer).at(ihigh_layer) = config_var.max_phi / config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).at(ihigh_layer);
                config_var.high_resolution.number_of_pixels_flatten.push_back(config_var.high_resolution.number_of_pixels_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_dphi_flatten.push_back(config_var.high_resolution.layer_dphi_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_deta_flatten.push_back(config_var.high_resolution.layer_deta_ECAL.at(ilow_layer).at(ihigh_layer));
                long double r_out = r_inn + config_var.high_resolution.resolution_width_of_ECAL_layers_in_X0.at(ilow_layer).at(ihigh_layer) * config_var.Material_ECAL->GetRadlen();
                config_var.high_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(ihigh_layer) = r_inn;
                config_var.high_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(ihigh_layer) = r_out - 0.5 * (r_out - r_inn);
                config_var.high_resolution.layer_out_radius_ECAL.at(ilow_layer).at(ihigh_layer) = r_out;
                config_var.high_resolution.layer_inn_radius_flatten.push_back(config_var.high_resolution.layer_inn_radius_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_mid_radius_flatten.push_back(config_var.high_resolution.layer_mid_radius_ECAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_out_radius_flatten.push_back(config_var.high_resolution.layer_out_radius_ECAL.at(ilow_layer).at(ihigh_layer));
                r_inn = r_out;
            }
        }
        r_inn = r_inn + config_var.Layer_gap;
        nLow_Layers = config_var.high_resolution.number_of_pixels_HCAL.size();
        for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
        {
            int nHigh_Layers = config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).size();
            kNLayers += nHigh_Layers;
            for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
            {
                config_var.high_resolution.layer_deta_HCAL.at(ilow_layer).at(ihigh_layer) = 2 * config_var.max_eta_endcap / config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).at(ihigh_layer);
                config_var.high_resolution.layer_dphi_HCAL.at(ilow_layer).at(ihigh_layer) = config_var.max_phi / config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).at(ihigh_layer);
                long double r_out = r_inn + config_var.high_resolution.resolution_width_of_HCAL_layers_in_Lambda_int.at(ilow_layer).at(ihigh_layer) * config_var.Material_HCAL->GetNuclearInterLength(); // Material_H->GetRadlen();
                config_var.high_resolution.number_of_pixels_flatten.push_back(config_var.high_resolution.number_of_pixels_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_dphi_flatten.push_back(config_var.high_resolution.layer_dphi_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_deta_flatten.push_back(config_var.high_resolution.layer_deta_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(ihigh_layer) = r_inn;
                config_var.high_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(ihigh_layer) = r_out - 0.5 * (r_out - r_inn);
                config_var.high_resolution.layer_out_radius_HCAL.at(ilow_layer).at(ihigh_layer) = r_out;
                config_var.high_resolution.layer_inn_radius_flatten.push_back(config_var.high_resolution.layer_inn_radius_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_mid_radius_flatten.push_back(config_var.high_resolution.layer_mid_radius_HCAL.at(ilow_layer).at(ihigh_layer));
                config_var.high_resolution.layer_out_radius_flatten.push_back(config_var.high_resolution.layer_out_radius_HCAL.at(ilow_layer).at(ihigh_layer));
                r_inn = r_out;
            }
        }
    }
    
    config_var.high_resolution.kNLayers = kNLayers;

    config_var.run_hadron_test = configs.get( "run_hadron_test", false ).asBool();
    config_var.run_piZero_test = configs.get( "run_piZero_test", false ).asBool();
    config_var.run_jets_test   = configs.get( "run_jets_test", false ).asBool();
    size_t n_tests = 0;
    for ( bool test : { config_var.run_hadron_test,
		        config_var.run_piZero_test,
		config_var.run_jets_test } ) {
	if ( test)
	    ++n_tests;
    }
    if ( n_tests > 1 )
	G4cerr << "Warning: chose more than one test to perform, but tests are meant to be mutually exclusive." << G4endl;
    
}

G4Material *Config_reader_func::Material_build(std::string name)
{
    G4NistManager *nistManager = G4NistManager::Instance();
    long double a, z, density_liq, fracMass;
    std::vector<std::string> sMaterial_vector;
    std::vector<G4Material *> Material_vector;
    G4Material *Material;

    Json::Value &characters = configs["Geometry_definition"]["Material_for_" + name];
    Fill_1D_vector(characters, sMaterial_vector);
    int icustom_element = 0;
    for (unsigned long imat_ecal = 0; imat_ecal < sMaterial_vector.size(); imat_ecal++)
    {
        G4Material *element = nistManager->FindOrBuildMaterial(sMaterial_vector[imat_ecal]);
        if (!element)
        {
            element = new G4Material(sMaterial_vector[imat_ecal], 
                                     z = configs["Geometry_definition"]["Characteristic_of_custom_material_for_" + name][icustom_element][0].asDouble(),
                                     a = configs["Geometry_definition"]["Characteristic_of_custom_material_for_" + name][icustom_element][1].asDouble() * g / mole,
                                     density_liq = configs["Geometry_definition"]["Characteristic_of_custom_material_for_" + name][icustom_element][2].asDouble() * g / cm3);
        icustom_element++;
        }
        Material_vector.push_back(element);
    }
    long double numerator = 0;
    long double denominator = 0;
    for (unsigned long irho_ecal = 0; irho_ecal < configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"].size(); irho_ecal++)
    {
        numerator += Material_vector[irho_ecal]->GetDensity() * configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"][(int)irho_ecal].asDouble();
        denominator += configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"][(int)irho_ecal].asDouble();
    }
    Material = new G4Material("Material_" + name, numerator / denominator, sMaterial_vector.size());
    for (unsigned long imat_ecal = 0; imat_ecal < Material_vector.size(); imat_ecal++)
    {
        Material->AddMaterial(Material_vector[imat_ecal],
                              fracMass = (Material_vector[imat_ecal]->GetDensity() * configs["Geometry_definition"][name + "_material_mixing_in_volume_proportion"][(int)imat_ecal].asDouble()) / (numerator)*100 * perCent);
    }
    return Material;
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<std::vector<long double>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<long double> vec_elem;
        vec_elem.push_back(Json_list[i].asDouble());
        vec.push_back(vec_elem);
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<float> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        vec.push_back(Json_list[i].asDouble());
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<int> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        vec.push_back(Json_list[i].asInt());
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<std::vector<int>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<int> vec_elem;
        vec_elem.push_back(Json_list[i].asUInt());
        vec.push_back(vec_elem);
    }
}
void Config_reader_func::Fill_1D_vector(const Json::Value &Json_list, std::vector<std::string> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        vec.push_back(Json_list[i].asString());
    }
}
void Config_reader_func::Fill_2D_vector(const Json::Value &Json_list, std::vector<std::vector<int>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<int> vec_elem;
        vec_elem.clear();
        for (unsigned int j = 0; j < Json_list[i].size(); j++)
        {
            vec_elem.push_back(Json_list[i][j].asUInt());
        }
        vec.push_back(vec_elem);
    }
}
void Config_reader_func::Fill_2D_vector(const Json::Value &Json_list, std::vector<std::vector<long double>> &vec)
{
    vec.clear();
    for (unsigned int i = 0; i < Json_list.size(); i++)
    {
        std::vector<long double> vec_elem;
        vec_elem.clear();
        for (unsigned int j = 0; j < Json_list[i].size(); j++)
        {
            vec_elem.push_back(Json_list[i][j].asDouble());
        }
        vec.push_back(vec_elem);
    }
}
