#include "Detector_analysis_var.hh"

Detector_analysis_var::Detector_analysis_var()
{
    Config_reader_var config_obj = Config_reader_var::GetInstance();
    Geometry = config_obj.low_resolution;
    rad_length.resize(Geometry.kNLayers + 2);
    int_length.resize(Geometry.kNLayers + 2);
}

void Detector_analysis_var::set_eta(float Eta)
{
    eta = Eta;
}

float Detector_analysis_var::get_eta()
{
    return eta;
}

void Detector_analysis_var::add_lengths(float x0, float lambda, int layer)
{
    rad_length.at(layer) += x0;
    int_length.at(layer) += lambda;
}

void Detector_analysis_var::set_tree_branches(TTree *outTree, std::string type_of_running)
{
    if (type_of_running == "Detector_Geantino")
    {
        outTree->Branch("geantino_init_eta", &eta);
        outTree->Branch("geantino_radiation_length_layers", "vector<float>", &rad_length);
        outTree->Branch("geantino_interaction_length_layers", "vector<float>", &int_length);
    }
    if (type_of_running == "Detector_Response")
    {
        outTree->Branch("response_charge_transverse_energy", &charge_transverse_energy);
        outTree->Branch("response_neutral_transverse_energy", &neutral_transverse_energy);
        outTree->Branch("response_charge_energy", &charge_energy);
        outTree->Branch("response_neutral_energy", &neutral_energy);
    }
}

void Detector_analysis_var::transverse_energy_calculation(std::vector<std::vector<std::vector<Cell>>> fcell_array) 
{
    for (int ilay = 0; ilay < (int) fcell_array.size(); ilay++)
    {
        for (int ieta = 0; ieta < (int)fcell_array.at(ilay).size(); ieta++)
        {
            for (int iphi = 0; iphi < (int)fcell_array.at(ilay).at(ieta).size(); iphi++)
            {
                Cell &cell = fcell_array.at(ilay).at(ieta).at(iphi);
                charge_transverse_energy+=cell.get_charge_energy()/cosh(cell.get_eta_pos());
                neutral_transverse_energy+=cell.get_neutral_energy()/cosh(cell.get_eta_pos());
                charge_energy +=cell.get_charge_energy();
                neutral_energy +=cell.get_neutral_energy();
            }
        }
    }
}

void Detector_analysis_var::clear()
{
    rad_length.clear();
    int_length.clear();
    rad_length.resize(Geometry.kNLayers + 2);
    int_length.resize(Geometry.kNLayers + 2);
    eta += 0.001;
    charge_transverse_energy = 0;
    neutral_transverse_energy = 0;
    charge_energy  = 0;
    neutral_energy = 0;
}