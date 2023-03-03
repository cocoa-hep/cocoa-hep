#include "Cell_var.hh"
#include "DetectorConstruction.hh"

#include <math.h>

Cell::Cell()
{
    Reset();
}
Cell::Cell(const Cell &orig) = default;

void Cell::Reset()
{

    total_energy = 0;
    charge_energy = 0;
    neutral_energy = 0;
    noise = 0;
    sigma = 0;
    layer_idx = 0;
    eta_idx = 0;
    phi_idx = 0;
    label = 0;
    is_shared = false;
    first_local_max_weight = 0;
    second_local_max_weight = 0;
    first_local_max_label = 0;
    second_local_max_label = 0;
    size_of_seeds_list = 0;
    cell_link_superres.clear();
    _1st_local_max_indexes.clear();
    _2nd_local_max_indexes.clear();
    list_of_particles_dep_energy.clear();
    list_of_conv_electrons_dep_energy.clear();
    eta_pos = 0;
    theta_pos = 0;
    phi_pos = 0;
    x_pos = 0;
    y_pos = 0;
    z_pos = 0;
    abs_eta_min = 0;
    abs_eta_max = 0;
    abs_phi_min = 0;
    abs_phi_max = 0;
    pflow_remnant = 0;
    distance_to_track = 0;
    energy_density = 0;
    ring_number = 0;
    is_1st_local_max_assoisited_with_track = false;
    is_2nd_local_max_assoisited_with_track = false;
    truth_label = -1;
}

void Cell::subtract_energys_fraction_from_pflow_remnant(float fraction)
{
    if (is_shared)
    {
        if (is_1st_local_max_assoisited_with_track && is_2nd_local_max_assoisited_with_track)
        {
            energy_density -= fraction * energy_density;
            pflow_remnant -= fraction * pflow_remnant;
        }
        else if (is_1st_local_max_assoisited_with_track && !is_2nd_local_max_assoisited_with_track)
        {
            energy_density -= fraction * energy_density * first_local_max_weight;
            pflow_remnant -= fraction * pflow_remnant * first_local_max_weight;
            if ((1 - fraction * first_local_max_weight) != 0)
            {
                second_local_max_weight *= 1 / (1 - fraction * first_local_max_weight);
            }
            else
            {
                second_local_max_weight = 0;
            }
            first_local_max_weight = 1 - second_local_max_weight;
        }
        else if (!is_1st_local_max_assoisited_with_track && is_2nd_local_max_assoisited_with_track)
        {
            energy_density -= fraction * energy_density * second_local_max_weight;
            pflow_remnant -= fraction * pflow_remnant * second_local_max_weight;
            second_local_max_weight -= fraction * second_local_max_weight;
            if ((1 - fraction * second_local_max_weight) != 0)
            {
                first_local_max_weight *= 1 / (1 - fraction * second_local_max_weight);
            }
            else
            {
                first_local_max_weight = 0;
            }
            second_local_max_weight = 1 - second_local_max_weight;
        }
    }
    else
    {
        energy_density -= fraction * energy_density;
        pflow_remnant -= fraction * pflow_remnant;
    }
}

int Cell::modulo(int value, int mo)
{
    int mod = value % (int)mo;
    if (value < 0)
    {
        mod += mo;
    }
    return mod;
}

float Cell::get_signal_to_noise()
{
    return sigma;
}

int Cell::get_layer()
{
    return layer_idx;
}

int Cell::get_eta()
{
    return eta_idx;
}

int Cell::get_phi()
{
    return phi_idx;
}

float Cell::get_1st_local_max_weight()
{
    return first_local_max_weight;
}
float Cell::get_2nd_local_max_weight()
{
    return second_local_max_weight;
}

float Cell::get_abs_min_eta()
{
    return abs_eta_min;
}

float Cell::get_abs_max_eta()
{
    return abs_eta_max;
}

float Cell::get_abs_min_phi()
{
    return abs_phi_min;
}
float Cell::get_abs_max_phi()
{
    return abs_phi_max;
}

void Cell::set_1st_local_max_label(int LocalLabel)
{
    first_local_max_label = LocalLabel;
}
int Cell::get_1st_local_max_label()
{
    return first_local_max_label;
}
void Cell::set_2nd_local_max_label(int LocalLabel)
{
    second_local_max_label = LocalLabel;
}
int Cell::get_2nd_local_max_label()
{
    return second_local_max_label;
}
void Cell::set_1st_local_max_indexes(int lay, int eta, int phi)
{
    _1st_local_max_indexes.clear();
    _1st_local_max_indexes.push_back(lay);
    _1st_local_max_indexes.push_back(eta);
    _1st_local_max_indexes.push_back(phi);
}
// void Cell::set_1st_local_max_indexes(std::vector<int> &coord)
// {
//     _1st_local_max_indexes = coord;
// }
void Cell::get_1st_local_max_indexes(std::vector<int> &coord)
{
    coord = _1st_local_max_indexes;
}
void Cell::set_2nd_local_max_indexes(int lay, int eta, int phi)
{
    _2nd_local_max_indexes.clear();
    _2nd_local_max_indexes.push_back(lay);
    _2nd_local_max_indexes.push_back(eta);
    _2nd_local_max_indexes.push_back(phi);
}
// void Cell::set_2nd_local_max_indexes(std::vector<int> &coord)
// {
//     _2nd_local_max_indexes = coord;
// }
void Cell::get_2nd_local_max_indexes(std::vector<int> &coord)
{
    coord = _2nd_local_max_indexes;
}

void Cell::set_noise_signal(float en)
{
    total_energy -= noise;
    noise = en;
    total_energy += noise;
    pflow_remnant = total_energy;
}
void Cell::add_charge_energy(float en)
{
    charge_energy += en;
    total_energy += en;
    pflow_remnant = total_energy;
}
float Cell::get_charge_energy()
{
    return charge_energy;
}
void Cell::add_neutral_energy(float en)
{
    neutral_energy += en;
    total_energy += en;
    pflow_remnant = total_energy;
}
float Cell::get_neutral_energy()
{
    return neutral_energy;
}
void Cell::set_signal_to_noise(float LocalSigma)
{
    sigma = LocalSigma;
}

void Cell::add_particle(Particle_dep_in_cell p, bool isConversionElectron)
{

    std::vector<Particle_dep_in_cell> *particle_list_target = &list_of_particles_dep_energy;
    if ( isConversionElectron ) {
	particle_list_target = &list_of_conv_electrons_dep_energy;
    }
    
    int sizeptr = particle_list_target->size();
    bool check = true;
    for (int part = 0; part < sizeptr; part++)
    {
        if (particle_list_target->at(part).particle_pos_in_true_list == p.particle_pos_in_true_list)
        {
            particle_list_target->at(part).Energy += p.Energy;
            check = false;
        }
    }
    if (check)
    {
        particle_list_target->push_back(p);
    }
}

void Cell::get_particles(std::vector<Particle_dep_in_cell> &particles)
{
    particles = list_of_particles_dep_energy;
}

void Cell::get_conv_electrons(std::vector<Particle_dep_in_cell> &conv_electrons)
{
    conv_electrons = list_of_conv_electrons_dep_energy;
}

void Cell::set_indexes(int lay, int ind_e, int ind_p, bool is_high)
{
    Config_reader_var &config_var = Config_reader_var::GetInstance();
    float deta;
    float dphi;
    float R_inn_0, R_mid_central;
        
    if (is_high)
    {
	R_inn_0       = config_var.high_resolution.layer_inn_radius_flatten[0];
        deta          = config_var.high_resolution.layer_deta_flatten.at(lay);
        dphi          = config_var.high_resolution.layer_dphi_flatten.at(lay);
        R_mid_central = config_var.high_resolution.layer_mid_radius_flatten.at(lay);
    }
    else
    {
	R_inn_0       = config_var.low_resolution.layer_inn_radius_flatten[0];
        deta          = config_var.low_resolution.layer_deta_flatten.at(lay);
        dphi          = config_var.low_resolution.layer_dphi_flatten.at(lay);
        R_mid_central = config_var.low_resolution.layer_mid_radius_flatten.at(lay);
    }

    layer_idx = lay;
    eta_idx = ind_e;
    phi_idx = ind_p;

    eta_pos = -1 * config_var.max_eta_endcap + ((float(eta_idx) + 0.5) * deta);   

    abs_eta_min = -1 * config_var.max_eta_endcap + deta * (float)eta_idx;
    abs_eta_max = -1 * config_var.max_eta_endcap + deta * ((float)eta_idx + 1);

    theta_pos = 2.0 * TMath::ATan(TMath::Exp(-eta_pos));

    phi_pos = ((float(phi_idx) + 0.5) * dphi);
    if (phi_pos > TMath::Pi())
        phi_pos -= 2.0 * TMath::Pi();
    //!
    abs_phi_min = dphi * phi_idx;
    abs_phi_max = dphi * (phi_idx + 1);

    float depth = R_mid_central - R_inn_0;
    
    if ( fabs( eta_pos ) <= config_var.max_eta_barrel ) {
	
	//
	// Barrel cell
	//
	
	float R = R_inn_0 + depth / cosh( eta_pos );
	
        x_pos = R * cos(phi_pos);
        y_pos = R * sin(phi_pos);
        z_pos = R / tan(theta_pos);
	
    } else {
	
	//
	// Endcap cell
	//

	float theta_endcap_max = EtaToTheta( config_var.max_eta_barrel );
	float dz_ip_endcap     = R_inn_0 / tan( theta_endcap_max );
	float d_ip             = dz_ip_endcap / fabs( cos( theta_pos ) ) + depth; // distance between the cell and the IP
	
	x_pos = d_ip * sin( theta_pos ) * cos( phi_pos );
	y_pos = d_ip * sin( theta_pos ) * sin( phi_pos );
	z_pos = d_ip * cos( theta_pos );

    }

}

float Cell::get_noise_signal()
{
    return noise;
}

void Cell::set_label(int LocalLabel)
{
    label = LocalLabel;
}
int Cell::get_label()
{
    return label;
}

void Cell::set_is_cell_shared(bool shar)
{
    is_shared = shar;
}
bool Cell::get_is_cell_shared()
{
    return is_shared;
}
float Cell::get_total_energy()
{
    return total_energy;
}

float Cell::get_pflow_remnant()
{
    return pflow_remnant;
}
void Cell::set_pflow_remnant(float pflowremnant)
{
    pflow_remnant = pflowremnant;
}

void Cell::set_distance_to_track(float R)
{
    distance_to_track = R;
}

float Cell::get_distance_to_track()
{
    return distance_to_track;
}

void Cell::set_energy_dencity(float enden)
{
    energy_density = enden;
}

float Cell::get_energy_dencity()
{
    return energy_density;
}

void Cell::set_number_of_ring(int num)
{
    ring_number = num;
}

float Cell::get_number_of_ring()
{
    return ring_number;
}

void Cell::set_is_1st_local_max_assoisited_with_track(bool as)
{
    is_1st_local_max_assoisited_with_track = as;
}

bool Cell::get_is_1st_local_max_assoisited_with_track()
{
    return is_1st_local_max_assoisited_with_track;
}

void Cell::set_is_2nd_local_max_assoisited_with_track(bool as)
{
    is_2nd_local_max_assoisited_with_track = as;
}

bool Cell::get_is_2nd_local_max_assoisited_with_track()
{
    return is_2nd_local_max_assoisited_with_track;
}

void Cell::set_truth_label()
{
    int n_particle_size = list_of_particles_dep_energy.size();

    float max_parE = -500;
    int highest_par_label = -1;
    for (int ip = 0; ip < n_particle_size; ip++)
    {
        if (list_of_particles_dep_energy.at(ip).Energy > max_parE)
        {
            max_parE = list_of_particles_dep_energy.at(ip).Energy;
            highest_par_label = list_of_particles_dep_energy.at(ip).PDG_ID;
        }
    }

    truth_label = highest_par_label;
}

int Cell::get_truth_label()
{
    return truth_label;
}

float Cell::get_x()
{
    return x_pos;
}
float Cell::get_y()
{
    return y_pos;
}
float Cell::get_z()
{
    return z_pos;
}
float Cell::get_eta_pos()
{
    return eta_pos;
}
float Cell::get_phi_pos()
{
    return phi_pos;
}

void Cell::set_size_of_seeds_list(int size)
{
    size_of_seeds_list = size;
}

int Cell::get_size_of_seeds_list()
{
    return size_of_seeds_list;
}

void Cell::set_cell_link_superres(std::vector<std::vector<int>> cell_link)
{
    cell_link_superres = cell_link;
}

void Cell::get_cell_link_superres(std::vector<std::vector<int>> &cell_link)
{
    cell_link = cell_link_superres;
}

void Cell::set_local_max_weights(std::vector<float> & _1st_local_max_cluster, float energy_1st_local_max_cluster ,std::vector<float> &_2nd_local_max_cluster, float energy_2nd_local_max_cluster)
{
    
    float distance_to_1st = sqrtf(sqr(x_pos - _1st_local_max_cluster.at(0)) +
                                  sqr(y_pos - _1st_local_max_cluster.at(1)) +
                                  sqr(z_pos - _1st_local_max_cluster.at(2)));
    
    float distance_to_2nd = sqrtf(sqr(x_pos - _2nd_local_max_cluster.at(0)) +
                                  sqr(y_pos - _2nd_local_max_cluster.at(1)) +
                                  sqr(z_pos - _2nd_local_max_cluster.at(2)));
    float exp_r = TMath::Exp(distance_to_1st - distance_to_2nd);
    first_local_max_weight = energy_1st_local_max_cluster/(energy_1st_local_max_cluster + exp_r * energy_2nd_local_max_cluster);
    second_local_max_weight = 1 - first_local_max_weight;

}

