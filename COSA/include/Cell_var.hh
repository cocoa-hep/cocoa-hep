#ifndef __CELL_VAR_H__
#define __CELL_VAR_H__

#include <cmath>
#include "TMath.h"
#include <algorithm>
#include <vector>
#include "Config_reader_var.hh"


struct Particle_dep_in_cell
{
    int PDG_ID;
    int particle_pos_in_true_list;
    float Energy;
};

class Cell
{
private:
    float total_energy;
    float charge_energy;
    float neutral_energy;
    float noise;
    float sigma;
    int layer_idx;
    int eta_idx;
    int phi_idx;
    int label;
    bool is_shared;
    float first_local_max_weight;
    float second_local_max_weight;
    int first_local_max_label;
    int second_local_max_label;
    int size_of_seeds_list;
    std::vector<std::vector<int>> cell_link_superres;
    std::vector<int> _1st_local_max_indexes;
    std::vector<int> _2nd_local_max_indexes;
    std::vector<Particle_dep_in_cell> list_of_particles_dep_energy;
    std::vector<Particle_dep_in_cell> list_of_conv_electrons_dep_energy;
    float eta_pos;
    float theta_pos;
    float phi_pos;
    float x_pos;
    float y_pos;
    float z_pos;
    float abs_eta_min;
    float abs_eta_max;
    float abs_phi_min;
    float abs_phi_max;

    float pflow_remnant;
    float distance_to_track;
    float energy_density;
    float ring_number;
    bool is_1st_local_max_assoisited_with_track;
    bool is_2nd_local_max_assoisited_with_track;

    int truth_label;

public:
    Cell();
    Cell(const Cell &orig);
    Cell(float ener, float chener, float nuener, float LocalSigma, int lay, int ind_e, int ind_p, int LocalLabel, bool shar);
    Cell *first_local_max;
    Cell *second_local_max;
    int position_in_list;

    void Reset();
    int modulo(int value, int m);

    void set_cell_link_superres(std::vector<std::vector<int>> cell_link);
    void get_cell_link_superres(std::vector<std::vector<int>> &cell_link);

    void set_truth_label();
    int get_truth_label();

    void set_label(int LocalLabel);
    int get_label();

    void set_is_cell_shared(bool shar);
    bool get_is_cell_shared();

    float get_total_energy();

    float get_charge_energy();
    void add_charge_energy(float en);

    float get_neutral_energy();
    void add_neutral_energy(float en);

    float get_noise_signal();
    void set_noise_signal(float en);

    void add_particle(Particle_dep_in_cell p, bool isConversionElectron = false );
    void get_particles(std::vector<Particle_dep_in_cell> &particles);
    void get_conv_electrons(std::vector<Particle_dep_in_cell> &conv_electrons);

    int getNParticles() { return list_of_particles_dep_energy.size(); }
    int getNConvElectrons() { return list_of_conv_electrons_dep_energy.size(); }

    float get_signal_to_noise();
    void set_signal_to_noise(float LocalSigma);

    void set_indexes(int l, int e, int phi, bool is_high);
    int get_layer();
    int get_eta();
    int get_phi();

    // void subtract_energys_fraction_from_pflow_remnant(float en, std::vector<Topo_clust> &TList);
    void subtract_energys_fraction_from_pflow_remnant(float en);
    float get_pflow_remnant();
    void set_pflow_remnant(float pflowremnant);

    void set_distance_to_track(float R);
    float get_distance_to_track();

    void set_energy_dencity(float enden);
    float get_energy_dencity();

    void set_number_of_ring(int num);
    float get_number_of_ring();

    void set_is_1st_local_max_assoisited_with_track(bool assos);
    bool get_is_1st_local_max_assoisited_with_track();

    void set_is_2nd_local_max_assoisited_with_track(bool assos);
    bool get_is_2nd_local_max_assoisited_with_track();

    float get_x();
    float get_y();
    float get_z();
    float get_eta_pos();
    float get_phi_pos();
    float get_theta_pos();

    void set_1st_local_max_label(int LocalLabel);
    int get_1st_local_max_label();
    void set_2nd_local_max_label(int LocalLabel);
    int get_2nd_local_max_label();

    void set_1st_local_max_indexes(int lay, int eta, int phi);
    // void set_1st_local_max_indexes(std::vector<int> &coord);
    void get_1st_local_max_indexes(std::vector<int> &coord);

    void set_2nd_local_max_indexes(int lay, int eta, int phi);
    // void set_2nd_local_max_indexes(std::vector<int> &coord);
    void get_2nd_local_max_indexes(std::vector<int> &coord);

    // void set_1st_local_max_weight(float FWeight);
    float get_1st_local_max_weight();

    // void set_2nd_local_max_weight(float SWeight);
    float get_2nd_local_max_weight();

    float get_abs_min_eta();
    float get_abs_max_eta();
    float get_abs_min_phi();
    float get_abs_max_phi();

    int get_size_of_seeds_list();
    void set_size_of_seeds_list(int size);
    void set_local_max_weights(std::vector<float> & _1st_local_max_cluster, float energy_1st_local_max_cluster ,std::vector<float> &_2nd_local_max_cluster, float energy_2nd_local_max_cluster);
};

#endif // __CELL_VAR_H__
