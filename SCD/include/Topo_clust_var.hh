#ifndef __TOPO_CLUST_VAR_H__
#define __TOPO_CLUST_VAR_H__

#include "Cell_var.hh"
#include "Config_reader_var.hh"
#include <vector>

class Topo_clust
{
public:
    Topo_clust();
    Topo_clust(int label, Geometry_definition Geometry);
    // Topo_clust(const Topo_clust &orig) = default;
    int label;
    float total_energy;
    float EM_energy;
    float abs_total_energy;
    float noise;
    float neutral_energy;
    float charge_energy;
    float eta_com;
    float phi_com;
    float R_com;
    float x_com;
    float y_com;
    float z_com;
    std::vector<float> xyz_com;
    float sigma_eta;
    float sigma_phi;
    float px;
    float py;
    float pz;
    int truth_link;
    float cell_energy;
    std::vector<std::pair<float,int>> closest_tracks;
    int charge;

    void add_cell(Cell &cell);
    void subtract_cell(Cell &cell, float fraction);

private:
    Geometry_definition geometry;
    float cos_phi;
    float sin_phi;
    float sigma_eta_buf;
    float sigma_phi_buf;
    int n_cell;
    void COM_update(Cell &cell, float weight, int add_subtract);
    void cell_update(Cell &cell, int add_subtract, float fraction);
    std::vector<float> Energy_Share = {0, 0, 0, 0, 0};
    std::map<int, int> particle_type = {
        {22, 0},
        {11, 1},
        {-11, 1},
        {13, 2},
        {-13, 2},
        {130, 3},
        {310, 3},
        {2112, 3},
        {-2112, 3},
        {-3322, 3},
        {3322, 3},
        {3122, 3},
        {-3122, 3},
        {211, 4},
        {-211, 4},
        {321, 4},
        {-321, 4},
        {2212, 4},
        {-2212, 4},
        {3112, 4},
        {-3112, 4},
        {3312, 4},
        {-3312, 4},
        {3222, 4},
        {-3222, 4},
        {3334, 4},
        {-3334, 4},
    };
};
#endif // __TOPO_CLUST_VAR_H__
