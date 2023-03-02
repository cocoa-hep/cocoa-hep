#ifndef __PARTICLE_FLOW_FUNC_H__
#define __PARTICLE_FLOW_FUNC_H__

#include "vector"
#include "algorithm"
#include "Track_var.hh"
#include "Topo_clust_var.hh"
#include "Cell_var.hh"
#include "Config_reader_var.hh"
#include "CSVReader.hh"
#include "Particle_flow_data.hh"
#include "Particle_flow_var.hh"
#include "Debug_Particle_Flow_func.hh"
#include <tuple>

struct Ring_substraction
{
    int ring_number = 0;
    int layer = 0;
    float energy = 0;
    float energy_dencity = 0;
};

inline bool sort_by_energy_dencity(const Ring_substraction &vect_1, const Ring_substraction &vect_2)
{
    return vect_1.energy_dencity > vect_2.energy_dencity;
}
// inline bool sort_by_pt(const Track_struct &track_1, const Track_struct &track_2)
// {
//     return track_1.pt > track_2.pt;
// }

class Particle_flow_func
{
private:
    float sum_e_dep = 0;
    float E_div_p_template;
    float sigma_E_div_p_template;
    float delta_Rprime_threshold;
    
    int layer_max;
    std::string Type_Of_Runnig;

    std::vector<Cell *> cells_in_topoclust;
    std::vector<Topo_clust> Topo_List;
    Geometry_definition geometry;
    Particle_flow_alg_var pflow_param;
public:
    void track_match_to_topoclusters(Track_struct &track);
    void recovering_shower_split(Track_struct &track);
    void cell_by_cell_subst(Track_struct &Track);
    static float cell_volume(Cell &cell, Geometry_definition geometry_);
    static std::tuple<float, float> energy_div_momentum_template(float pt, int lhed, float eta);
    void subst_from_cell(Cell &cell, float fraction, std::vector<Particle_flow_var> pflow_objs);
    void modified_topoclusters(std::vector<Particle_flow_var> &topoclust_pflow_obj);
    static float phi_distance(float phi_1, float phi_2);
    static float R_distance_func(float phi_1, float eta_1, float phi_2, float eta_2, float sigma_phi = 1, float sigma_eta = 1);

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
        {211, 4},
        {-211, 4},
        {321, 4},
        {-321, 4},
        {2212, 4},
        {-2212, 4}};
    Particle_flow_func(std::vector<Track_struct> &track_list, std::vector<Topo_clust> &topo_list, std::vector<Cell *> &calo_cells,
                       std::vector<Particle_flow_var> &pflow_list, Geometry_definition Geometry, Particle_flow_alg_var Pflow_Param);
    template<typename Charge_particle> static inline int LHED(Charge_particle track, Geometry_definition geometry_, std::vector<Cell *> &cells_in_topoclust_)
    {
        std::vector<double> eta = track.eta;
        std::vector<double> phi = track.phi;
        Config_reader_var &config_var = Config_reader_var::GetInstance();
        float constant = config_var.particle_flow.Moliere_radius;
        std::vector<float> energy_dencity(geometry_.kNLayers - 1, 0);
        int num_ECAL_layers = geometry_.layer_deta_ECAL.size();
        std::vector<float> deepth(num_ECAL_layers + 1, 0);
        for (int ideep = 0; ideep < num_ECAL_layers; ideep++)
        {
            if (ideep != 0)
                deepth.at(ideep) = deepth.at(ideep - 1);
            deepth.at(ideep) += (geometry_.layer_out_radius_ECAL.at(ideep).at(0) - geometry_.layer_inn_radius_ECAL.at(ideep).at(0)) / config_var.Material_ECAL->GetNuclearInterLength();
        }
        deepth.back() = deepth.at(num_ECAL_layers - 1);
        for (int ihlay = 0; ihlay < (int)(geometry_.resolution_width_of_HCAL_layers_in_Lambda_int.size()); ihlay++)
        {
            deepth.back() += geometry_.resolution_width_of_HCAL_layers_in_Lambda_int.at(ihlay).at(0);
        }
        int size_cells_in_topoclust_ = cells_in_topoclust_.size();
        for (int icell = 0; icell < size_cells_in_topoclust_; icell++)
        {
            Cell *cell = cells_in_topoclust_.at(icell);
            if (cell->get_label() == track.index_of_closest_topoclusters.front())
            {
                int idx = cell->get_layer();
                float deta_min = fabs(cell->get_abs_min_eta() - eta.at(idx));
                float deta_max = fabs(cell->get_abs_max_eta() - eta.at(idx));
                float dphi_min = fabs(cell->get_abs_min_phi() - phi.at(idx));
                float dphi_max = fabs(cell->get_abs_max_phi() - phi.at(idx));
                if (dphi_min > M_PI)
                    dphi_min -= 2 * M_PI;
                if (dphi_max > M_PI)
                    dphi_max -= 2 * M_PI;
                float lhed = fabs(cell->get_total_energy()) / cell_volume(*cell, geometry_);
                float weight_var = 0.25 * fabs(erf(deta_min / (sqrt(2) * constant)) - erf(deta_max / (sqrt(2) * constant))) *
                                fabs(erf(dphi_min / (sqrt(2) * constant)) - erf(dphi_max / (sqrt(2) * constant)));
                cell->set_energy_dencity(lhed * weight_var);

                if (idx >= num_ECAL_layers)
                    idx = num_ECAL_layers;

                energy_dencity.at(idx + 1) += cell->get_energy_dencity();
            }
        }
        std::vector<float> rate_of_increase;
        for (int i = 0; i < num_ECAL_layers + 1; i++)
            rate_of_increase.push_back((energy_dencity.at(i + 1) - energy_dencity.at(i)) / deepth.at(i));
        float max_buf = 0.0;
        int _layer_max = -1;
        for (int i = 0; i < num_ECAL_layers + 1; i++)
        {
            if (rate_of_increase.at(i) >= max_buf)
            {
                max_buf = rate_of_increase.at(i);
                _layer_max = i;
            }
        }
        return _layer_max;
    }

};

#endif // __PARTICLE_FLOW_FUNC_H__