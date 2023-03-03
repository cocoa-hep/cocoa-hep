#ifndef __PION_INFO_FOR_PFLOW_H__
#define __PION_INFO_FOR_PFLOW_H__

#include "TTree.h"
#include <vector>
#include "Track_var.hh"
#include "Topo_clust_var.hh"
#include "Cell_var.hh"
#include "Pion_for_pflow_var.hh"
#include "Debug_Particle_Flow_data.hh"
#include "Config_reader_var.hh"
#include "Particle_flow_func.hh"

class Debug_Particle_Flow_func
{
private:
    template <typename T>
    int sgn(T val);
    // vector <vector <float> > ChPionsInfo;

public:
    Debug_Particle_Flow_func(std::vector<Track_struct> track_list,
                        std::vector<Topo_clust> &topo_list,
                        std::vector<Cell *> &topo_cells,
                        Debug_Particle_Flow_data &debug_data,
                        Geometry_definition Geometry,
                        Particle_flow_alg_var Pflow_Param,
                        std::string type_of_running);
    void Containment_shower(std::vector<Track_struct> &track_list,
                            Debug_Particle_Flow_data &debug_data);
    void PFlow_parameters(Debug_Particle_Flow_data &debug_data,
                          std::string type_of_running);
    void debug_track_match_to_topoclusters(Pion_for_pflow_var &pion);
    void debug_parameters(Pion_for_pflow_var &pion);
    void clear();
    Geometry_definition geometry;
    float E_div_p_template;
    float sigma_E_div_p_template;
    float delta_Rprime_threshold;
    int layer_max;
    std::string Type_Of_Runnig;
    std::vector<Topo_clust> Topo_List;
    Config_reader_var &config_var = Config_reader_var::GetInstance();
    Particle_flow_alg_var pflow_param;
    std::vector<Cell *> &cells_in_topoclust;
    std::vector<Topo_clust> &topoclusters_list;
};

#endif // __PION_INFO_FOR_PFLOW_H__