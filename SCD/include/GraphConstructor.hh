#ifndef __GRAPHCONSTRUCTOR_H__
#define __GRAPHCONSTRUCTOR_H__


#include <vector>
#include <stdio.h>
#include <iostream>
#include <map>
#include <queue>
#include "TMatrixD.h"
#include "TVector3.h"
#include "Cells_data.hh"
#include "Tracks_data.hh"
#include "Full_trajectory_info_data.hh"
#include "Track_var.hh"
#include "Config_reader_var.hh"
#include "Graph_construction_data.hh"

class GraphConstructor
{
private:
    Config_reader_var& config_json_var = Config_reader_var::GetInstance();
public:
    GraphConstructor(std::vector<Cell*> &low_cells_in_topoclusters, std::vector<Cell*> &high_cells_in_topoclusters, std::vector<Track_struct> &tracks_list, std::vector<G4int> &particle_to_track, Graph_construction_data &graph_obj,
		     std::vector<float> *_particle_dep_energies = NULL);
    GraphConstructor(std::vector<Cell*> &low_cells_in_topoclusters, std::vector<Track_struct> &tracks_list, std::vector<G4int> &particle_to_track, Graph_construction_data &graph_obj, std::vector<float> *_particle_dep_energies = NULL);
    void fill_cell_to_cell_edges(std::vector<Cell*> &cell,
                                 std::vector<int> &cell_to_cell_edge_start, std::vector<int> &cell_to_cell_edge_end);
    void fill_track_to_cell_edges(std::vector<Cell*> &cell,
                                  std::vector<Track_struct> &track_vec, std::vector<int> &track_to_cell_edge_start, std::vector<int> &track_to_cell_edge_end);
    void fill_particle_to_node_edges(std::vector<Cell*> &cell,
                                     std::vector<int> particle_to_track, std::vector<std::vector<float>> &particle_to_node_idx, std::vector<std::vector<float>> &particle_to_node_weight,
				     std::vector<float> *_particle_dep_energies = NULL);
    void fill_superres_cell_to_cell_edges(std::vector<Cell*> &high_cells, std::vector<Cell*> &low_cells,
														std::vector<int> &high_cell_to_low_cell_edge, std::vector<int> &low_cell_to_high_cell_edge);
};



#endif // __GRAPHCONSTRUCTOR_H__
