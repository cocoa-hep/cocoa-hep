#ifndef __TOPO_CLUST_FUNC_H__
#define __TOPO_CLUST_FUNC_H__

#include "Topo_clusts_data.hh"
#include "Topo_clust_var.hh"
#include "Cell_var.hh"
#include "Cells_data.hh"
#include <algorithm>

class Topo_clust_func
{

private:
    void find_seed_cells(std::vector<Cell*> &seeds);
    void cluster_maker(std::vector<Cell*> &seeds);
    void cluster_maker_neighbor(std::vector<Cell*> &seeds, int list_size);
    void cluster_maker_add_cell(Cell &cell, Cell &neighbor_cell, std::vector<Cell*> &seed);
    void global_label_change(int old_label, int new_label);
    void find_local_max(std::vector<Cell*> &local_maxs);
    void cluster_split(std::vector<Cell*> &local_maxs, std::vector<Cell*> &shares);
    void cluster_split_neighbor(std::vector<Cell*> &local_maxs, std::vector<Cell*> &shares, int clust_label, int list_size);
    void cluster_split_add_cell(Cell &cell, std::vector<Cell*> &local_max_cells_in_clust, std::vector<Cell*> &shares, Cell &neighbor_cell);
    void cluster_share_neighbor(std::vector<Cell*> &shares, int clust_label, int list_size);
    int number_of_labels();
    void fill_clusters_list(std::vector<Cell*> &shares, std::vector<Topo_clust> &topo_clusts_list);
    bool is_cell_local_max(Cell &cell);
    void insert_to_order_cell_list(Cell &cell, std::vector<Cell*> &cells_vector);
    void cluster_share_neighbor(Cell &cell, std::vector<Cell*> &shares, int clust_label);
    std::vector<std::vector<std::vector<Cell>>> &Cells_Array;
    int Nlayers;
    Geometry_definition geometry;
    Topological_clustering_var topo_config;

public:
    Topo_clust_func(std::vector<std::vector<std::vector<Cell>>> &cells_array, Geometry_definition Geometry, Topological_clustering_var Topo_Config, std::string class_of_debug = "Standard");
    void topoclustering(std::vector<Topo_clust> &topo_clusts_list);
};
#endif // __TOPO_CLUST_FUNC_H__