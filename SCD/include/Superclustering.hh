#ifndef __SUPERCLUSTERING_H__
#define __SUPERCLUSTERING_H__


#include <vector>
#include "Config_reader_var.hh"
#include "Cell_var.hh"
#include "Cells_data.hh"
#include "Track_var.hh"
#include "Tracks_data.hh"
#include "Topo_clust_func.hh"
#include "Topo_clust_var.hh"
#include "Topo_clusts_data.hh"
#include "Supercluster.hh"

inline bool compare_super_e(const Supercluster &super_1, const Supercluster &super_2)
{
    return super_1.total_energy > super_2.total_energy;
}

class Superclustering
{
    private:        
        std::vector<Cell *> cells_in_topoclust;
        std::vector<Topo_clust> Topo_List;
        std::vector<Track_struct> Track_list;
        //std::vector<std::pair<Topo_clust,Track_struct>> seed_track_pairs;
        std::vector<Supercluster> Super_list;
        Geometry_definition geometry;
        void sort_by_energy();
    public:
        Superclustering(std::vector<Track_struct> &_track_list, std::vector<Topo_clust> &_topo_list, std::vector<Cell *> &_cell_list, Geometry_definition Geometry);
        void find_seed_clusters();
        void add_neighbor_clusters();
};

#endif //__SUPERCLUSTERING_H__