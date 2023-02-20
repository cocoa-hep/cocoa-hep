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

inline bool compare_topo_e(const Topo_clust &topo_1, const Topo_clust &topo_2)
{
    return topo_1.EM_energy > topo_2.EM_energy;
}

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
        std::vector<ConversionVertex> ConvVertex_list;
        void sort_by_energy(std::vector<Supercluster> &Super_list);
    public:
        Superclustering(const std::vector<Track_struct> &_track_list, const std::vector<Topo_clust> &_topo_list, const std::vector<Cell *> &_cell_list, std::vector<Supercluster> &Super_list);
        void topo_match_to_tracks();
        void find_conversion_vertices();
        void find_seed_clusters(std::vector<Supercluster> &Super_list);
        void add_neighbor_clusters(std::vector<Supercluster> &Super_list);
        void calc_merged_clusters(std::vector<Supercluster> &Super_list);
};

#endif //__SUPERCLUSTERING_H__
