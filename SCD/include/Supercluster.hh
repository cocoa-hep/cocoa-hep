#ifndef __SUPERCLUSTER_H__
#define __SUPERCLUSTER_H__


#include <vector>
#include "Topo_clust_var.hh"
#include "Topo_clust_func.hh"
#include "Track_var.hh"

inline bool compare_e(const Topo_clust &topo_1, const Topo_clust &topo_2)
{
    return topo_1.total_energy > topo_2.total_energy; //TODO: use only energy in the ECAL
}

class Supercluster
{
    private:
        Track_struct track_assoc;
        Topo_clust topo_seed;
        std::vector<Topo_clust> topo_members;
        void sort_by_energy();
        void update();
    public:
        int n_clusters;
        float total_energy;
        Supercluster();
        void add_cluster(Topo_clust topo);
        void set_seed(Topo_clust seed);
        void set_track(Track_struct track);
        Topo_clust get_seed();
};

#endif //__SUPERCLUSTER_H__