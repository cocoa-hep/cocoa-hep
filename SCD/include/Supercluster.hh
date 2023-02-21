#ifndef __SUPERCLUSTER_H__
#define __SUPERCLUSTER_H__


#include <vector>
#include "Topo_clust_var.hh"
#include "Topo_clust_func.hh"
#include "Track_var.hh"


class ConversionVertex
{
    private:
        std::vector<int> track_list;
        std::vector<int> topo_list;
        float x;
        float y;
        float z;
    public:
        ConversionVertex();
        bool try_add_track(Track_struct);
        bool try_add_topo(Topo_clust);
        bool contains(Track_struct);
        bool contains(Topo_clust);
        int get_friend(Topo_clust);
        int get_friend(Track_struct);
};


inline bool compare_e(const Topo_clust &topo_1, const Topo_clust &topo_2)
{
    return topo_1.EM_energy > topo_2.EM_energy;
}

class Supercluster
{
    private:
        Track_struct track_assoc;
        Track_struct track_conv;
        Topo_clust topo_seed;
        ConversionVertex conv_vertex;
        std::vector<Topo_clust> topo_members;
        void sort_by_energy();
        void update();
    public:
        int n_clusters;
        float total_energy;
        float distance_seed_track;
        std::vector<float> com;
        float eta;
        float phi;
        Supercluster();
        void add_cluster(Topo_clust topo);
        void set_seed(Topo_clust seed);
        void set_track(Track_struct track);
        void set_conv_track(Track_struct track);
        void set_conv_vertex(ConversionVertex vtx);
        std::vector<Topo_clust> get_clusters();
        Topo_clust get_merged_cluster();
        Track_struct get_track();
        Track_struct get_conv_track();
        ConversionVertex get_conv_vertex();
};

#endif //__SUPERCLUSTER_H__