#ifndef __SUPERCLUSTERING_DATA_H__
#define __SUPERCLUSTERING_DATA_H__

#include "Superclustering.hh"
#include "Supercluster.hh"
#include "TTree.h"

class Superclustering_data
{
    private:

    public:
    Superclustering_data();
    ~Superclustering_data();
    static Superclustering_data &GetInstance()
    {
        static Superclustering_data instance;
        return instance;
    };
    void set_tree_branches(TTree *outTree);
    void fill_supercluster_data();
    void clear();
    std::vector<Supercluster> super_list;
    std::vector<float> super_e;
    std::vector<float> super_eta;
    std::vector<float> super_phi;
    std::vector<float> seed_e;
    std::vector<float> track_eta;
    std::vector<float> track_phi;
    std::vector<float> track_pt;
    std::vector<int> super_pdgid;
    std::vector<int> super_N;
    std::vector<std::vector<float>> topo_idx_list;
    std::vector<int> track_idx;
    std::vector<int> conv_track_idx;
};

#endif