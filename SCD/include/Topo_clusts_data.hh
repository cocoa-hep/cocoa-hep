#ifndef __TOPO_CLUSTS_DATA_H__
#define __TOPO_CLUSTS_DATA_H__

#include "Topo_clust_var.hh"
#include "TTree.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

class Topo_clusts_data
{
public:
    Topo_clusts_data(std::string prefix = "");
    std::vector<Topo_clust> topo_clusts_list;
    std::vector<float> topo_idx;
    std::vector<float> topo_bary_eta;
    std::vector<float> topo_bary_phi;
    std::vector<float> topo_bary_R;
    std::vector<float> topo_bary_sigma_eta;
    std::vector<float> topo_bary_sigma_phi;
    std::vector<float> topo_e;
    void set_tree_branches(TTree *outTree);
    void fill_topo_var();
    void make_pseudo_jet_particles();
    static Topo_clusts_data &GetInstance()
    {
        static Topo_clusts_data topo_list;
        return topo_list;
    };
    static Topo_clusts_data &GetInstance_neutral()
    {
        static Topo_clusts_data topo_list("neutral");
        return topo_list;
    };
    static Topo_clusts_data &GetInstance_charge()
    {
        static Topo_clusts_data topo_list("charge");
        return topo_list;
    };
    static Topo_clusts_data &GetInstance_noise()
    {
        static Topo_clusts_data topo_list("noise");
        return topo_list;
    };
    void clear();
    std::vector <fastjet::PseudoJet> jets_objects;
private:
    int number_of_cluster;
    std::string Prefix;
};

#endif // __TOPO_CLUSTS_DATA_H__