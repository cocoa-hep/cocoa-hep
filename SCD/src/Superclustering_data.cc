#include "Superclustering_data.hh"
#include "Superclustering.hh"
#include "Supercluster.hh"


Superclustering_data::Superclustering_data(){}

Superclustering_data::~Superclustering_data(){}

void Superclustering_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch("supercluster_e",     "vector<float>",         &super_e);
    outTree->Branch("supercluster_seed_e","vector<float>",         &seed_e);
    outTree->Branch("supercluster_topos", "vector<vector<float>>", &topo_idx_list);
    outTree->Branch("supercluster_track", "vector<int>",           &track_idx);
    outTree->Branch("supercluster_conv_track", "vector<int>",      &conv_track_idx);
    outTree->Branch("supercluster_pdgid", "vector<int>",           &super_pdgid);
    outTree->Branch("supercluster_N",     "vector<int>",           &super_N);
}

void Superclustering_data::fill_supercluster_data()
{
    int n_superclusters = super_list.size();

    //Find out total size needed
    int needed_size = 0;
    for(int isuper=0; isuper<n_superclusters; isuper++)
        needed_size += super_list.at(isuper).get_clusters().size();

    topo_idx_list.resize(needed_size);
    for(int isuper=0; isuper<n_superclusters; isuper++)
    {
        super_e.push_back(super_list.at(isuper).total_energy);
        seed_e.push_back(super_list.at(isuper).get_clusters()[0].total_energy);

        std::vector<Topo_clust> topos = super_list.at(isuper).get_clusters();
        const int n_topos = topos.size();
        for(int itopo=0; itopo<topos.size(); itopo++){
            topo_idx_list.at(isuper).push_back(topos.at(itopo).label - 1);
        }

        track_idx.push_back(super_list.at(isuper).get_track().position_in_list);
        conv_track_idx.push_back(super_list.at(isuper).get_conv_track().position_in_list);
        super_pdgid.push_back(super_list.at(isuper).get_track().pdgcode);
        super_N.push_back(n_topos);
    }
}

void Superclustering_data::clear()
{
    super_list.clear();
    super_e.clear();
    seed_e.clear();
    super_pdgid.clear();
    super_N.clear();
    topo_idx_list.clear();
    track_idx.clear();
    conv_track_idx.clear();
}