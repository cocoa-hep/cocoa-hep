#include "Superclustering_data.hh"
#include "Superclustering.hh"
#include "Supercluster.hh"


Superclustering_data::Superclustering_data(){}

Superclustering_data::~Superclustering_data(){}

void Superclustering_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch("supercluster_e",     "vector<vector<float>>", &super_e);
    outTree->Branch("supercluster_seed",  "vector<vector<int>>",   &seed_idx);
    outTree->Branch("supercluster_track", "vector<vector<int>>",   &track_idx);
}

void Superclustering_data::fill_supercluster_data()
{
    int n_superclusters = super_list.size();
    for(int isuper=0; isuper<n_superclusters; isuper++)
    {
        super_e.push_back(super_list.at(isuper).total_energy);
        seed_idx.push_back(super_list.at(isuper).get_seed().label - 1);
        track_idx.push_back(super_list.at(isuper).get_track().position_in_list);
    }
}

void Superclustering_data::clear()
{
    super_list.clear();
    super_e.clear();
    seed_idx.clear();
    track_idx.clear();
}