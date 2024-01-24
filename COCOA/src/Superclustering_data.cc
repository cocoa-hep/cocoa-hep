#include "Superclustering_data.hh"
#include "Superclustering.hh"
#include "Supercluster.hh"


Superclustering_data::Superclustering_data(){}

Superclustering_data::~Superclustering_data(){}

void Superclustering_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch("supercluster_e",          "vector<float>",         &super_e);
    outTree->Branch("supercluster_phi",        "vector<float>",         &super_phi);
    outTree->Branch("supercluster_eta",        "vector<float>",         &super_eta);
    outTree->Branch("supercluster_seed_e",     "vector<float>",         &seed_e);
    outTree->Branch("supercluster_track_phi",  "vector<float>",         &track_phi);
    outTree->Branch("supercluster_track_eta",  "vector<float>",         &track_eta);
    outTree->Branch("supercluster_track_pt",   "vector<float>",         &track_pt);
    outTree->Branch("supercluster_topos",      "vector<vector<float>>", &topo_idx_list);
    outTree->Branch("supercluster_track",      "vector<int>",           &track_idx);
    outTree->Branch("supercluster_conv_track", "vector<int>",           &conv_track_idx);
    outTree->Branch("supercluster_pdgid",      "vector<int>",           &super_pdgid);
    outTree->Branch("supercluster_N",          "vector<int>",           &super_N);
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
        super_phi.push_back(super_list.at(isuper).phi);
        super_eta.push_back(super_list.at(isuper).eta);
        seed_e.push_back(super_list.at(isuper).get_clusters()[0].total_energy);

        std::vector<Topo_clust> topos = super_list.at(isuper).get_clusters();
        const int n_topos = topos.size();
        for(size_t itopo=0; itopo<topos.size(); itopo++){
            topo_idx_list.at(isuper).push_back(topos.at(itopo).label - 1);
        }

        Track_struct assoc_track = super_list.at(isuper).get_track();
        float track_phi_i   = atan2(assoc_track.py,assoc_track.px);
        float track_p_i     = sqrtf(assoc_track.px*assoc_track.px + assoc_track.py*assoc_track.py + assoc_track.pz*assoc_track.pz);
        float track_theta_i = acos(assoc_track.pz/track_p_i);
        float track_eta_i   = -1*log(tan(track_theta_i/2.));

        track_phi.push_back(track_phi_i);
        track_eta.push_back(track_eta_i);
        track_pt.push_back(super_list.at(isuper).get_track().pt);
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
    super_phi.clear();
    super_eta.clear();
    seed_e.clear();
    track_eta.clear();
    track_phi.clear();
    track_pt.clear();
    super_pdgid.clear();
    super_N.clear();
    topo_idx_list.clear();
    track_idx.clear();
    conv_track_idx.clear();
}
