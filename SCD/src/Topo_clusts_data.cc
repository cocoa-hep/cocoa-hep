#include "Topo_clusts_data.hh"

Topo_clusts_data::Topo_clusts_data(std::string prefix)
{
    Prefix = prefix;
    topo_clusts_list.clear();
}

void Topo_clusts_data::set_tree_branches(TTree *outTree)
{
    if (Prefix!="")
    {
        outTree->Branch(TString(Prefix + "_number_of_clusters_idx"), &number_of_cluster,TString(Prefix + "_number_of_clusters_idx/I"));
    }
    else
    {
        outTree->Branch("topo_idx",             "vector<float>", &topo_idx);
        outTree->Branch("topo_bary_eta",        "vector<float>", &topo_bary_eta);
        outTree->Branch("topo_bary_phi",        "vector<float>", &topo_bary_phi);
        outTree->Branch("topo_bary_rho",        "vector<float>", &topo_bary_R);
        outTree->Branch("topo_bary_sigma_eta",  "vector<float>", &topo_bary_sigma_eta);
        outTree->Branch("topo_bary_sigma_phi",  "vector<float>", &topo_bary_sigma_phi);
        outTree->Branch("topo_e",               "vector<float>", &topo_e);
    }
}

void Topo_clusts_data::make_pseudo_jet_particles()
{
    int size_topo_clusts_list = topo_clusts_list.size();
    for (int itopo = 0; itopo < size_topo_clusts_list; itopo++)
    {
        fastjet::PseudoJet topoclust(topo_clusts_list.at(itopo).px,
                                     topo_clusts_list.at(itopo).py,
                                     topo_clusts_list.at(itopo).pz,
                                     topo_clusts_list.at(itopo).total_energy);
        topoclust.set_user_index(itopo);
        jets_objects.push_back(topoclust);
    }
}

void Topo_clusts_data::fill_topo_var()
{
    number_of_cluster = topo_clusts_list.size();
    
    for (int itopo = 0; itopo < number_of_cluster; itopo++)
    {
        topo_idx.push_back(topo_clusts_list.at(itopo).label);
        topo_bary_eta.push_back(topo_clusts_list.at(itopo).eta_com);
        topo_bary_phi.push_back(topo_clusts_list.at(itopo).phi_com);
        topo_bary_R.push_back(topo_clusts_list.at(itopo).R_com);
        topo_bary_sigma_eta.push_back(topo_clusts_list.at(itopo).sigma_eta);
        topo_bary_sigma_phi.push_back(topo_clusts_list.at(itopo).sigma_phi);
        topo_e.push_back(topo_clusts_list.at(itopo).total_energy);
    }
}

void Topo_clusts_data::clear()
{
    topo_clusts_list.clear();
    number_of_cluster = 0;
    jets_objects.clear();

    topo_idx.clear();
    topo_bary_eta.clear();
    topo_bary_phi.clear();
    topo_bary_R.clear();
    topo_bary_sigma_eta.clear();
    topo_bary_sigma_phi.clear();
    topo_e.clear();
}
