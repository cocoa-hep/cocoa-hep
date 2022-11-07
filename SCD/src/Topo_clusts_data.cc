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
}

void Topo_clusts_data::clear()
{
    topo_clusts_list.clear();
    number_of_cluster = 0;
    jets_objects.clear();
}
