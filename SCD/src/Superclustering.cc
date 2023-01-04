#include "Superclustering.hh"
#include "Supercluster.hh"
#include "Topo_clust_var.hh"


bool in_eta_phi_window(const Topo_clust &topo_1, const Topo_clust &topo_2)
{

    float deta = fabs(topo_1.eta_com - topo_2.eta_com);
    float dphi = fabs(topo_1.phi_com - topo_2.phi_com);
    if (dphi > M_PI)
        dphi -= 2 * M_PI;

    bool passeta = (deta < 0.125);
    bool passphi = (dphi < 0.300);
 
    return (passeta && passphi);
}

Superclustering::Superclustering(std::vector<Track_struct> &_track_list, std::vector<Topo_clust> &_Topo_List, std::vector<Cell *> &_cell_list, Geometry_definition Geometry)
{
    cells_in_topoclust = _cell_list;
    Topo_List = _Topo_List;
    Track_list = _track_list;
    geometry = Geometry;

    //std::vector<Topo_clust> seeds;
    //seeds.reserve(50000000);
    find_seed_clusters();
    sort_by_energy();
    add_neighbor_clusters();

    for(int isuper=0; isuper<Super_list.size(); isuper++)
    {
        G4cout << "n_clusters   = " << Super_list.at(isuper).n_clusters   << G4endl;
        G4cout << "total_energy = " << Super_list.at(isuper).total_energy << G4endl;    
    }
}

void Superclustering::find_seed_clusters()
{
    int size_track_list = Track_list.size();
    for (int itrack = 0; itrack < size_track_list; itrack++)
    {
        if (Track_list.at(itrack).index_of_closest_topoclusters.size() > 0)
        {
            int closest = Track_list.at(itrack).index_of_closest_topoclusters[0];
            if (closest > -1)
            {
                Supercluster sc;
                sc.set_track(Track_list.at(itrack));
                sc.set_seed(Topo_List.at(closest));
                Super_list.push_back(sc);
            }
        }
    }
}

void Superclustering::sort_by_energy()
{
    std::sort(Super_list.begin(),Super_list.end(),compare_super_e);
}

void Superclustering::add_neighbor_clusters()
{
    //std::vector<std::pair<float,Topo_clust>>;
    std::vector<int> taken_topos;

    for(int isuper=0; isuper<Super_list.size(); isuper++)
    {
        Topo_clust seed = Super_list.at(isuper).get_seed();

        for(int itopo=0; itopo<Topo_List.size(); itopo++)
        {
            Topo_clust test = Topo_List.at(itopo);

            if (std::count(taken_topos.begin(), taken_topos.end(), test.label))
                continue;

            if (seed.label == test.label)
                continue;

            if (in_eta_phi_window(test,seed)){
                Super_list.at(isuper).add_cluster(test);
                taken_topos.push_back(test.label);
            }
        }
    }
}