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

bool in_eta_phi_window(const Track_struct &track, const Topo_clust &topo)
{

    float deta = fabs(track.eta[1] - topo.eta_com);
    float dphi =     track.phi[1] - topo.phi_com; //signed
    if (fabs(dphi) > M_PI)
        dphi -= 2 * M_PI;

    float q = track.charge;
    bool passeta = (deta   < 0.05);
    bool passphi = (dphi*q < 0.05) && (dphi*q > -0.10);
 
    return (passeta && passphi);
}

Superclustering::Superclustering(std::vector<Track_struct> &_track_list, std::vector<Topo_clust> &_Topo_List, std::vector<Cell *> &_cell_list, std::vector<Supercluster> &Super_list)
{
    cells_in_topoclust = _cell_list;
    Topo_List = _Topo_List;
    Track_list = _track_list;

    find_seed_clusters(Super_list);
    sort_by_energy(Super_list);
    add_neighbor_clusters(Super_list);
}

void Superclustering::find_seed_clusters(std::vector<Supercluster> &Super_list)
{

    // //Need to start with energy-sorted list of topoclusters, but keep track of original indices
    // std::vector<std::map<int,Topo_clust>> Topo_map_sorted;
    // for(int itopo=0; itopo<size_topo_list; itopo++)
    // {
    //     Topo_map_sorted.push_back(std::pair<itopo,Topo_List.at(itopo)>);
    // }
    std::sort(Topo_List.begin(),Topo_List.end(),compare_topo_e);

    std::vector<int> used_tracks;
    int size_topo_list = Topo_List.size();
    for(int itopo=0; itopo<size_topo_list; itopo++)
    {
        if(Topo_List.at(itopo).total_energy < 1000) //TODO: change to EM energy only. Minimum 1 GeV required to qualify for seed candidate
            continue;
        int ilabel = Topo_List.at(itopo).label;
        std::vector<std::pair<float,int>> closest_tracks;

        int size_track_list = Track_list.size();
        for (int itrack = 0; itrack < size_track_list; itrack++)
        {
            if(std::find(used_tracks.begin(),used_tracks.end(),itrack) != used_tracks.end() && used_tracks.size() > 0)
                continue;

            std::vector<int>    closest_topos  = Track_list.at(itrack).index_of_closest_topoclusters;
            std::vector<double> distance_topos = Track_list.at(itrack).Rprime_to_closest_topoclusters;

            auto itr  = std::find(closest_topos.begin(),closest_topos.end(),ilabel);
            if(itr != closest_topos.end())
            {
                int place      = std::distance(closest_topos.begin(), itr);
                float distance = distance_topos.at(place);
                if(in_eta_phi_window(Track_list.at(itrack),Topo_List.at(itopo)))
                {
                    closest_tracks.push_back(std::pair(distance_topos.at(place),itrack));
                }
            }
        }

        if(closest_tracks.size()>0)
        {
            std::sort(closest_tracks.begin(),closest_tracks.end()); //Sort based on first element in pair (distance)
            Supercluster sc;
            sc.set_seed(Topo_List.at(itopo));
            sc.set_track(Track_list.at(closest_tracks.at(0).second));
            sc.distance_seed_track = closest_tracks.at(0).first;
            Super_list.push_back(sc);
            used_tracks.push_back(closest_tracks.at(0).second);
        }
    }
}

void Superclustering::sort_by_energy(std::vector<Supercluster> &Super_list)
{
    std::sort(Super_list.begin(),Super_list.end(),compare_super_e);
}

void Superclustering::add_neighbor_clusters(std::vector<Supercluster> &Super_list)
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