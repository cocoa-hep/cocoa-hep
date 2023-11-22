#include "Superclustering.hh"
#include "Supercluster.hh"
#include "Topo_clust_var.hh"
#include "Cell_var.hh"


bool in_eta_phi_window(const Topo_clust &topo_1, const Topo_clust &topo_2, std::string wp="tight")
{

    float deta = fabs(topo_1.eta_com - topo_2.eta_com);
    float dphi = fabs(topo_1.phi_com - topo_2.phi_com);
    if (dphi > M_PI)
        dphi -= 2 * M_PI;

    bool passeta = (wp=="tight") ? (deta < 0.075) : (deta < 0.125);
    bool passphi = (wp=="tight") ? (dphi < 0.125) : (dphi < 0.300);
 
    return (passeta && passphi);
}

bool in_eta_phi_window(const Track_struct &track, const Topo_clust &topo)
{

    float deta = fabs(track.eta[1] - topo.eta_com);
    float dphi =      track.phi[1] - topo.phi_com; //signed
    if (fabs(dphi) > M_PI)
        dphi -= 2 * M_PI;

    float q = track.charge;
    bool passeta = (deta   < 0.05);
    bool passphi = (dphi*q < 0.05) && (dphi*q > -0.10);
 
    return (passeta && passphi);
}

float R_distance_func(float phi_1, float eta_1, float phi_2, float eta_2)
{
    float deta = eta_1 - eta_2;
    float dphi = fabs(phi_1 - phi_2);
    if (dphi > M_PI)
        dphi -= 2 * M_PI;
    return sqrtf(sqr(dphi) + sqr(deta));
}

Superclustering::Superclustering(const std::vector<Track_struct> &_track_list, const std::vector<Topo_clust> &_topo_list, const std::vector<Cell *> &_cell_list, std::vector<Supercluster> &Super_list)
{
    cells_in_topoclust = _cell_list;
    Topo_List = _topo_list;
    Track_list = _track_list;

    topo_match_to_tracks();
    find_conversion_vertices();
    find_seed_clusters(Super_list);
    sort_by_energy(Super_list);
    add_neighbor_clusters(Super_list);
    calc_merged_clusters(Super_list);
}

void Superclustering::topo_match_to_tracks()
{
    int size_topo_list = Topo_List.size();
    for(int itopo=0; itopo<size_topo_list; itopo++)
    {
        int size_track_list = Track_list.size();

        for (int itrack = 0; itrack < size_track_list; itrack++)
        {
            if(!in_eta_phi_window(Track_list.at(itrack),Topo_List.at(itopo)))
                continue;

            float R_distance = R_distance_func(Topo_List.at(itopo).phi_com, Topo_List.at(itopo).eta_com, Track_list.at(itrack).phi.at(1), Track_list.at(itrack).eta.at(1));

            std::pair<float,int> dist_itrack_pair;
            dist_itrack_pair.first = R_distance;
            dist_itrack_pair.second = itrack;

            if (Topo_List.at(itopo).closest_tracks.empty())
            {
                Topo_List.at(itopo).closest_tracks.push_back(dist_itrack_pair);
            }
            else if (R_distance < Topo_List.at(itopo).closest_tracks.front().first)
            {
                Topo_List.at(itopo).closest_tracks.front() = dist_itrack_pair;
            }
        }

        std::sort(Topo_List.at(itopo).closest_tracks.begin(),Topo_List.at(itopo).closest_tracks.end()); //Sort based on first element in pair (distance)
    }
}

void Superclustering::find_conversion_vertices()
{
    std::vector<int> used_tracks;
    int size_track_list = Track_list.size();
    for (int itrack = 0; itrack < size_track_list; itrack++)
    {
        if(!Track_list.at(itrack).is_conversion_track)
            continue;

        int n_added = 0;
        ConversionVertex vtx = ConversionVertex();

        //Find the pair of tracks that came from this conversion
        if(vtx.try_add_track(Track_list.at(itrack)))
        {
            n_added++;
            for (int jtrack = 0; jtrack < size_track_list; jtrack++)
            {
                if(vtx.try_add_track(Track_list.at(jtrack)))
                {
                    n_added++;
                    break;
                }
            }
        }

        if(n_added==2)
        {
            //Add topoclusters based on match to conversion tracks
            int added_topos = 0;
            int size_topo_list = Topo_List.size();

            for(int itopo=0; itopo<size_topo_list; itopo++)
            {
                if(vtx.try_add_topo(Topo_List.at(itopo)))
                    added_topos++;
            }

            ConvVertex_list.push_back(vtx);
        }
    }
}

void Superclustering::find_seed_clusters(std::vector<Supercluster> &Super_list)
{

    std::sort(Topo_List.begin(),Topo_List.end(),compare_topo_e);

    std::vector<int> used_tracks;
    int size_topo_list = Topo_List.size();
    for(int itopo=0; itopo<size_topo_list; itopo++)
    {
        Topo_clust topo = Topo_List.at(itopo);
        if(topo.EM_energy < 1000) //Minimum 1 GeV required to qualify for seed candidate
            continue;

        int size_closest_tracks = topo.closest_tracks.size();
        for (int itrack = 0; itrack < size_closest_tracks; itrack++)
        {
            int track_pos_in_list = topo.closest_tracks.at(itrack).second;

            // Check that this track is either from a photon conversion or from a primary electron (i.e. we are not doing electron ID)
            if(!Track_list.at(track_pos_in_list).is_conversion_track && !(abs(Track_list.at(track_pos_in_list).pdgcode)==11))
                continue;

            // Check if this track was already matched to a seed cluster (with higher energy)
            if((std::find(used_tracks.begin(),used_tracks.end(),track_pos_in_list) != used_tracks.end()) && (used_tracks.size() > 0))
                continue;

            // Treatment for topoclusters matched to conversion tracks...
            int matched_conv_vertex = -1;
            int track_friend_idx    = -1;
            if(Track_list.at(track_pos_in_list).is_conversion_track)
            {
                //Check if the other track ("friend") of this conversion vertex was already used in a (higher pT) supercluster
                bool friend_is_matched_to_seed = false;
                int size_ConvVertex_list = ConvVertex_list.size();
                for(int ivtx = 0; ivtx<size_ConvVertex_list; ivtx++)
                {
                    track_friend_idx = ConvVertex_list.at(ivtx).get_friend(Track_list.at(track_pos_in_list));
                    if(track_friend_idx > -1)
                    {
                        matched_conv_vertex = ivtx;
                        //Look if the friend is among the tracks that have already been used for superclusters
                        if((std::find(used_tracks.begin(),used_tracks.end(),track_friend_idx) != used_tracks.end()) && (used_tracks.size() > 0))
                            friend_is_matched_to_seed = true;
                        else
                            friend_is_matched_to_seed = false;
                        break;
                    }
                }

                //In this case, do not create a new supercluster
                if(friend_is_matched_to_seed)
                    continue;
            }

            Supercluster sc;
            sc.set_seed(topo);
            sc.set_track(Track_list.at(track_pos_in_list));
            if(matched_conv_vertex > -1)
            {
                sc.set_conv_track(Track_list.at(track_friend_idx));
                sc.set_conv_vertex(ConvVertex_list.at(matched_conv_vertex));
            }
            sc.distance_seed_track = topo.closest_tracks.at(itrack).first;

            Super_list.push_back(sc);
            used_tracks.push_back(track_pos_in_list);
        }
    }
}

void Superclustering::sort_by_energy(std::vector<Supercluster> &Super_list)
{
    std::sort(Super_list.begin(),Super_list.end(),compare_super_e);
}

void Superclustering::add_neighbor_clusters(std::vector<Supercluster> &Super_list)
{
    std::vector<int> taken_topos;

    //Loop over superclusters (still only seeds)
    for(size_t isuper=0; isuper<Super_list.size(); isuper++)
    {
        Topo_clust seed = Super_list.at(isuper).get_clusters()[0];

        //Look for candidate neighbor clusters by checking criteria
        for(size_t itopo=0; itopo<Topo_List.size(); itopo++)
        {
            Topo_clust test = Topo_List.at(itopo);

            //Check EM energy            
            if (test.EM_energy <= 400)
                continue;

            //Check EM fraction
            if ((test.EM_energy/test.total_energy) <= 0.5)
                continue;

            //Check already-taken list
            if (std::count(taken_topos.begin(), taken_topos.end(), test.label) > 0)
                continue;

            //Check if test is the seed topocluster
            if (seed.label == test.label)
                continue;

            //Check angular window
            if (in_eta_phi_window(test,seed,"tight")){
                Super_list.at(isuper).add_cluster(test);
                taken_topos.push_back(test.label);
            }
            else if(!Super_list.at(isuper).get_track().is_conversion_track){
                //For primary electrons (not photon conversions), 
                if(test.closest_tracks.size()==0)
                    continue;
                //if test shares closest track with seed
                if(test.closest_tracks.at(0).second != Super_list.at(isuper).get_track().position_in_list)
                    continue;
                //check in wider window if test shares closest track with seed
                if(in_eta_phi_window(test,seed,"loose")){
                    Super_list.at(isuper).add_cluster(test);
                    taken_topos.push_back(test.label);
                }
            }
            else //i.e. for a conversion track...
            {
                //Check if associated with same conversion vertex as seed
                if(Super_list.at(isuper).get_conv_vertex().contains(test))
                {
                    Super_list.at(isuper).add_cluster(test);
                    taken_topos.push_back(test.label);
                }
            }
        }
    }
}

void Superclustering::calc_merged_clusters(std::vector<Supercluster> &Super_list)
{

    //Loop over superclusters (still only seeds)
    for(size_t isuper=0; isuper<Super_list.size(); isuper++)
    {
        float sum_cell_e_times_x = 0;
        float sum_cell_e_times_y = 0;
        float sum_cell_e_times_z = 0;
        float sum_cell_e         = 0;

        std::vector<Topo_clust> topo_members = Super_list.at(isuper).get_clusters();

        for(size_t itopo=0; itopo<topo_members.size(); itopo++)
        {
            int topo_label = topo_members.at(itopo).label;        
            int n_cells    = cells_in_topoclust.size();

            for (int cell_i = 0; cell_i < n_cells; ++cell_i)
            {
                if (cells_in_topoclust.at(cell_i)->get_label() != topo_label)
                    continue;
                
                Cell *cell = cells_in_topoclust.at(cell_i);
                if (cell->get_layer() < 3)
                {
                    float en = cell->get_total_energy();
                    sum_cell_e_times_x += en*cell->get_x();
                    sum_cell_e_times_y += en*cell->get_y();
                    sum_cell_e_times_z += en*cell->get_z();
                    sum_cell_e         += en;
                }
            }
        }

        sum_cell_e_times_x /= sum_cell_e;
        sum_cell_e_times_y /= sum_cell_e;
        sum_cell_e_times_z /= sum_cell_e;

        std::vector<float> com = {sum_cell_e_times_x,sum_cell_e_times_y,sum_cell_e_times_z};
        Super_list.at(isuper).com = com;

        float R                   = sqrtf(com[0]*com[0] + com[1]*com[1] + com[2]*com[2]);
	    float theta               =  acos(com[2]/R);
	    Super_list.at(isuper).phi =  atan2(com[1],com[0]);
	    Super_list.at(isuper).eta = -1*log(tan(theta/2.));
    }
}
