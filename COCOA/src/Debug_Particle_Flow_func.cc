#include "Debug_Particle_Flow_func.hh"

bool sortbysec(const std::pair<float, float> &a, const std::pair<float, float> &b)
{
    return (a.first > b.first);
}

template <typename T>
int Debug_Particle_Flow_func::sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

Debug_Particle_Flow_func::Debug_Particle_Flow_func(std::vector<Track_struct> track_list,
                                                   std::vector<Topo_clust> &topo_list,
                                                   std::vector<Cell *> &topo_cells,
                                                   Debug_Particle_Flow_data &debug_data,
                                                   Geometry_definition Geometry,
                                                   Particle_flow_alg_var Pflow_Param,
                                                   std::string type_of_running) : cells_in_topoclust(topo_cells), topoclusters_list(topo_list)
{
    pflow_param = Pflow_Param;
    geometry = Geometry;
    Containment_shower(track_list, debug_data);
    PFlow_parameters(debug_data, type_of_running);
}
void Debug_Particle_Flow_func::Containment_shower(std::vector<Track_struct> &track_list,
                                                  Debug_Particle_Flow_data &debug_data)
{
    int size_track_list = track_list.size();
    for (int itrack = 0; itrack < size_track_list; itrack++)
    {
        Track_struct pion_trk = track_list[itrack];
        if (fabs(pion_trk.pdgcode) == 211)
        {
            if (pion_trk.Is_Track_Useable())
            {
                std::vector<std::pair<float, float>> Cluster_conter;
                std::vector<int> clusters_labels;
                std::vector<float> clustEnergy;
                std::vector<int> ind_of_clust_contain_dep;
                std::vector<float> energy_deposit_in_clust;
                int size_cells_in_topoclust = cells_in_topoclust.size();
                for (int icell = 0; icell < size_cells_in_topoclust; icell++)
                {
                    Cell *cell = cells_in_topoclust.at(icell);

                    std::vector<Particle_dep_in_cell> list_of_particles_dep_energy;
                    cell->get_particles(list_of_particles_dep_energy);
                    int size_list_of_particles_dep_energy = list_of_particles_dep_energy.size();
                    for (int iparticle = 0; iparticle < size_list_of_particles_dep_energy; iparticle++)
                    {
                        Particle_dep_in_cell particle_dep_energy = list_of_particles_dep_energy.at(iparticle);
                        if (particle_dep_energy.particle_pos_in_true_list == pion_trk.nFinal_State_Particles)
                        {
                            int size_clusters_labels = clusters_labels.size();
                            bool IsNewClust = true;
                            for (int clust = 0; clust < size_clusters_labels; clust++)
                            {
                                
                                if (clusters_labels[clust] == cell->get_label())
                                {
                                    clustEnergy[clust] += list_of_particles_dep_energy[iparticle].Energy;
                                    IsNewClust = false;
                                }
                            }
                            if (IsNewClust)
                            {
                                clusters_labels.push_back(cell->get_label());
                                clustEnergy.push_back(list_of_particles_dep_energy[iparticle].Energy);
                            }
                        }
                    }
                }
                int size_clusters_labels = clusters_labels.size();
                float total_dep_energy_by_pion = 0;
                float totalenergyClusters = 0;
                for (int clust = 0; clust < size_clusters_labels; clust++)
                {
                    Cluster_conter.push_back(std::make_pair(clustEnergy[clust], clusters_labels[clust]));
                    total_dep_energy_by_pion += clustEnergy[clust];
                    totalenergyClusters += topoclusters_list.at(clusters_labels[clust] - 1).total_energy;
                }
                sort(Cluster_conter.begin(), Cluster_conter.end(), sortbysec);

                for (int itopo = 0; itopo < size_clusters_labels; itopo++)
                {
                    ind_of_clust_contain_dep.push_back((int)Cluster_conter.at(itopo).second);
                    energy_deposit_in_clust.push_back(Cluster_conter.at(itopo).first);
                }
                if (total_dep_energy_by_pion > 0.2 * pion_trk.energy)
                {
                    int n_90 = 0;
                    float energy = 0;
                    for (int clust = 0; clust < size_clusters_labels; clust++)
                    {
                        if (energy < 0.9 * total_dep_energy_by_pion)
                        {
                            energy += Cluster_conter[clust].first;
                            n_90 += 1;
                        }
                    }
                    Pion_for_pflow_var pion(track_list.at(itrack));
                    pion.energies_deposited_in_clusters = energy_deposit_in_clust;
                    pion.indexes_of_clusters_contained_dep = ind_of_clust_contain_dep;
                    pion.total_dep_energies_by_pions = total_dep_energy_by_pion;
                    float rho = Cluster_conter.front().first / (topoclusters_list.at(Cluster_conter[0].second - 1).total_energy - topoclusters_list.at(Cluster_conter.front().second - 1).noise);
                    float epsilon = Cluster_conter.front().first / total_dep_energy_by_pion;
                    pion.n_90 = n_90;
                    pion.epsilon_lead = epsilon;
                    pion.sigma_eta = topoclusters_list.at(Cluster_conter.front().second - 1).sigma_eta;
                    pion.sigma_phi = topoclusters_list.at(Cluster_conter.front().second - 1).sigma_phi;
                    pion.rho = rho;
                    debug_data.charge_pions.push_back(pion);
                }
            }
        }
    }
}

void Debug_Particle_Flow_func::PFlow_parameters(Debug_Particle_Flow_data &debug_data,
                                                std::string type_of_running)
{
    E_div_p_template = 0;
    sigma_E_div_p_template = 0;
    delta_Rprime_threshold = 0;
    layer_max = -1;
    int size_charge_pions = debug_data.charge_pions.size();
    Topo_List = topoclusters_list;
    Type_Of_Runnig = type_of_running;
    if (config_var.Type_of_running.find("PFlow_debug") != std::string::npos)
    {
        if (Type_Of_Runnig == "PFlow_debug_E_p_template")
        {
            int size_topo_list = topoclusters_list.size();
            for (int ipion = 0; ipion < size_charge_pions; ipion++)
            {
                Pion_for_pflow_var &pion = debug_data.charge_pions.at(ipion);
                if (abs(pion.pdgcode) != 11 && pion.Is_Track_Useable())
                {
                    for (int itopo = 0; itopo < size_topo_list; itopo++)
                    {
                        float R_dist = Particle_flow_func::R_distance_func(pion.phi.at(1), pion.eta.at(1),
                                                                           topoclusters_list.at(itopo).phi_com, topoclusters_list.at(itopo).eta_com);
                        if (R_dist < 0.4)
                        {
                            E_div_p_template += topoclusters_list.at(itopo).total_energy;
                        }
                    }
                    E_div_p_template /= pion.p_perigee;
                    debug_track_match_to_topoclusters(debug_data.charge_pions.at(ipion));
                    layer_max = Particle_flow_func::LHED(pion, geometry, cells_in_topoclust);
                }
                pion.Eref = E_div_p_template;
                pion.LHED = layer_max;
            }
        }
        else
        {
            for (int ipion = 0; ipion < size_charge_pions; ipion++)
            {
                if (abs(debug_data.charge_pions.at(ipion).pdgcode) != 11 && debug_data.charge_pions.at(ipion).Is_Track_Useable())
                {
                    debug_track_match_to_topoclusters(debug_data.charge_pions.at(ipion));
                    debug_parameters(debug_data.charge_pions.at(ipion));
                }
            }
        }
    }
}
void Debug_Particle_Flow_func::debug_track_match_to_topoclusters(Pion_for_pflow_var &pion)
{
    int size_topo_list = Topo_List.size();
    for (int itopo = 0; itopo < size_topo_list; itopo++)
    {
        float R_dist = Particle_flow_func::R_distance_func(pion.phi.at(1), pion.eta.at(1), Topo_List.at(itopo).phi_com, Topo_List.at(itopo).eta_com,
                                                           Topo_List.at(itopo).sigma_phi, Topo_List.at(itopo).sigma_eta);
        if (pion.Rprime_to_closest_topoclusters.empty())
        {
            pion.Rprime_to_closest_topoclusters.push_back(R_dist);
            pion.index_of_closest_topoclusters.push_back(itopo + 1);
        }
        else
        {
            bool check = true;
            int size_closest_topoclusters = pion.index_of_closest_topoclusters.size();
            for (int iclosest = 0; iclosest < size_closest_topoclusters; iclosest++)
            {
                if (R_dist < pion.Rprime_to_closest_topoclusters.at(iclosest))
                {
                    pion.Rprime_to_closest_topoclusters.insert(pion.Rprime_to_closest_topoclusters.begin() + iclosest, R_dist);
                    pion.index_of_closest_topoclusters.insert(pion.index_of_closest_topoclusters.begin() + iclosest, itopo + 1);
                    check = false;
                    break;
                }
            }
            if (check)
            {
                pion.Rprime_to_closest_topoclusters.push_back(R_dist);
                pion.index_of_closest_topoclusters.push_back(itopo + 1);
            }
        }
    }
}
void Debug_Particle_Flow_func::debug_parameters(Pion_for_pflow_var &pion)
{
    int size_topo_list = Topo_List.size();
    
    float topo_clust_energy = Topo_List.at(pion.index_of_closest_topoclusters.front() - 1).total_energy;
    float R_dist = Particle_flow_func::R_distance_func(Topo_List.at(pion.index_of_closest_topoclusters.front() - 1).phi_com,
                                                       Topo_List.at(pion.index_of_closest_topoclusters.front() - 1).eta_com,
                                                       pion.phi.at(1), pion.eta.at(1));
    pion.is_leading_clust_eq_to_closest = pion.index_of_closest_topoclusters.front() == pion.indexes_of_clusters_contained_dep.front();

    pion.E_p_prime = topo_clust_energy / pion.p_perigee;
    pion.deltaRPrime_prime = pion.Rprime_to_closest_topoclusters.front();
    if (size_topo_list > 1)
    {
        pion.E_p_second = Topo_List.at(pion.index_of_closest_topoclusters.at(1) - 1).total_energy / pion.p_perigee;
        pion.deltaRPrime_second = pion.Rprime_to_closest_topoclusters.at(1);
    }
    else
    {
        pion.E_p_second = -100;
        pion.deltaRPrime_second = -100;
    }
    pion.deltaR = R_dist;
    layer_max = Particle_flow_func::LHED(pion, geometry, cells_in_topoclust);
    std::tie(E_div_p_template, sigma_E_div_p_template) = Particle_flow_func::energy_div_momentum_template(pion.pt_perigee, layer_max, pion.eta.at(0));
    float S_discriminant = (topo_clust_energy - pion.p_perigee * E_div_p_template) / (sigma_E_div_p_template * pion.p_perigee);

    int num_assosiated_clust_pion = pion.indexes_of_clusters_contained_dep.size();
    float epsilon_matched = 0;
    float rho_matched = 0;
    for (int ind = 0; ind < num_assosiated_clust_pion; ind++)
    {
        if (pion.index_of_closest_topoclusters[0] == pion.indexes_of_clusters_contained_dep.at(ind))
        {
            epsilon_matched = pion.energies_deposited_in_clusters.at(ind) / pion.total_dep_energies_by_pions;                                          //ChPionsInfo[pion][ind+1]/ChPionsInfo[pion][5];
            rho_matched = pion.energies_deposited_in_clusters.at(ind) / Topo_List.at(pion.indexes_of_clusters_contained_dep.at(ind) - 1).total_energy; //ChPionsInfo[pion][ind+1]/(TList[(int)ChPionsInfo[pion][ind]-1][0]);
            break;
        }
    }
    pion.S_prime = S_discriminant;
    pion.epsilon_matched = epsilon_matched;
    pion.rho_matched = rho_matched;
}