#include "Particle_flow_func.hh"
#include "DetectorConstruction.hh"

float Particle_flow_func::phi_distance(float phi_1, float phi_2)
{
    float dphi = fabs(phi_1 - phi_2);
    if (dphi > M_PI)
        dphi -= 2 * M_PI;
    return dphi;
}

float Particle_flow_func::R_distance_func(float phi_1, float eta_1, float phi_2, float eta_2, float sigma_phi, float sigma_eta)
{
    float d_phi = phi_distance(phi_1, phi_2);
    float d_eta = eta_1 - eta_2;
    return sqrtf(sqr(d_phi / sigma_phi) + sqr(d_eta / sigma_eta));
}
//* Addapd only for Low resolution
Particle_flow_func::Particle_flow_func(std::vector<Track_struct> &track_list, std::vector<Topo_clust> &topo_list,
                                       std::vector<Cell *> &topo_cells, std::vector<Particle_flow_var> &pflow_list,
                                       Geometry_definition Geometry, Particle_flow_alg_var Pflow_Param) : cells_in_topoclust(topo_cells)
{
    pflow_param = Pflow_Param;
    geometry = Geometry;
    E_div_p_template = 0;
    sigma_E_div_p_template = 0;
    Type_Of_Runnig = "Standard";
    layer_max = -1;
    int size_track_list = track_list.size();
    Topo_List = topo_list;

    for (int itrack = 0; itrack < size_track_list; itrack++)
    {
        if ((abs(track_list.at(itrack).pdgcode)!=13) && track_list.at(itrack).Is_Track_Useable())
        {
            track_match_to_topoclusters(track_list.at(itrack));
            if (track_list.at(itrack).Rprime_to_closest_topoclusters.at(0) < delta_Rprime_threshold)
            {
                if (abs(track_list.at(itrack).pdgcode)!=11)
                {
                    layer_max = LHED(track_list.at(itrack), geometry, cells_in_topoclust);
                    track_list.at(itrack).SetLHED( layer_max );
                    recovering_shower_split(track_list.at(itrack));
                    cell_by_cell_subst(track_list.at(itrack));
                }
                else
                {
                    for (auto tc_idx : track_list.at(itrack).index_of_closest_topoclusters)
		    {
                        Topo_List.at(tc_idx-1).charge = track_list.at(itrack).charge;
                    }
                }
            }
        }
    }

    for (int itopo = 0; itopo < (int)Topo_List.size(); itopo++)
    {
        if (Topo_List.at(itopo).total_energy > 0)
        {
            pflow_list.push_back(Topo_List.at(itopo));
        }
    }

    for (int itrack = 0; itrack < size_track_list; itrack++)
    {
        if (abs(track_list.at(itrack).pdgcode)!=11 && track_list.at(itrack).Is_Track_Useable() )
            pflow_list.emplace_back(track_list.at(itrack));
    }    
}

void Particle_flow_func::track_match_to_topoclusters(Track_struct &track)
{
    int size_topo_list = Topo_List.size();
    int size_momentum_delta_Rprime_threshold = pflow_param.momentum_delta_Rprime_threshold.size();
    delta_Rprime_threshold = pflow_param.delta_Rprime_threshold.back();
    for (int ipar = 0; ipar < size_momentum_delta_Rprime_threshold; ipar++)
    {
        if (track.p_perigee <= pflow_param.momentum_delta_Rprime_threshold.at(ipar))
            delta_Rprime_threshold = pflow_param.delta_Rprime_threshold.at(ipar);
    }

    for (int itopo_clust = 0; itopo_clust < size_topo_list; itopo_clust++)
    {
        if (Topo_List.at(itopo_clust).total_energy / track.p_perigee > pflow_param.E_div_p_threshold)
        {
            float R_distance = R_distance_func(track.phi.at(1), track.eta.at(1),
                                               Topo_List.at(itopo_clust).phi_com,
                                               Topo_List.at(itopo_clust).eta_com,
                                               Topo_List.at(itopo_clust).sigma_phi,
                                               Topo_List.at(itopo_clust).sigma_eta);
            if (track.Rprime_to_closest_topoclusters.empty())
            {
                track.Rprime_to_closest_topoclusters.push_back(R_distance);
                track.index_of_closest_topoclusters.push_back(Topo_List.at(itopo_clust).label);
            }
            else if (R_distance < track.Rprime_to_closest_topoclusters.front())
            {
                track.Rprime_to_closest_topoclusters.front() = R_distance;
                track.index_of_closest_topoclusters.front() = Topo_List.at(itopo_clust).label;
            }
        }
    }
    if (track.Rprime_to_closest_topoclusters.empty())
    {
        track.Rprime_to_closest_topoclusters.push_back(1000);
        track.index_of_closest_topoclusters.push_back(-1);
    }
}



float Particle_flow_func::cell_volume(Cell &cell, Geometry_definition geometry_)
{

    Config_reader_var &config_var = Config_reader_var::GetInstance();
    
    long double r_inn_0               = geometry_.layer_inn_radius_flatten[0];
    int         iCellLayer            = cell.get_layer();
    long double r_inn_layer           = geometry_.layer_inn_radius_flatten[iCellLayer];
    long double depth_cell            = geometry_.layer_out_radius_flatten[iCellLayer] - r_inn_layer;
    long double depth_previous_layers = r_inn_layer - geometry_.layer_inn_radius_flatten[0];
    	
    long double theta_min_barrel = EtaToTheta( config_var.max_eta_barrel );
    long double cell_eta         = cell.get_eta_pos();
    long double dEta             = geometry_.layer_deta_flatten[iCellLayer];
    long double dPhi             = geometry_.layer_dphi_flatten[iCellLayer];
    long double theta_cell_max   = EtaToTheta( fabs( cell_eta ) - 0.5 * dEta );
    long double theta_cell_min   = EtaToTheta( fabs( cell_eta ) + 0.5 * dEta );
    
    long double volume = 0.0;
    if ( fabs( cell_eta ) < config_var.max_eta_barrel ) {
	//
	// Barrel cells
	//
	long double r_inn_cell  = r_inn_0 + sin( theta_cell_max ) * depth_previous_layers;
	volume                  = 1.0 / 3.0 * dPhi * ( tan( 0.5 * M_PI - theta_cell_min ) - tan( 0.5 * M_PI - theta_cell_max ) );
	volume                 *= powf( r_inn_cell + sin( theta_cell_max ) * depth_cell , 3) - powf(r_inn_cell, 3);
    } else {
	//
	// Endcap cells
	//
	long double z_min  = r_inn_0 / tan( theta_min_barrel ) + cos( theta_cell_min ) * depth_previous_layers;
	volume             = 1.0 / 6.0 * dPhi * ( powf( tan( theta_cell_max ), 2) - powf( tan( theta_cell_min ), 2) );
	volume            *= powf( z_min + cos( theta_cell_max ) * depth_cell, 3 ) - powf( z_min, 3 );
	
    }
    return volume;

}

void Particle_flow_func::recovering_shower_split(Track_struct &track)
{
    int size_topo_list = Topo_List.size();
    std::tie(E_div_p_template, sigma_E_div_p_template) = energy_div_momentum_template(track.pt, layer_max, track.eta.at(0));
    float sigma_e_dep = sigma_E_div_p_template * track.p_perigee;
    float e_dep = track.p_perigee * E_div_p_template;
    float S_discriminant = (Topo_List.at(track.index_of_closest_topoclusters.front() - 1).total_energy - e_dep)/sigma_e_dep;
    if ( S_discriminant < pflow_param.S_discriminant_threshold)
    {
        for (int itopo = 0; itopo < size_topo_list; itopo++)
        {
            float R_distance = R_distance_func(track.phi.at(1), track.eta.at(1),
                                               Topo_List.at(itopo).phi_com, Topo_List.at(itopo).eta_com);
            if ((track.index_of_closest_topoclusters.at(0) - 1 != itopo) && (R_distance < 0.4))//Todo
            {
                track.index_of_closest_topoclusters.push_back(Topo_List.at(itopo).label);
                R_distance = R_distance_func(track.phi.at(1), track.eta.at(1),
                                             Topo_List.at(itopo).phi_com, Topo_List.at(itopo).eta_com,
                                             Topo_List.at(itopo).sigma_phi, Topo_List.at(itopo).sigma_eta);
                track.Rprime_to_closest_topoclusters.push_back(R_distance);
            }
        }
    }
}

std::tuple<float, float> Particle_flow_func::energy_div_momentum_template(float pt, int lhed, float eta)
{
    CSVReader &csv_reader = CSVReader::GetInstance();
    return {csv_reader.hist_pflow_mean_template.GetBinContent(csv_reader.hist_pflow_mean_template.FindBin(pt, fabs(eta), lhed) ),//+ 0.2
            csv_reader.hist_pflow_std_template.GetBinContent(csv_reader.hist_pflow_std_template.FindBin(pt, fabs(eta), lhed))};
}

void Particle_flow_func::cell_by_cell_subst(Track_struct &track)
{
    // if (layer_max)
    int num_matched_clusters = track.index_of_closest_topoclusters.size();
    float E_clust_tot = 0;
    std::vector<Cell *> cell_list;
    float rest_energy = 0;
    int size_cells_in_topoclust = cells_in_topoclust.size();
    float sigma_e_dep = sigma_E_div_p_template * track.p_perigee;
    float e_dep = track.p_perigee * E_div_p_template;
    //
    // Sum up energies of cells associated with TCs matched to the track.
    // Also add these cells to 'cell_list'.
    //
    // At this point Cell_var::get_pflow_remnant returns the total cell energy.
    //
    for (int icell = 0; icell < size_cells_in_topoclust; icell++)
    {
        Cell *cell = cells_in_topoclust.at(icell);

        float R_distance = R_distance_func(track.phi.at(layer_max), track.eta.at(layer_max),
                                           cell->get_phi_pos(), cell->get_eta_pos());
        cell->set_distance_to_track(R_distance);
        if (cell->get_is_cell_shared())
        {
            for (int itopo_match = 0; itopo_match < num_matched_clusters; itopo_match++)
            {
                if (cell->get_1st_local_max_label() == track.index_of_closest_topoclusters.at(itopo_match))
                {
                    cell->set_is_1st_local_max_assoisited_with_track(true);
                    E_clust_tot += cell->get_pflow_remnant() * cell->get_1st_local_max_weight();
                    cell_list.push_back(cell);
                    rest_energy += cell->get_pflow_remnant() * cell->get_1st_local_max_weight();
                }
                else if (cell->get_2nd_local_max_label() == track.index_of_closest_topoclusters.at(itopo_match))
                {
                    cell->set_is_2nd_local_max_assoisited_with_track(true);
                    E_clust_tot += cell->get_pflow_remnant() * cell->get_2nd_local_max_weight();
                    cell_list.push_back(cell);
                    rest_energy += cell->get_pflow_remnant() * cell->get_2nd_local_max_weight();
                }
            }
        }
        else
        {
            for (int itopo_match = 0; itopo_match < num_matched_clusters; itopo_match++)
            {
                if (cell->get_label() == track.index_of_closest_topoclusters.at(itopo_match))
                {
                    E_clust_tot += cell->get_pflow_remnant();
                    Cell *cell_for_list = cells_in_topoclust.at(icell);
                    cell_list.push_back(cell_for_list);
                    rest_energy += cell_for_list->get_pflow_remnant();
                    break;
                }
            }
        }
    }
    
    //
    // sum_e_dep : sum of all mean expected deposited energies of tracks matched to a TC
    //

    if (e_dep >= E_clust_tot)
    {
	//
	// Total deposited energy below the mean exp. deposited energy : subtract all the energy from the cells in question
	//
        sum_e_dep += e_dep;
        for (int icell = 0; icell < size_cells_in_topoclust; icell++)
        {
            for (int itopo_match = 0; itopo_match < num_matched_clusters; itopo_match++)
            {
                if ((cells_in_topoclust.at(icell)->get_label() == track.index_of_closest_topoclusters.at(itopo_match)) && (cells_in_topoclust.at(icell)->get_pflow_remnant()>0))
                {
                    cells_in_topoclust.at(icell)->subtract_energys_fraction_from_pflow_remnant(1); //! Topo_list
                    Topo_List.at(cells_in_topoclust.at(icell)->get_label() - 1).subtract_cell(*cells_in_topoclust.at(icell), 1);
                }
            }
        }
        rest_energy = 0;
    }
    else
    {
        int num_ring = 1;
        std::vector<Ring_substraction> rings;
        float size_cell_list = cell_list.size();
        for (int ilay = 0; ilay < (int)geometry.kNLayers; ilay++)
        {
            float factor = 1.2 * geometry.layer_deta_flatten.at(ilay);
            for (int iradius = 0; iradius < (int)(geometry.number_of_pixels_flatten.at(ilay) / 2) + 1; iradius++)
            {
                Ring_substraction ring;
                ring.layer = ilay;
                ring.ring_number = num_ring;
                int num_cells_in_ring = 0;
                for (int icell = 0; icell < size_cell_list; icell++)
                {
                    if ((cell_list.at(icell)->get_layer() == ilay) &&
                        (sqr((iradius + 1) * factor) > sqr(cell_list.at(icell)->get_distance_to_track())) &&
                        (sqr(iradius * factor) <= sqr(cell_list.at(icell)->get_distance_to_track())))
                    {
                        num_cells_in_ring++;
                        if (!cell_list.at(icell)->get_is_cell_shared())
                        {
                            ring.energy += cell_list.at(icell)->get_pflow_remnant();
                            ring.energy_dencity += cell_list.at(icell)->get_energy_dencity();
                        }
                        if (cell_list.at(icell)->get_is_1st_local_max_assoisited_with_track())
                        {
                            ring.energy += cell_list.at(icell)->get_1st_local_max_weight() * cell_list.at(icell)->get_pflow_remnant();
                            ring.energy_dencity += cell_list.at(icell)->get_1st_local_max_weight() * cell_list.at(icell)->get_energy_dencity();
                        }
                        if (cell_list.at(icell)->get_is_2nd_local_max_assoisited_with_track())
                        {
                            ring.energy += cell_list.at(icell)->get_2nd_local_max_weight() * cell_list.at(icell)->get_pflow_remnant();
                            ring.energy_dencity += cell_list.at(icell)->get_2nd_local_max_weight() * cell_list.at(icell)->get_energy_dencity();
                        }
                        cell_list.at(icell)->set_number_of_ring(num_ring);
                    }
                }
                num_ring++;
                rings.push_back(ring);
            }
        }
        int size_rings = rings.size();
        sort(rings.begin(), rings.end(), sort_by_energy_dencity);
        for (int iring = 0; iring < size_rings; iring++)
        {
            if (e_dep <= rings.at(iring).energy)
            {
                sum_e_dep += e_dep;
                float ch_frac = e_dep / rings.at(iring).energy;
                e_dep = 0;
                for (int icell = 0; icell < size_cell_list; icell++)
                {
                    if ((cell_list.at(icell)->get_number_of_ring() == rings.at(iring).ring_number) && (cell_list.at(icell)->get_pflow_remnant()>0))
                    {
                        rest_energy -= ch_frac * cell_list.at(icell)->get_pflow_remnant();
                        cell_list.at(icell)->subtract_energys_fraction_from_pflow_remnant(ch_frac);
                        Topo_List.at(cell_list.at(icell)->get_label() - 1).subtract_cell(*cell_list.at(icell), ch_frac);
                    }
                }
                break;
            }
            else
            {
                e_dep -= rings.at(iring).energy;
                sum_e_dep += rings.at(iring).energy;
                for (int icell = 0; icell < size_cell_list; icell++)
                {
                    if ((cell_list.at(icell)->get_number_of_ring() == rings.at(iring).ring_number) && (cell_list.at(icell)->get_pflow_remnant()>0))
                    {
                        rest_energy -= cell_list.at(icell)->get_pflow_remnant();
                        cell_list.at(icell)->subtract_energys_fraction_from_pflow_remnant(1);
                        Topo_List.at(cell_list.at(icell)->get_label() - 1).subtract_cell(*cell_list.at(icell), 1);
                    }
                }
            }
        }

        if (rest_energy <= pflow_param.factor_sigma_E_div_p_template * sigma_e_dep)
        {
            sum_e_dep += rest_energy;
            for (int icell = 0; icell < size_cell_list; icell++)
            {
                if (cell_list.at(icell)->get_pflow_remnant()>0)
                {
                    cells_in_topoclust.at(cell_list.at(icell)->position_in_list)->subtract_energys_fraction_from_pflow_remnant(1);
                    Topo_List.at(cell_list.at(icell)->get_label() - 1).subtract_cell(*cell_list.at(icell), 1);
                }
            }
            rest_energy = 0;
        }
    }
}
