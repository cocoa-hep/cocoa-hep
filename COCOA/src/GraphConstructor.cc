#include "GraphConstructor.hh"

GraphConstructor::GraphConstructor(std::vector<Cell*> &low_cells_in_topoclusters, std::vector<Cell*> &high_cells_in_topoclusters, std::vector<Track_struct> &tracks_list, std::vector<G4int> &particle_to_track, Graph_construction_data &graph_obj,
				     std::vector<float> *_particle_dep_energies)
{
	fill_cell_to_cell_edges(low_cells_in_topoclusters, graph_obj.cell_to_cell_edge_start, graph_obj.cell_to_cell_edge_end);
	fill_track_to_cell_edges(low_cells_in_topoclusters, tracks_list, graph_obj.track_to_cell_edge_start, graph_obj.track_to_cell_edge_end);
	fill_particle_to_node_edges(low_cells_in_topoclusters, particle_to_track, graph_obj.particle_to_node_idx, graph_obj.particle_to_node_weight, _particle_dep_energies);
	fill_superres_cell_to_cell_edges(high_cells_in_topoclusters, low_cells_in_topoclusters,
									 graph_obj.high_cell_to_low_cell_edge, graph_obj.low_cell_to_high_cell_edge);//! Probably some mistake
}

GraphConstructor::GraphConstructor(std::vector<Cell*> &low_cells_in_topoclusters, std::vector<Track_struct> &tracks_list, std::vector<G4int> &particle_to_track, Graph_construction_data &graph_obj,
				     std::vector<float> *_particle_dep_energies)
{
	fill_cell_to_cell_edges(low_cells_in_topoclusters, graph_obj.cell_to_cell_edge_start, graph_obj.cell_to_cell_edge_end);
	fill_track_to_cell_edges(low_cells_in_topoclusters, tracks_list, graph_obj.track_to_cell_edge_start, graph_obj.track_to_cell_edge_end);
	fill_particle_to_node_edges(low_cells_in_topoclusters, particle_to_track, graph_obj.particle_to_node_idx, graph_obj.particle_to_node_weight, _particle_dep_energies);
}

void GraphConstructor::fill_cell_to_cell_edges(std::vector<Cell*> &cell,
											   std::vector<int> &cell_to_cell_edge_start, std::vector<int> &cell_to_cell_edge_end)
{

	std::map<std::pair<int, int>, int> edge_list;
	int n_cells = cell.size();
	const int N = n_cells;
	TMatrixD dr_array(N, N);

	// 1) Fill NxN matrix with dR of cell pairs
	for (int cell_i = 0; cell_i < n_cells; ++cell_i)
	{
		TVector3 v_cell_i(cell.at(cell_i)->get_x(), cell.at(cell_i)->get_y(), cell.at(cell_i)->get_z());
		int layer_i = cell.at(cell_i)->get_layer();
		for (int cell_j = cell_i + 1; cell_j < n_cells; ++cell_j)
		{
			TVector3 v_cell_j(cell.at(cell_j)->get_x(), cell.at(cell_j)->get_y(), cell.at(cell_j)->get_z());
			int layer_j = cell.at(cell_j)->get_layer();
			float dr = 999;

			if( std::abs(layer_i - layer_j) <= config_json_var.graph_construction.max_layer_sep ){
				dr = v_cell_i.DeltaR(v_cell_j);
			}

			dr_array[cell_i][cell_j] = dr;
			dr_array[cell_j][cell_i] = dr;
		}
	}

    // 2) Double for-loop over cell i, cell j and fill priority queue based on dR of pairs
	for (int cell_i = 0; cell_i < n_cells; ++cell_i)
	{

		//int next_layer_cell = -1;
		//float next_layer_cell_dr = 100.0;

		//int layer_below_cell = -1;
		//float layer_below_cell_dr = 100.0;

		int layer_i = cell.at(cell_i)->get_layer();

        // List of neighbor cell_j candidates and their dr w.r.t. cell_i, ordered by dr (LARGEST dr served first)
        // Note that in case of std::pair, priority_queue orders based on first element in pair (i.e. dr in our case)
		std::priority_queue<std::pair<float, int>> pq;        //cells in same layer
        std::priority_queue<std::pair<float,int>> pq_inter;   //cells in separate layers

		for (int cell_j = 0; cell_j < n_cells; ++cell_j)
		{
			if (cell_j == cell_i)
			{
				continue;
			}

			int layer_j = cell.at(cell_j)->get_layer();

			// If too many layers apart...
			if( std::abs(layer_i - layer_j) > config_json_var.graph_construction.max_layer_sep ){
				continue;
			}

        	// Look up the dr of this pair from NxN matrix
			float dr = dr_array[cell_i][cell_j];

         	// If too far separated in dR...
			if(dr > config_json_var.graph_construction.max_dr[layer_i]){
				continue;
			}

			// same-layer edges:
			if(layer_j == layer_i)
			{
				std::pair<float,int> cell_dr_pair = std::make_pair(dr,cell_j);
				// If the queue is already at or above capacity, and if this candidate pair is closer than the largest pair in the queue...
				if(pq.size() == config_json_var.graph_construction.max_samelayer_edges[layer_i] && (dr < pq.top().first) )
				{
					pq.push(cell_dr_pair); // insert this pair into queue
					pq.pop();              // and remove highest-priority (i.e. largest dr) pair
				}
				else if(pq.size() < config_json_var.graph_construction.max_samelayer_edges[layer_i] ) //just add the pair without checking dr
				{
					pq.push(cell_dr_pair);
				}
				else continue;
			}
			else // inter-layer edges:
			{
				if(pq_inter.size() == config_json_var.graph_construction.max_interlayer_edges[layer_i] && (dr < pq_inter.top().first) )
				{
					std::pair<float,int> cell_dr_pair_inter = std::make_pair(dr,cell_j);
					pq_inter.push(cell_dr_pair_inter); // insert this pair into queue
					pq_inter.pop();                    // and remove highest-priority (i.e. largest dr) pair
				}
				else if(pq_inter.size() < config_json_var.graph_construction.max_interlayer_edges[layer_i] )
				{ //just add the pair without checking dr
					std::pair<float,int> cell_dr_pair_inter = std::make_pair(dr,cell_j);
					pq_inter.push(cell_dr_pair_inter);
				}
				else continue;
			}

		} // end loop over cell j

		//int added_edges = 0;
		//int found_edges = 0;
		size_t n_edges_to_add = pq.size();
        if(n_edges_to_add > config_json_var.graph_construction.max_samelayer_edges[layer_i]){
        	std::cout << "WARNING: n_edges_to_add = " << n_edges_to_add << " for layer " << layer_i << " (expected " << config_json_var.graph_construction.max_samelayer_edges[layer_i] << ")" << std::endl;
        }

        // Add (incoming) edges to cell_i based on the ordering of the priority queue
		for (size_t edge_i = 0; edge_i < n_edges_to_add; edge_i++)
		{
			std::pair<float,int> cell_j_pair = pq.top(); //start with furthest away
			pq.pop(); // won't be needing you anymore

			int j = cell_j_pair.second;
			int i = cell_i;

			std::pair<int,int> cand_edge = std::make_pair(i,j);

			// If this pair already exists, don't add again
			if(edge_list.find(cand_edge)!=edge_list.end())
			{
				//std::cout << "found {i,j} = " << i << ", " << j << " (skipping!)" << std::endl;
				continue;
			}
			//added_edges = 1;
			edge_list.insert(std::pair<std::pair<int,int>,int> (cand_edge,1));
        }

        size_t n_edges_to_add_inter = pq_inter.size();

        // Do same for inter-layer edges, update same edge_list
        for(size_t edge_i=0; edge_i < n_edges_to_add_inter; edge_i++)
		{
			std::pair<float,int> cell_j_pair = pq_inter.top();
			pq_inter.pop();

			std::pair<int,int> cand_edge = std::make_pair(cell_i,cell_j_pair.second);

			if(edge_list.find(cand_edge)!=edge_list.end())
			{
				continue;
			}
			edge_list.insert(std::pair<std::pair<int,int>,int> (cand_edge,1));
        }

		//if (added_edges < 1)
		//{
		//	std::cout << " cell " << cell_i << " has no edges " << found_edges << std::endl;
		//}
	} // end loop over cell i

	for (auto &[a, b] : edge_list)
	{
		//start = source node index (j), end = destination node index (i)
		cell_to_cell_edge_start.push_back( a.second);
		cell_to_cell_edge_end.push_back( a.first );
	}
}

void GraphConstructor::fill_track_to_cell_edges(std::vector<Cell*> &cell,
												std::vector<Track_struct> &track, std::vector<int> &track_to_cell_edge_start, std::vector<int> &track_to_cell_edge_end)
{
	int n_tracks = track.size();
	int n_cells = cell.size();

	for (int track_i = 0; track_i < n_tracks; ++track_i)
	{
		if (!track.at(track_i).Is_Track_Useable()){
			continue;
		}
		for (int lay=0; lay<6; ++lay)
		{
			TVector3 v_track;
			//float theta_thislayer = 2*atan(exp(-1.*track.at(track_i).eta[lay]));
			//v_track.SetMagThetaPhi(1, theta_thislayer, track.at(track_i).phi[lay]);
			v_track.SetPtEtaPhi(track.at(track_i).pt,track.at(track_i).eta[lay],track.at(track_i).phi[lay]);

			std::priority_queue<std::pair<float, int>> pq;
			for (int cell_i = 0; cell_i < n_cells; ++cell_i)
			{
				TVector3 v_cell(cell.at(cell_i)->get_x(), cell.at(cell_i)->get_y(), cell.at(cell_i)->get_z());

				int cell_l = cell.at(cell_i)->get_layer();

                if(cell_l != lay){
                  continue;
                }
			
				float dr = v_cell.DeltaR(v_track); //dr w.r.t. projection of this track in the layer cell_l
				if(dr > config_json_var.graph_construction.max_celltrack_dr[cell_l])
				{
					continue;
				}

                if(pq.size() >= config_json_var.graph_construction.max_celltrack_edges[cell_l] && (dr < pq.top().first) )
				{
                    std::pair<float,int> cell_dr_pair = std::make_pair(dr,cell_i);
                    pq.push(cell_dr_pair);
                    pq.pop();
                }
				else if(pq.size() < config_json_var.graph_construction.max_celltrack_edges[cell_l])
				{
                    std::pair<float,int> cell_dr_pair = std::make_pair(dr,cell_i);
                    pq.push(cell_dr_pair);
                }
			}

			size_t n_edges_to_add = pq.size();
			for (size_t edge_i = 0; edge_i < n_edges_to_add; edge_i++)
			{
				std::pair<float, int> edge_pair = pq.top();

				track_to_cell_edge_start.push_back(track_i);
				track_to_cell_edge_end.push_back(edge_pair.second);

				pq.pop();
			}
		}
	}
}

//input are the cell to particle associations, and the fraction of energy in each cell for that given particle. Then we need to add the tracks

void GraphConstructor::fill_particle_to_node_edges( std::vector<Cell*> &cell,
						    std::vector<int> particle_to_track,
						    std::vector<std::vector<float>> &particle_to_node_idx,
						    std::vector<std::vector<float>> &particle_to_node_weight,
						    std::vector<float> *_particle_dep_energies )
{
	//remember to initialize particle_to_node_weight and particle_to_node_idx in the main particle loop
	int n_particles = particle_to_track.size();
	int n_cells = cell.size();
	particle_to_node_idx.resize(n_particles);
	particle_to_node_weight.resize(n_particles);
	std::vector<float> cell_energy_particle_j(n_particles, 0);
	for (int cell_i = 0; cell_i < n_cells; cell_i++)
	{
		// int n_parent_cell_i = fcell_parent_list.at(cell_i).size();
		std::vector<Particle_dep_in_cell> cell_particle;
		cell.at(cell_i)->get_particles(cell_particle);
		int n_parent_cell_i = cell_particle.size();
		for (int parent_j_cell_i = 0; parent_j_cell_i < n_parent_cell_i; parent_j_cell_i++)
		{
			int particle_idx = cell_particle.at(parent_j_cell_i).particle_pos_in_true_list;
			particle_to_node_idx.at(particle_idx).push_back(cell_i);
			float energy = cell_particle.at(parent_j_cell_i).Energy;
			particle_to_node_weight.at(particle_idx).push_back(energy);
			cell_energy_particle_j.at(particle_idx) = cell_energy_particle_j.at(particle_idx) + energy;
		}
	}
	if ( _particle_dep_energies )
	    *_particle_dep_energies = cell_energy_particle_j;
	//normalize cells energy now
	for (int particle_i = 0; particle_i < n_particles; particle_i++)
	{
		float energy = cell_energy_particle_j.at(particle_i);
		int ncells = particle_to_node_weight.at(particle_i).size();
		float attenuation = 1.;
		if (particle_to_track.at(particle_i) != -1)
			attenuation = 0.5;
		energy = energy / attenuation;
		for (int cell_j = 0; cell_j < ncells; cell_j++)
		{
			particle_to_node_weight.at(particle_i).at(cell_j) = particle_to_node_weight.at(particle_i).at(cell_j) / energy;
		}
	}
	//now add a loop for the tracks - shift the index by n_cells
	for (int particle_i = 0; particle_i < n_particles; particle_i++)
	{
		float track_idx = particle_to_track.at(particle_i);
		if (track_idx != -1)
		{
			track_idx = track_idx + n_cells;
			track_idx = (float)track_idx;
			particle_to_node_idx.at(particle_i).push_back(track_idx);
			particle_to_node_weight.at(particle_i).push_back(0.5);
		}
	}
}

void GraphConstructor::fill_superres_cell_to_cell_edges(std::vector<Cell*> &high_cells, std::vector<Cell*> &low_cells,
														std::vector<int> &high_cell_to_low_cell_edge, std::vector<int> &low_cell_to_high_cell_edge)
{
	int n_high_cell = high_cells.size();
	int n_low_cell = low_cells.size();
	for (int ilow_cell = 0; ilow_cell < n_low_cell; ilow_cell++)
	{
		Cell *low_cell = low_cells.at(ilow_cell);
		for (int ihigh_cell = 0; ihigh_cell < n_high_cell; ihigh_cell++)
		{
			Cell *high_cell = high_cells.at(ihigh_cell);
			std::vector<std::vector<int>> cell_link;
			high_cell->get_cell_link_superres(cell_link);
			if ((low_cell->get_layer() == cell_link.front().at(0)) &&
				(low_cell->get_eta() == cell_link.front().at(1)) &&
				(low_cell->get_phi() == cell_link.front().at(2)))
			{
				low_cell_to_high_cell_edge.push_back(ilow_cell);
				high_cell_to_low_cell_edge.push_back(ihigh_cell);
			}
		}
	}
}
