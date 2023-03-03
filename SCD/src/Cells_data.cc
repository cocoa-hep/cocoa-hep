#include "Cells_data.hh"


Cells_data::Cells_data(bool is_high)
{
	high = is_high;
	if (high)
		Number_Pixel_Flatten = config_var.high_resolution.number_of_pixels_flatten;
	else
		Number_Pixel_Flatten = config_var.low_resolution.number_of_pixels_flatten;
	ilay_num = Number_Pixel_Flatten.size();
	
	for (int ilow = 0; ilow < ilay_num; ilow++)
	{
		std::vector<std::vector<Cell>> vect_cell(Number_Pixel_Flatten.at(ilow),
												 std::vector<Cell>(Number_Pixel_Flatten.at(ilow)));
		fCell_array.push_back(vect_cell);
		for (int ieta = 0; ieta < Number_Pixel_Flatten.at(ilow); ieta++)
		{
			for (int iphi = 0; iphi < Number_Pixel_Flatten.at(ilow); iphi++)
			{
				fCell_array.at(ilow).at(ieta).at(iphi).Reset();
				fCell_array.at(ilow).at(ieta).at(iphi).set_indexes(ilow, ieta, iphi, high);
			}
		}
	}
}

void Cells_data::clear()
{
	Cells_in_topoclusters.clear();
	fCell_array.clear();
	cell_pflow_object_idx.clear();
	cell_layer.clear();
	cell_x.clear();
	cell_y.clear();
	cell_z.clear();
	cell_eta.clear();
	cell_phi.clear();
	cell_e.clear();
	cell_che.clear();
	cell_nue.clear();
	cell_topo_idx.clear();
	cell_parent_idx.clear();
	cell_conv_el_idx.clear();
	cell_parent_list.clear();
	cell_parent_energy.clear();
	for (int ilow = 0; ilow < ilay_num; ilow++)
	{
		std::vector<std::vector<Cell>> vect_cell(Number_Pixel_Flatten.at(ilow),
												 std::vector<Cell>(Number_Pixel_Flatten.at(ilow)));
		fCell_array.push_back(vect_cell);
		for (int ieta = 0; ieta < Number_Pixel_Flatten.at(ilow); ieta++)
		{
			for (int iphi = 0; iphi < Number_Pixel_Flatten.at(ilow); iphi++)
			{
				fCell_array.at(ilow).at(ieta).at(iphi).Reset();
				fCell_array.at(ilow).at(ieta).at(iphi).set_indexes(ilow, ieta, iphi, high);
			}
		}
	}
}

void Cells_data::ChangeLabelForCells(int OldLabel, int NewLabel)
{
	for (int ilow = 0; ilow < ilay_num; ilow++)
	{
		for (int ieta = 0; ieta < Number_Pixel_Flatten.at(ilow); ieta++)
		{
			for (int iphi = 0; iphi < Number_Pixel_Flatten.at(ilow); iphi++)
			{
				Cell &local_cell = fCell_array.at(ilow).at(ieta).at(iphi);
				int cell_lab = local_cell.get_label();
				if (cell_lab == OldLabel)
					local_cell.set_label(NewLabel);
				else if (cell_lab > OldLabel)
					local_cell.set_label(cell_lab - 1);
				local_cell.get_eta_pos();
			}
		}
	}
}

void Cells_data::set_tree_branches(TTree *outTree)
{
	//* cell branches
    if ( config_var.doPFlow )
	    outTree->Branch("cell_pflow_object_idx"   , "vector<int>",   &cell_pflow_object_idx);
	outTree->Branch("cell_layer", "vector<int>", &cell_layer);
	outTree->Branch("cell_x", "vector<float>", &cell_x);
	outTree->Branch("cell_y", "vector<float>", &cell_y);
	outTree->Branch("cell_z", "vector<float>", &cell_z);
	outTree->Branch("cell_eta", "vector<float>", &cell_eta);
	outTree->Branch("cell_phi", "vector<float>", &cell_phi);
	outTree->Branch("cell_e", "vector<float>", &cell_e);
	outTree->Branch("cell_che", "vector<float>", &cell_che);
	outTree->Branch("cell_nue", "vector<float>", &cell_nue);
	outTree->Branch("cell_topo_idx", "vector<int>", &cell_topo_idx);
	outTree->Branch("cell_parent_idx", "vector<int>", &cell_parent_idx);
	outTree->Branch("cell_conv_el_idx", "vector<int>", &cell_conv_el_idx);
	outTree->Branch("cell_parent_list", "vector< vector <float> >", &cell_parent_list);
	outTree->Branch("cell_parent_energy", "vector< vector <float> >", &cell_parent_energy);
}

void Cells_data::fill_cell_var()
{    
	Full_trajectory_info_data &fsp_obj = Full_trajectory_info_data::GetInstance();
	Particle_flow_data &pflow_obj = Particle_flow_data::GetInstance();
	int size_fsp_obj = fsp_obj.fAllTrajectoryInfo.size();
	int size_pflow = pflow_obj.pflow_list.size();
	int cell_counter = 0;
	int size_Cells_in_topoclusters = Cells_in_topoclusters.size();
	for (int icell = 0; icell < size_Cells_in_topoclusters; icell++)
	{
		Cell *local_cell = Cells_in_topoclusters.at(icell);
		cell_layer.push_back(local_cell->get_layer());
		cell_x.push_back(local_cell->get_x());
		cell_y.push_back(local_cell->get_y());
		cell_z.push_back(local_cell->get_z());
		cell_eta.push_back(local_cell->get_eta_pos());
		cell_phi.push_back(local_cell->get_phi_pos());
		cell_e.push_back(local_cell->get_total_energy());
		cell_che.push_back(local_cell->get_charge_energy());
		cell_nue.push_back(local_cell->get_neutral_energy());
		cell_topo_idx.push_back(local_cell->get_label());
		for (int ipflow = 0; ipflow < size_pflow; ipflow++)
		{
			if ((pflow_obj.pflow_list.at(ipflow).charge == 0) && (local_cell->get_label() ==pflow_obj.pflow_list.at(ipflow).label))
			{
				cell_pflow_object_idx.push_back(ipflow);
			}
		}

		std::vector<float> valueE;
		std::vector<float> valueP;
		std::vector<Particle_dep_in_cell> ParticleVec;
		ParticleVec.reserve(5000);
		local_cell->get_particles(ParticleVec);
		//* (Positiion in list of stabel particles, energy deposition)

		for (std::vector<Particle_dep_in_cell>::iterator it1 = ParticleVec.begin(); it1 != ParticleVec.end(); ++it1)
		{
			valueE.push_back(it1->Energy);
			valueP.push_back(it1->particle_pos_in_true_list);
		}
		int n_parents = valueE.size();
		int parent_idx = -1;

		if (n_parents == 0)
		{
			cell_parent_idx.push_back(-1);
		}
		else
		{
			std::vector<float>::iterator max_parent;
			max_parent = std::max_element(valueE.begin(), valueE.end());
			int max_idx = std::distance(valueE.begin(), max_parent);
			parent_idx = valueP.at(max_idx);
			cell_parent_idx.push_back(parent_idx);
		}

		std::vector<Particle_dep_in_cell> conv_electrons;
		conv_electrons.reserve(1000);
		local_cell->get_conv_electrons(conv_electrons);
		
		if ( !conv_electrons.size() )
		    cell_conv_el_idx.push_back(-1);
		else {
		    float energyDepositedMax  = 0.0;
		    int   conv_el_maxEdep_idx = 0;
		    for( size_t iConvEl = 0; iConvEl < conv_electrons.size(); ++iConvEl ) {
			if ( conv_electrons[iConvEl].Energy > energyDepositedMax ) {
			    energyDepositedMax  = conv_electrons[iConvEl].Energy;
			    conv_el_maxEdep_idx = conv_electrons[iConvEl].particle_pos_in_true_list;
			}
		    }
		    cell_conv_el_idx.push_back( conv_el_maxEdep_idx );
		}
		
		std::vector<float> cell_parents;
		std::vector<float> cell_i_parent_energy;

		for (int parent_i = 0; parent_i < n_parents; parent_i++)
		{
			float parent_energy = valueE.at(parent_i);

			float parent_idx_float = valueP.at(parent_i);
			if (parent_energy / local_cell->get_total_energy() > 0.01)
			{
				cell_parents.push_back(parent_idx_float);
				cell_i_parent_energy.push_back(parent_energy);
			}
		}
		cell_parent_list.push_back(cell_parents);
		cell_parent_energy.push_back(cell_i_parent_energy);
		cell_counter++;
	}
}

void Cells_data::add_cell_info(int ilay, int ieta, int iphi, float ch_en, float nu_en, Particle_dep_in_cell particle, Particle_dep_in_cell* conv_el)
{
	fCell_array.at(ilay).at(ieta).at(iphi).add_charge_energy(ch_en);
	fCell_array.at(ilay).at(ieta).at(iphi).add_neutral_energy(nu_en);
	fCell_array.at(ilay).at(ieta).at(iphi).add_particle(particle);
	if ( conv_el )
	    fCell_array.at(ilay).at(ieta).at(iphi).add_particle( *conv_el, true );
	
}

void Cells_data::fill_cells_in_topoclusters()
{
	float sum = 0;
	int pos_in_list = 0;
    for (int ilay = 0; ilay < ilay_num; ilay++)
    {
        for (int ieta = 0; ieta < (int) Number_Pixel_Flatten.at(ilay); ieta++)
        {
            for (int iphi = 0; iphi < (int) Number_Pixel_Flatten.at(ilay); iphi++)
            {
                Cell &local_cell = fCell_array.at(ilay).at(ieta).at(iphi);
                if (local_cell.get_label() != 0)
                {
                    local_cell.position_in_list = pos_in_list;
                    Cells_in_topoclusters.push_back(&fCell_array.at(ilay).at(ieta).at(iphi));
                    pos_in_list++;
					sum+=local_cell.get_total_energy();
                }
            }
        }
    }
}
