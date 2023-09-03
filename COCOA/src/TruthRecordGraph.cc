#include "TruthRecordGraph.hh"
#include "OutputRunAction.hh"

void TruthRecordGraph::clear()
{
	m_interesting_particles.clear();
	m_final_state_particles.clear();
	node_idx.clear();
	node_pdg_id.clear();
	node_phi.clear();
	node_eta.clear();
	node_pt.clear();
	node_m.clear();
	node_isfinal.clear();
	node_prodx.clear();
	node_prody.clear();
	node_prodz.clear();
	final_idx.clear();
	node_decx.clear();
	node_decy.clear();
	node_decz.clear();
	edge_start.clear();
	edge_end.clear();
}

bool TruthRecordGraph::isCharmHadron(int pdgid)
{
	std::vector<int>::iterator it = std::find(CharmHadrons.begin(), CharmHadrons.end(), fabs(pdgid));
	if (it != CharmHadrons.end())
	{
		return true;
	}
	return false;
}

bool TruthRecordGraph::isBottomHadron(int pdgid)
{
	std::vector<int>::iterator it = std::find(BottomHadrons.begin(), BottomHadrons.end(), fabs(pdgid));
	if (it != BottomHadrons.end())
	{
		return true;
	}
	return false;
}

void TruthRecordGraph::add_to_vector(HepMC::GenParticle *to_add, std::vector<HepMC::GenParticle *> &add_to)
{
	std::vector<HepMC::GenParticle *>::iterator it = std::find(add_to.begin(), add_to.end(), to_add);
	if (it == add_to.end())
	{
		add_to.push_back(to_add);
	}
}

bool TruthRecordGraph::is_parent_same(HepMC::GenParticle *particle)
{

	if (!(particle->production_vertex()))
	{
		return false;
	}

	int pdgid = particle->pdg_id();

	for (HepMC::GenVertex::particle_iterator
			 parent = particle->production_vertex()->particles_begin(HepMC::parents);
		 parent != particle->production_vertex()->particles_end(HepMC::parents); ++parent)
	{

		int parent_pdgid = (*parent)->pdg_id();
		if (parent_pdgid == pdgid)
		{
			return true;
		}
	}
	return false;
}

HepMC::GenParticle *TruthRecordGraph::find_next_level_parent(HepMC::GenParticle *particle)
{

	for (HepMC::GenVertex::particle_iterator
			 parent = particle->production_vertex()->particles_begin(HepMC::parents);
		 parent != particle->production_vertex()->particles_end(HepMC::parents); ++parent)
	{

		float r = prod_radius(*parent);

		if (r < m_max_radius)
		{
			return (*parent);
		}
		else
		{
			return find_next_level_parent(*parent);
		}
	}
	return nullptr;
}

float TruthRecordGraph::prod_radius(HepMC::GenParticle *particle)
{
	HepMC::FourVector pos_start = (particle)->production_vertex()->position();
	//HepMC::FourVector pos_end = (particle)->end_vertex()->position();

	float r = TMath::Sqrt(pos_start.x() * pos_start.x() + pos_start.y() * pos_start.y());
	return r;
}

HepMC::GenParticle *TruthRecordGraph::check_prod_location(HepMC::GenParticle *particle)
{
	float r = prod_radius(particle);

	if (r > m_max_radius)
	{
		return find_next_level_parent(particle);
	}
	return particle;
}

void TruthRecordGraph::add_all_moving_parents(HepMC::GenParticle *to_add, std::vector<HepMC::GenParticle *> &add_to)
{

	if (!(to_add->production_vertex()))
	{
		return;
	}

	for (HepMC::GenVertex::particle_iterator parent = to_add->production_vertex()->particles_begin(HepMC::parents);
		 parent != to_add->production_vertex()->particles_end(HepMC::parents);
		 ++parent)
	{

		int parent_pdgid = (*parent)->pdg_id();
		bool is_charm = isCharmHadron(parent_pdgid);
		bool is_b = isBottomHadron(parent_pdgid);

		bool is_top = abs(parent_pdgid) == 6;
		bool is_up = parent_pdgid==2;
		bool is_W = abs(parent_pdgid) == 24;
		// bool is_gluon = abs(parent_pdgid)==21;
		/* Has child correct PDG code? */
		if (is_charm || is_b)
		{
			add_to_vector(*parent, add_to);
		}

		if (is_top || is_W || is_up)
		{
			if (!is_parent_same(*parent))
			{
				add_to_vector(*parent, add_to);
			}
		}

		if (!((*parent)->production_vertex()))
		{
			add_to_vector(*parent, add_to);
		}

		if (((*parent)->end_vertex()) && ((*parent)->production_vertex()))
		{
			HepMC::FourVector pos_start = (*parent)->production_vertex()->position();
			HepMC::FourVector pos_end = (*parent)->end_vertex()->position();
			TVector3 pos_start_vec(pos_start.x(), pos_start.y(), pos_start.z());
			TVector3 pos_end_vec(pos_end.x(), pos_end.y(), pos_end.z());

			float flight_distance = (pos_start_vec - pos_end_vec).Mag();
			if (flight_distance > 0.000001)
			{
				add_to_vector(*parent, add_to);
			}
		}

		add_all_moving_parents(*parent, add_to);
	}
}

bool TruthRecordGraph::is_parent_of(HepMC::GenParticle *parent, HepMC::GenParticle *potential_child)
{
	if (!(parent->end_vertex()))
	{
		return false;
	}

	bool found_match_for_potential_child = false;

	for (HepMC::GenVertex::particle_iterator
			 child = parent->end_vertex()->particles_begin(HepMC::children);
		 child != parent->end_vertex()->particles_end(HepMC::children); ++child)
	{

		if (*child == potential_child)
		{
			return true;
		}
		else
		{
			if (is_parent_of(*child, potential_child))
			{
				found_match_for_potential_child = true;
			}
		}
	}
	
	if (found_match_for_potential_child)
	{
		return true;
	}

	return false;
}

void TruthRecordGraph::clean_daugthers(std::vector<HepMC::GenParticle *> &interesting_particles,
									   std::vector<int> &all_daughters, std::vector<int> &direct_daughters)
{

	int n_daughters = all_daughters.size();

	for (int i = 0; i < n_daughters; ++i)
	{
		bool has_parent_in_daughters = false;
		for (int j = 0; j < n_daughters; ++j)
		{
			if (i == j)
			{
				continue;
			}

			if (is_parent_of(interesting_particles.at(all_daughters.at(j)), interesting_particles.at(all_daughters.at(i))))
			{
				has_parent_in_daughters = true;
			}
		}
		if (!has_parent_in_daughters)
		{
			direct_daughters.push_back(all_daughters.at(i));
		}
	}
}

void TruthRecordGraph::find_daughters(HepMC::GenParticle *parent, std::vector<HepMC::GenParticle *> &interesting_particles, std::vector<int> &direct_daughters)
{

	if (!(parent->end_vertex()))
	{
		return;
	}

	for (HepMC::GenVertex::particle_iterator
			 child = parent->end_vertex()->particles_begin(HepMC::children);
		 child != parent->end_vertex()->particles_end(HepMC::children); ++child)
	{
		std::vector<HepMC::GenParticle *>::iterator it = std::find(interesting_particles.begin(), interesting_particles.end(), *child);
		if (it != interesting_particles.end())
		{
			int idx = std::distance(interesting_particles.begin(), it);
			direct_daughters.push_back(idx);
		}
		else
		{
			find_daughters(*child, interesting_particles, direct_daughters);
		}
	}
}

void TruthRecordGraph::fill_truth_graph()
{
	size_t n_particles = m_interesting_particles.size();
	for (size_t part_i = 0; part_i < n_particles; part_i++)
	{
		HepMC::GenParticle *particle = m_interesting_particles.at(part_i);
		std::vector<int> all_direct_daughters;
		std::vector<int> direct_daughters;

		find_daughters(particle, m_interesting_particles, all_direct_daughters);
		clean_daugthers(m_interesting_particles, all_direct_daughters, direct_daughters);

		for (auto d_i : direct_daughters)
		{
			edge_start.push_back(part_i);
			edge_end.push_back(d_i);
		}
	}

	for (size_t part_i = 0; part_i < n_particles; part_i++)
	{

		bool has_mothers_or_daughters = false;

		std::vector<int>::iterator it = std::find(edge_start.begin(), edge_start.end(), part_i);
		if (it != edge_start.end())
		{
			has_mothers_or_daughters = true;
		}
		it = std::find(edge_end.begin(), edge_end.end(), part_i);
		if (it != edge_end.end())
		{
			has_mothers_or_daughters = true;
		}

		if (has_mothers_or_daughters)
		{

			node_idx.push_back(part_i);
			HepMC::GenParticle *particle = m_interesting_particles.at(part_i);

			int node_final_state_idx = -1;

			std::vector<HepMC::GenParticle *>::iterator find_it = std::find(m_final_state_particles.begin(), m_final_state_particles.end(), particle);
			if (find_it != m_final_state_particles.end())
			{
				node_final_state_idx = abs(std::distance(m_final_state_particles.end(), find_it)) - 1;
			}

			node_pdg_id.push_back(particle->pdg_id());
			node_phi.push_back(particle->momentum().phi());
			node_eta.push_back(particle->momentum().eta());
			node_pt.push_back(particle->momentum().perp());
			node_m.push_back(particle->momentum().m());
			node_isfinal.push_back(particle->status());
			final_idx.push_back(node_final_state_idx);
		}
	}
}

void TruthRecordGraph::set_tree_branches(TTree *outTree)
{
	outTree->Branch("node_idx", &node_idx);
	outTree->Branch("node_pdg_id", &node_pdg_id);
	outTree->Branch("node_phi", &node_phi);
	outTree->Branch("node_eta", &node_eta);
	outTree->Branch("node_pt", &node_pt);
	outTree->Branch("node_m", &node_m);
	outTree->Branch("node_isfinal", &node_isfinal);
	outTree->Branch("node_prodx", &node_prodx);
	outTree->Branch("node_prody", &node_prody);
	outTree->Branch("node_prodz", &node_prodz);
	outTree->Branch("node_final_state_idx", &final_idx); 
	outTree->Branch("node_decx", &node_decx);
	outTree->Branch("node_decy", &node_decy);
	outTree->Branch("node_decz", &node_decz);
	outTree->Branch("node_edge_start", &edge_start);
	outTree->Branch("node_edge_end", &edge_end);
}
