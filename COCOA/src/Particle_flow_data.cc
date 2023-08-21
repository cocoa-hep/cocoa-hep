#include "Particle_flow_data.hh"

Particle_flow_data::Particle_flow_data()
{;}

Particle_flow_data::~Particle_flow_data()
{
}

void Particle_flow_data::make_pseudo_jet_particles()
{
    int size_pflow_list = pflow_list.size();
    for (int ipflow = 0; ipflow < size_pflow_list; ipflow++)
    {
        fastjet::PseudoJet particle_flow(pflow_list.at(ipflow).px,
                                                  pflow_list.at(ipflow).py,
                                                  pflow_list.at(ipflow).pz,
                                                  pflow_list.at(ipflow).energy);
        particle_flow.set_user_index(ipflow);
        jets_objects.push_back(particle_flow);
    }
}

void Particle_flow_data::clear()
{
    pflow_list.clear();
    pflow_eta.clear();
    pflow_phi.clear();
    pflow_px.clear();
    pflow_py.clear();
    pflow_pz.clear();
    pflow_e.clear();
    pflow_charge.clear();
    pflow_parent_idx.clear();
    jets_objects.clear();
}

void Particle_flow_data::fill_cell_var()
{
    int size_pflow_list = pflow_list.size();
    for (int ipflow = 0; ipflow < size_pflow_list; ipflow++)
    {
        pflow_eta.push_back(pflow_list.at(ipflow).eta);
        float phi = pflow_list.at(ipflow).phi;
        if (phi < -M_PI) phi += 2*M_PI;
        if (phi > M_PI) phi -= 2*M_PI;
        pflow_phi.push_back(phi);
        pflow_px.push_back(pflow_list.at(ipflow).px);
        pflow_py.push_back(pflow_list.at(ipflow).py);
        pflow_pz.push_back(pflow_list.at(ipflow).pz);
        pflow_e.push_back(pflow_list.at(ipflow).energy);
        pflow_charge.push_back(pflow_list.at(ipflow).charge);
        pflow_parent_idx.push_back(pflow_list.at(ipflow).truth_link);
    }
}


void Particle_flow_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch("pflow_eta"               , "vector<float>", &pflow_eta);
    outTree->Branch("pflow_phi"               , "vector<float>", &pflow_phi);
    outTree->Branch("pflow_px"                , "vector<float>", &pflow_px);
    outTree->Branch("pflow_py"                , "vector<float>", &pflow_py);
    outTree->Branch("pflow_pz"                , "vector<float>", &pflow_pz);
    outTree->Branch("pflow_e"                 , "vector<float>", &pflow_e);
    outTree->Branch("pflow_charge"            , "vector<float>", &pflow_charge);
    outTree->Branch("pflow_parent_idx"        , "vector<int>",   &pflow_parent_idx);
}




