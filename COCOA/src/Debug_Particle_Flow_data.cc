#include "Debug_Particle_Flow_data.hh"


Debug_Particle_Flow_data::Debug_Particle_Flow_data(/* args */)
{
}

Debug_Particle_Flow_data::~Debug_Particle_Flow_data()
{
}
void Debug_Particle_Flow_data::set_tree_branches(TTree *outTree, std::string Type_Of_Runnig)
{
    type_of_runnig = Type_Of_Runnig;
    outTree->Branch("pion_in_tracks_list"   , "vector<int>",      &pion_in_track_list);
	outTree->Branch("pion_eta", "vector<float>",                    &pion_eta_list);
	outTree->Branch("pion_pt", "vector<float>",                      &pion_pt_list);
	
    if (type_of_runnig == "PFlow_debug_E_p_template")
    {
        outTree->Branch("E_div_p_template", &Eref, "Eref/F");
		outTree->Branch("LHED", &LHED, "LHED/F");
        outTree->Branch("pion_sigma_phi", "vector<float>",                    &sigma_phi_list);
	    outTree->Branch("pion_sigma_eta", "vector<float>",                    &sigma_eta_list);
    }
    else if (type_of_runnig == "PFlow_debug")
    {
        outTree->Branch("pion_epsilon_lead", "vector<float>",          &epsilon_lead_list);
        outTree->Branch("pion_n_90", "vector<int>",                      &n_90_list);
	    outTree->Branch("pion_rho", "vector<float>",                      &rho_list);
        outTree->Branch("E_p_prime_list", "vector<float>", &E_p_prime_list);
		outTree->Branch("E_p_second_list", "vector<float>", &E_p_second_list);
		outTree->Branch("deltaRPrime_prime_list", "vector<float>", &deltaRPrime_prime_list);
		outTree->Branch("deltaRPrime_second_list", "vector<float>", &deltaRPrime_second_list);
		outTree->Branch("deltaR_list", "vector<float>", &deltaR_list);
		outTree->Branch("S_prime_list", "vector<float>", &S_prime_list);
		outTree->Branch("is_leading_clust_eq_to_closest", "vector<int>", &is_leading_clust_eq_to_closest);
		outTree->Branch("epsilon_matched_list","vector<float>",&epsilon_matched_list);
		outTree->Branch("rho_matched_list","vector<float>",&rho_matched_list);
    }

}
void Debug_Particle_Flow_data::clear()
{
    charge_pions.clear();
    pion_in_track_list.clear();
    epsilon_lead_list.clear();
    pion_eta_list.clear();
    pion_pt_list.clear();
    n_90_list.clear();
    rho_list.clear();
    sigma_phi_list.clear();
    sigma_eta_list.clear();
    Eref = 0;
    LHED = 0;
    E_p_prime_list.clear();
    E_p_second_list.clear();
    deltaRPrime_prime_list.clear();
    deltaRPrime_second_list.clear();
    is_leading_clust_eq_to_closest.clear();
    deltaR_list.clear();
    S_prime_list.clear();
    epsilon_matched_list.clear();
    rho_matched_list.clear();
}

void Debug_Particle_Flow_data::fill_var()
{
    for (Pion_for_pflow_var &ipion: charge_pions)
    {
        pion_in_track_list.push_back(ipion.position_in_list);
        epsilon_lead_list.push_back(ipion.epsilon_lead);
        pion_eta_list.push_back(ipion.eta.at(0));
        pion_pt_list.push_back(ipion.pt_perigee);
        n_90_list.push_back(ipion.n_90);
        rho_list.push_back(ipion.rho);
        sigma_phi_list.push_back(ipion.sigma_phi);
        sigma_eta_list.push_back(ipion.sigma_eta);
        if (type_of_runnig == "PFlow_debug_E_p_template")
        {
            Eref = ipion.Eref;
            LHED = ipion.LHED;
        }
        else if (type_of_runnig == "PFlow_debug")
        {
            E_p_prime_list.push_back(ipion.E_p_prime);
            E_p_second_list.push_back(ipion.E_p_second);
            deltaRPrime_prime_list.push_back(ipion.deltaRPrime_prime);
            deltaRPrime_second_list.push_back(ipion.deltaRPrime_second);
            is_leading_clust_eq_to_closest.push_back(ipion.is_leading_clust_eq_to_closest);
            deltaR_list.push_back(ipion.deltaR);
            S_prime_list.push_back(ipion.S_prime);
            epsilon_matched_list.push_back(ipion.epsilon_matched);
            rho_matched_list.push_back(ipion.rho_matched);
        }
    }
}

