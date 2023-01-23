#ifndef __DEBUG_PARTICLE_FLOW_DATA_H__
#define __DEBUG_PARTICLE_FLOW_DATA_H__


#include <vector>
#include "TTree.h"
#include "Pion_for_pflow_var.hh"

class Debug_Particle_Flow_data
{

public:
    Debug_Particle_Flow_data(/* args */);
    ~Debug_Particle_Flow_data();
    std::vector<Pion_for_pflow_var> charge_pions;
    std::vector<int>   pion_in_track_list;
    std::vector<float> epsilon_lead_list;
    std::vector<float> pion_eta_list;
    std::vector<float> pion_pt_list;
    std::vector<int>   n_90_list;
    std::vector<float> rho_list;
    std::vector<float> sigma_phi_list;
    std::vector<float> sigma_eta_list;
    float Eref;
    float LHED;
    std::vector <float> E_p_prime_list;
    std::vector <float> E_p_second_list;
    std::vector <float> deltaRPrime_prime_list;
    std::vector <float> deltaRPrime_second_list;
    std::vector <int> is_leading_clust_eq_to_closest;
    std::vector <float> deltaR_list;
    std::vector <float> S_prime_list;
    std::vector <float> epsilon_matched_list;
    std::vector <float> rho_matched_list;
    static Debug_Particle_Flow_data &GetInstance()
    {
        static Debug_Particle_Flow_data pion_info;
        return pion_info;
    };
    void set_tree_branches(TTree *outTree, std::string Type_Of_Runnig);
    void clear();
    void fill_var();
    std::string type_of_runnig;
};


#endif // __DEBUG_PARTICLE_FLOW_DATA_H__