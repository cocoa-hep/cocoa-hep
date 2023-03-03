#ifndef __PION_FOR_PFLOW_VAR_H__
#define __PION_FOR_PFLOW_VAR_H__

#include "Track_var.hh"

class Pion_for_pflow_var
{
public:
    Pion_for_pflow_var(Track_struct track);
    int pdgcode;                //code of particle
    int nFinal_State_Particles; //where particle in the list of FinalStateParticles
    double energy;
    double px;
    double py;
    double pz;
    double absmom;
    std::vector<double> eta;         //eta coordinate of Particle in each layer
    std::vector<double> phi;         //phi coordinate of Particle in each layer
    std::vector<double> Rprime_to_closest_topoclusters;
    std::vector<int> index_of_closest_topoclusters;

    int position_in_list;

    bool Is_track_reconstracted;
    bool Is_inside_R;
    bool Is_reach_calorimeter;
    bool Is_Track_Useable();


    double charge;
    double alpha;
    double px_end_MF;
    double py_end_MF;
    double pz_perigee;
    double p_perigee;
    double px_perigee;
    double py_perigee;
    double pt_perigee;
    double energy_perigee;
    double mass_perigee;
    double x_end_MF;
    double y_end_MF;
    double z_end_MF;


    std::vector<int> indexes_of_clusters_contained_dep;
    std::vector<float> energies_deposited_in_clusters;
    float total_dep_energies_by_pions;
    int n_90;
    float epsilon_lead;
    float pion_eta;
    float pion_pt;
    float rho;
    float sigma_phi;
    float sigma_eta;
    float E_p_prime;
    float E_p_second;
    float deltaRPrime_prime;
    float deltaRPrime_second;
    float is_leading_clust_eq_to_closest;
    float deltaR;
    float S_prime;
    float epsilon_matched;
    float rho_matched;
    float Eref;
    float LHED;
};


#endif // __PION_FOR_PFLOW_VAR_H__