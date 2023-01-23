#include "Pion_for_pflow_var.hh"

Pion_for_pflow_var::Pion_for_pflow_var(Track_struct track)
{
    pdgcode = track.pdgcode;
    nFinal_State_Particles = track.nFinal_State_Particles;
    energy = track.energy;
    px = track.px;
    py = track.py;
    pz = track.pz;
    absmom = track.absmom;
    eta = track.eta;
    phi = track.phi;
    Rprime_to_closest_topoclusters = track.Rprime_to_closest_topoclusters;
    index_of_closest_topoclusters = track.index_of_closest_topoclusters;
    position_in_list = track.position_in_list;
    Is_track_reconstracted = track.Is_track_reconstracted;
    Is_inside_R = track.Is_inside_R;
    Is_reach_calorimeter = track.Is_reach_calorimeter;
    charge = track.charge;
    alpha = track.alpha;
    p_perigee = track.p_perigee;
    px_end_MF = track.px_end_MF;
    py_end_MF = track.py_end_MF;
    pz_perigee = track.pz_perigee;
    px_perigee = track.px_perigee;
    py_perigee = track.py_perigee;
    pt_perigee = track.pt_perigee;
    energy_perigee = track.energy_perigee;
    mass_perigee = track.mass_perigee;
    x_end_MF = track.x_end_MF;
    y_end_MF = track.y_end_MF;
    z_end_MF = track.z_end_MF;
    indexes_of_clusters_contained_dep.clear();
    energies_deposited_in_clusters.clear();
    total_dep_energies_by_pions = 0;
}

bool Pion_for_pflow_var::Is_Track_Useable()
{
    return Is_track_reconstracted && Is_inside_R && Is_reach_calorimeter;
}