#include "Particle_flow_var.hh"


Particle_flow_var::Particle_flow_var() 
{
    isTrack = false;
    charge = 0;
    label  = 0;
    eta = 0;
    phi = 0;
    energy = 0;
    px = 0;
    py = 0;
    pz = 0;
    truth_link = -1;
}

Particle_flow_var::Particle_flow_var(Track_struct track)
{
    isTrack = true;
    charge = sgn(track.q_p);
    label  = track.position_in_list;
    float theta = track.theta;
    eta = -1*log(tan(0.5*theta));
    phi = track.phiHelix;
    energy = track.energy_perigee;
    px = track.px_perigee;
    py = track.py_perigee;
    pz = track.pz_perigee;
    truth_link = 4;

}

Particle_flow_var::Particle_flow_var(Topo_clust topo_clust)
{
    isTrack = false;
    charge = topo_clust.charge;
    label = topo_clust.label;
    eta = topo_clust.eta_com;
    phi = topo_clust.phi_com;
    energy = topo_clust.total_energy;
    // float theta = 2*atan(exp(-1*eta));
    px = topo_clust.px;//energy*sin(theta)*cos(phi);
    py = topo_clust.py;//energy*sin(theta)*sin(phi);
    pz = topo_clust.pz;//energy*cos(theta);
    truth_link = topo_clust.truth_link;
}

Particle_flow_var::~Particle_flow_var()
{;}
