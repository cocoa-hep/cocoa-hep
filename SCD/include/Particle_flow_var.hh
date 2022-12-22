#ifndef __PARTICLE_FLOW_VAR_H__
#define __PARTICLE_FLOW_VAR_H__

#include "TTree.h"
#include "Topo_clust_var.hh"
#include "Track_var.hh"

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
class Particle_flow_var
{
private:
    /* data */
public:
    Particle_flow_var();
    Particle_flow_var(Track_struct track);
    Particle_flow_var(Topo_clust topo_clust);
    ~Particle_flow_var();
    void subt_cell(Cell &cell, float fraction);
    //if not isTrack && charge is not 0, then PFlow is TC from electron
    //if not isTrack && charge is 0, then PFlow is TC 
    //else isTrack, then PFlow is Track
    //if not isTrack, then label coresonds to topo_clust, if not then to track
    bool isTrack;
    float charge;
    int label;
    float eta;
    float phi;
    float px;
    float py;
    float pz;
    float energy;
    int truth_link;

};

#endif // __PARTICLE_FLOW_VAR_H__
