#ifndef __PARTICLE_FLOW_DATA_H__
#define __PARTICLE_FLOW_DATA_H__

#include "Particle_flow_var.hh"
#include "TTree.h"
#include "fastjet/ClusterSequence.hh"



class Particle_flow_data
{
private:
    std::vector<int> cell_pflow_object_idx;
    std::vector<float> pflow_eta;
    std::vector<float> pflow_phi;
    std::vector<float> pflow_px;
    std::vector<float> pflow_py;
    std::vector<float> pflow_pz;
    std::vector<float> pflow_e;
    std::vector<float> pflow_charge;
    std::vector<int> pflow_parent_idx;
public:
    Particle_flow_data();
    ~Particle_flow_data();
    std::vector<Particle_flow_var> pflow_list;
    static Particle_flow_data &GetInstance()
    {
        static Particle_flow_data particel_flow_data;
        return particel_flow_data;
    };
    void make_pseudo_jet_particles();
    void clear();
    void set_tree_branches(TTree *outTree);
    void fill_cell_var();
    std::vector <fastjet::PseudoJet> jets_objects;
};


#endif // __PARTICLE_FLOW_DATA_H__