#ifndef __JET_BUILDER_DATA_H__
#define __JET_BUILDER_DATA_H__

#include <vector>
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

class Jet_Builder_data
{
private:
    std::vector<float> jet_pt;
    std::vector<float> jet_eta;
    std::vector<float> jet_phi;
    std::vector<float> jet_m;
    std::vector<int>   jet_constituents_jetIndex; // jet index for each constituent. default -1: not associated to a jet.
    unsigned int       n_constituents;
    std::vector<int>   jet_truth_id; // pdg id of the truth level particle associated with the jet
public:
  
    Jet_Builder_data(std::string prefix);
    ~Jet_Builder_data();
    void        clear();
    void        set_tree_branches(TTree *outTree);
    inline void set_n_constituents( unsigned int n ) { n_constituents = n; }
    void        fill_cell_var( float truth_jet_radius = -1.0 );
    
    static Jet_Builder_data &Get_instance_pflow()
    {
        static Jet_Builder_data jets_data("pflow");
        return jets_data;
    };
    static Jet_Builder_data &Get_instance_true()
    {
        static Jet_Builder_data jets_data("true");
        return jets_data;
    };
    static Jet_Builder_data &Get_instance_topo()
    {
        static Jet_Builder_data jets_data("topo");
        return jets_data;
    };
    std::vector<fastjet::PseudoJet> jets;
    std::string Prefix;

};


#endif // __JET_BUILDER_DATA_H__
