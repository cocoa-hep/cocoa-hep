#ifndef __TRUTHRECORDGRAPH_H__
#define __TRUTHRECORDGRAPH_H__



#include "HepMC3/GenEvent.h"
#include "Config_reader_var.hh"
#include "TTree.h"
#include <vector>
#include <stdio.h>
#include <iostream>

class OutputRunAction;

class TruthRecordGraph {

public:
    Config_reader_var &config = Config_reader_var::GetInstance();
    std::vector<HepMC3::ConstGenParticlePtr> m_interesting_particles;
    std::vector<HepMC3::ConstGenParticlePtr> m_final_state_particles;
    float m_max_radius = config.r_inn_calo;
    static TruthRecordGraph &GetInstance()
    {
        static TruthRecordGraph fsp;
        return fsp;
    };
    void clear();
    void fill_truth_graph();
    void set_tree_branches(TTree *outTree);

    bool isCharmHadron(int pdgid);
    bool isBottomHadron(int pdgid);
    void add_to_vector(HepMC3::ConstGenParticlePtr to_add, std::vector<HepMC3::ConstGenParticlePtr > &add_to);
    void add_all_moving_parents(HepMC3::ConstGenParticlePtr to_add, std::vector<HepMC3::ConstGenParticlePtr > &add_to);
    bool is_parent_same(HepMC3::ConstGenParticlePtr particle);
    bool is_parent_of(HepMC3::ConstGenParticlePtr parent, HepMC3::ConstGenParticlePtr potential_child);
    void clean_daugthers(   std::vector<HepMC3::ConstGenParticlePtr > &interesting_particles,
                            std::vector<int> &all_daughters, std::vector<int> &direct_daughters);
    void find_daughters(HepMC3::ConstGenParticlePtr parent, std::vector<HepMC3::ConstGenParticlePtr > &interesting_particles, std::vector<int> &direct_daughters);
    HepMC3::ConstGenParticlePtr check_prod_location(HepMC3::ConstGenParticlePtr particle);
    HepMC3::ConstGenParticlePtr find_next_level_parent(HepMC3::ConstGenParticlePtr particle);
    float prod_radius(HepMC3::ConstGenParticlePtr particle);

    std::vector<int> CharmHadrons = {   411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435,
                                        4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414,
                                        4424, 4432, 4434, 4444, 441, 10441, 100441, 443, 10443, 20443, 100443, 30443, 9000443, 9010443, 9020443, 445, 100445};

    std::vector<int> BottomHadrons = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 2515, 2525, 531, 10531, 533, 10533, 20533, 2535, 541, 10541, 543, 10543, 20543, 2545};

private:
    std::vector<int> node_idx, node_pdg_id;
    std::vector<float> node_phi, node_eta, node_pt, node_m;
    std::vector<float> node_prodx, node_prody, node_prodz;
    std::vector<float> node_decx, node_decy, node_decz;
    std::vector<int> node_isfinal;
    std::vector<int> edge_start, edge_end;
    std::vector<int> final_idx;
};
#endif // __TRUTHRECORDGRAPH_H__
