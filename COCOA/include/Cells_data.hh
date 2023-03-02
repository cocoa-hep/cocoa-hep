#ifndef __CELLS_DATA_H__
#define __CELLS_DATA_H__

#include "Config_reader_var.hh"
#include "Cell_var.hh"
#include "TTree.h"
#include "Full_trajectory_info_data.hh"
#include "Particle_flow_data.hh"

class Cells_data
{
private:
    bool high;
    std::vector<int> Number_Pixel_Flatten;
    Config_reader_var &config_var = Config_reader_var::GetInstance();
    std::vector<int> cell_pflow_object_idx;
    std::vector<int> cell_layer;
    std::vector<float> cell_x;
    std::vector<float> cell_y;
    std::vector<float> cell_z;
    std::vector<float> cell_eta;
    std::vector<float> cell_phi;
    std::vector<float> cell_e;
    std::vector<float> cell_che;
    std::vector<float> cell_nue;
    std::vector<int> cell_topo_idx;
    std::vector<int> cell_parent_idx;
    std::vector<int> cell_conv_el_idx;
    std::vector<std::vector<float>> cell_parent_list;
    std::vector<std::vector<float>> cell_parent_energy;

public:
    Cells_data(bool is_high);

    static Cells_data &GetHigh()
    {
        static Cells_data cells_high(true);
        return cells_high;
    };
    static Cells_data &GetLow()
    {
        static Cells_data cells_low(false);
        return cells_low;
    };

    const std::vector<int>& GetConvElIdx() { return cell_conv_el_idx; }
    
    std::vector<std::vector<std::vector<Cell>>> fCell_array;
    std::vector<Cell*> Cells_in_topoclusters;
    void ChangeLabelForCells(int OldLabel, int NewLabel);
    void set_tree_branches(TTree *outTree);
    void fill_cell_var();
    void add_cell_info(int ilay, int ieta, int iphi, float ch_en, float nu_en, Particle_dep_in_cell particle, Particle_dep_in_cell* conv_el = nullptr);
    void clear();
    void fill_cells_in_topoclusters();
    int ilay_num;
    // std::vector< std::vector< std::vector<Cell > > > // fLow_cell_array;
};
#endif // __CELLS_DATA_H__
