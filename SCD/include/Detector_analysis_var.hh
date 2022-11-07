#ifndef __DETECTOR_ANALYSIS_VAR_H__
#define __DETECTOR_ANALYSIS_VAR_H__

#include <vector>
#include "Config_reader_var.hh"
#include "TTree.h"
#include "Cell_var.hh"

class Detector_analysis_var
{
private:
    float eta = 0;
    std::vector<float> rad_length;
    std::vector<float> int_length;
    std::string Type_Of_Running;
    float charge_transverse_energy = 0;
    float neutral_transverse_energy = 0;
    float neutral_energy = 0;
    float charge_energy = 0;
    Geometry_definition Geometry;
public:
    Detector_analysis_var();
    static Detector_analysis_var &GetInstance()
    {
        static Detector_analysis_var det_ana;
        return det_ana;
    };
    void set_eta(float Eta);
    float get_eta();
    void add_lengths(float x0, float lambda, int layer);
    void set_tree_branches(TTree *outTree, std::string type_of_running);
    void transverse_energy_calculation(std::vector<std::vector<std::vector<Cell>>> fcell_array);
    void clear();
};

#endif // __DETECTOR_ANALYSIS_VAR_H__