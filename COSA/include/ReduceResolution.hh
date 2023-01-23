#ifndef __DOWNRESOLUTION_H__
#define __DOWNRESOLUTION_H__

// #include "DataStorage.hh"
// #include "Cell_interface.hh"
#include "Cells_data.hh"
// #include "DetectorGeometryDefinitions.hh"
#include "Config_reader_var.hh"
#include "TRandom.h"


class ReduceResolution
{
public:

    ReduceResolution(std::vector<std::vector<std::vector<Cell>>> &CellArray_High, std::vector<std::vector<std::vector<Cell>>> &CellArray_Low);
    ReduceResolution();
    ~ReduceResolution();
    void link_apply(std::vector<std::vector<std::vector<Cell>>> &CellArray_High, std::vector<std::vector<std::vector<Cell>>> &CellArray_Low);
    void link_apply(std::vector<Cell> &high_cells_in_clusts, std::vector<Cell> &low_cells_in_clusts);
    void apply_noise(std::vector<std::vector<std::vector<Cell>>> &CellArray);

private:
    int init_low;
    int init_high;
    Config_reader_var &config_var = Config_reader_var::GetInstance();
    void sum_pixels(std::vector<std::vector<std::vector<Cell>>> &CellArray_High, std::vector<std::vector<std::vector<Cell>>> &CellArray_Low,
                    std::vector<std::vector<int>> Low_pixel, std::vector<std::vector<int>> High_pixel);
    void topo_label_apply(std::vector<std::vector<std::vector<Cell>>> &CellArray_High, std::vector<std::vector<std::vector<Cell>>> &CellArray_Low,
                          std::vector<std::vector<int>> Low_pixel, std::vector<std::vector<int>> High_pixel);
};
#endif // __DOWNRESOLUTION_H__