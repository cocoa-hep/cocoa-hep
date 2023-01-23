#ifndef __GRAPH_CONSTRUCTION_DATA_H__
#define __GRAPH_CONSTRUCTION_DATA_H__

#include "TTree.h"

class Graph_construction_data
{
private:
    
public:
    Graph_construction_data(/* args */);
    ~Graph_construction_data();
    static Graph_construction_data &GetLow()
    {
        static Graph_construction_data graph_low;
        return graph_low;
    };

    static Graph_construction_data &GetHigh()
    {
        static Graph_construction_data graph_high;
        return graph_high;
    };
    void set_tree_branches(TTree *outTree);
    // void fill_graph_var();
    void clear();
    std::vector<int>                track_to_cell_edge_start;
    std::vector<int>                track_to_cell_edge_end;
    std::vector<int>                cell_to_cell_edge_start;
    std::vector<int>                cell_to_cell_edge_end;
    std::vector<std::vector<float>> particle_to_node_idx;
    std::vector<std::vector<float>> particle_to_node_weight;
    std::vector<int> low_cell_to_high_cell_edge;
    std::vector<int> high_cell_to_low_cell_edge;
};


#endif // __GRAPH_CONSTRUCTION_DATA_H__