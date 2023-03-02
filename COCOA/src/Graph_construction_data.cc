#include "Graph_construction_data.hh"
#include "Config_reader_var.hh"

Graph_construction_data::Graph_construction_data(/* args */)
{
}

Graph_construction_data::~Graph_construction_data()
{
}

void Graph_construction_data::clear()
{
    track_to_cell_edge_start.clear();
    track_to_cell_edge_end.clear();
    cell_to_cell_edge_start.clear();
    cell_to_cell_edge_end.clear();
    particle_to_node_idx.clear();
    particle_to_node_weight.clear();
    low_cell_to_high_cell_edge.clear();
    high_cell_to_low_cell_edge.clear();
}

void Graph_construction_data::set_tree_branches(TTree *outTree)
{
    Config_reader_var &config_var = Config_reader_var::GetInstance();
    if ( config_var.Use_high_granularity )
    {
        outTree->Branch("low_cell_to_high_cell_edge", "vector<int>", &low_cell_to_high_cell_edge);
        outTree->Branch("high_cell_to_low_cell_edge", "vector<int>", &high_cell_to_low_cell_edge);
    }
    outTree->Branch("track_to_cell_edge_start", "vector<int>", &track_to_cell_edge_start);
    outTree->Branch("track_to_cell_edge_end", "vector<int>", &track_to_cell_edge_end);
    outTree->Branch("cell_to_cell_edge_start", "vector<int>", &cell_to_cell_edge_start);
    outTree->Branch("cell_to_cell_edge_end", "vector<int>", &cell_to_cell_edge_end);
    outTree->Branch("particle_to_node_idx", "vector<vector<float>>", &particle_to_node_idx);
    outTree->Branch("particle_to_node_weight", "vector<vector<float>>", &particle_to_node_weight);
}