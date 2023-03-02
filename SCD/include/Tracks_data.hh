#ifndef __TRACKS_DATA_H__
#define __TRACKS_DATA_H__

#include <vector>
#include "Track_var.hh"
#include "TTree.h"
#include "Config_reader_var.hh"
#include "Full_trajectory_info_data.hh"



class Tracks_data
{
    public:
        Tracks_data();
        // ~Tracks_data();
        void Clear();
        void set_tree_branches(TTree *outTree, int NLayers, std::string Type_of_running = "");
        void Fill_perigee_var();
        static Tracks_data &GetHigh() 
        {
            static Tracks_data track_list_high; 
            return track_list_high;
        };
        static Tracks_data &GetLow() 
        {
            static Tracks_data track_list_low; 
            return track_list_low; 
        };
        std::vector <int> track_pflow_object_idx;
        std::vector <Track_struct> Tracks_list;
        std::vector <float> PerigeeA0;
        std::vector <float> PerigeeZ0;
        std::vector <float> PerigeeTheta;
        std::vector <float> PerigeePhi;
        std::vector <float> PerigeeQ_P;
        std::vector <int>   track_reconstructed;
        std::vector <int>   track_accepted;
        std::vector <int>   TrckPDGID;
        std::vector <int>   TrckPosInRealList;
        std::map<TString, std::vector<float>* > track_extrap_branches;
        std::vector<int>    track_LHED;

    private:
        Config_reader_var& config_json_var = Config_reader_var::GetInstance();
        int nLayers;
};


#endif // __TRACKS_DATA_H__
