#ifndef __TRACKING_FUNC_H__
#define __TRACKING_FUNC_H__

#include "Config_reader_var.hh"
#include "Track_var.hh"
#include "G4SystemOfUnits.hh"
#include "Tracks_data.hh"
#include "Full_trajectory_info_data.hh"
#include "algorithm"

inline bool sort_by_pt(const Track_struct &track_1, const Track_struct &track_2)
{
    return track_1.pt > track_2.pt;
}

class Tracking
{
    public:
        Tracking(std::vector<FullTrajectoryInfo>    FSPs, bool smearing,std::vector<Track_struct> &Tracks_list, Geometry_definition Geometry);
    private:
        bool smearing;
        Config_reader_var& config_json_var = Config_reader_var::GetInstance();
        // Tracks_data& tracks_list = Tracks_data::GetHigh();
        // Tracks_data& tracks_list = Tracks_data::GetLow();
        void Trajectory_finder(FullTrajectoryInfo FSP, int nParticle, std::vector<Track_struct> &Tracks_list);
        // void TrajectoryInMF(Track_struct track, float endR, int i);
        void TrajectoryInMF(Track_struct &track, int lay) ;
        float ligh_speed = 299792458 * m/s;
        double theta_end_barrel;
        double theta_end_endcap;
        std::vector <long double> R_mid;
        int NLayers;
        std::vector <int> LayersPix;
        std::vector<long double> d_eta;
        std::vector<long double> d_phi;
        Geometry_definition geometry;
};
#endif // __TRACKING_FUNC_H__