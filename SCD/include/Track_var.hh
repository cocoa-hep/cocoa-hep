#ifndef __TRACK_VAR_H__
#define __TRACK_VAR_H__

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "Config_reader_var.hh"
#include "DetectorGeometryDefinitions.hh"
#include "CSVReader.hh"
// #include "DataStorage.hh"

class Track_struct
{

    // private:
    // Config_reader_var& config_json_var = Config_reader_var::GetInstance();
public:
    Track_struct();
    Track_struct(const Track_struct &orig);
    Track_struct(int Pdgcode, int NFinalStateParticles, double Energy,
                 double Mass, double Charge,
                 double Px, double Py, double Pz,
                 double InitX, double InitY, double InitZ,
                 bool _is_conversion_track=false);
    int pdgcode;                //code of particle
    int nFinal_State_Particles; //where particle in the list of FinalStateParticles
    double energy;
    double px;
    double py;
    double pz;
    double absmom;
    std::vector<double> eta;         //eta coordinate of Particle in each layer
    std::vector<double> phi;         //phi coordinate of Particle in each layer
    std::vector<double> x_mid_layer; //eta coordinate of Particle in each layer
    std::vector<double> y_mid_layer; //phi coordinate of Particle in each layer
    std::vector<double> z_mid_layer; //eta coordinate of Particle in each layer
    std::vector<int> ind_eta;        //index of cell in eta direction of Particle in each layer
    std::vector<int> ind_phi;        //index of cell in phi direction of Particle in each layer
    std::vector<double> Rprime_to_closest_topoclusters;
    std::vector<int> index_of_closest_topoclusters;
    double a0;
    double z0;
    double initX;
    double initY;
    double initZ;
    double theta;
    double phiHelix;
    double q_p;
    double mass;
    double pt;
    int position_in_list;

    bool Is_track_reconstracted;
    bool Is_inside_R;
    bool Is_reach_calorimeter;
    void IsProductInsideRadius();
    void IsTrackReconstructed();
    bool Is_Track_Useable();

    std::vector<std::vector<double>> sigmaA0;
    std::vector<std::vector<double>> sigmaZ0;
    std::vector<std::vector<double>> sigmaQP;
    std::vector<std::vector<double>> sigmaTHETA;
    std::vector<std::vector<double>> sigmaPHI;
    void smearing();
    // double sigmafind(std::vector < std::vector <double> > sigmavect, double ptvalue);

    double charge;
    double alpha;
    double px_end_MF;
    double py_end_MF;
    double pz_perigee;
    double p_perigee;
    double px_perigee;
    double py_perigee;
    double pt_perigee;
    double energy_perigee;
    double mass_perigee;
    double rho_perigee;
    double x_end_MF;
    double y_end_MF;
    double z_end_MF;
    double R_of_megfield_sqr;
    double a_coeff;
    double b_coeff;
    double theta_end_barrel;
    double theta_end_endcap;
    bool is_in_endcap;
    bool is_conversion_track;
    double rho;

    int  GetLHED() const { return m_lhed; }
    void SetLHED( int lhed ) { m_lhed = lhed; }

private:
    int m_lhed; // calo layer of highest energy density increase around the track
    
};

#endif // __TRACK_VAR_H__
