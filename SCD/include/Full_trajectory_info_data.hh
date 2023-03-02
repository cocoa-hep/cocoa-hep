#ifndef __FULL_TRAJECTORY_INFO_VAR_H__
#define __FULL_TRAJECTORY_INFO_VAR_H__

#include "vector"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "TTree.h"
#include "Tracks_data.hh"
#include "Track_var.hh"
#include "fastjet/ClusterSequence.hh"


struct FullTrajectoryInfo {
	G4int fTrackID;
	G4int fParentID;
	G4int fPDGCode;
	G4int fprocessId;
	G4double fPDGCharge;
	G4ThreeVector fMomentum;
	G4ThreeVector fMomentumDir;
	G4double    fEnergy ;
    G4double    fMass ;
	G4ThreeVector fVertexPosition;
	G4double fGlobalTime;
        
    float caloExtrapolMaxEkin;
    float caloExtrapolEta;
    float caloExtrapolPhi;
    
    float idExtrapolMaxEkin;
    float idExtrapolEta;
    float idExtrapolPhi;
        
	std::vector<G4int> vTrackID;
	std::vector<G4int> vParentID;
	std::vector<G4ThreeVector> vTrackMomentumDir;
	std::vector<G4ThreeVector> vTrackPos;
	std::vector<G4double> vTrackTime;
	std::vector<G4double> vTrackPDGID;

    bool is_conversion_track;

	// fastjet::PseudoJet particle_object;
};

class Full_trajectory_info_data
{
public:
    Full_trajectory_info_data();
	void make_pseudo_jet_particles();
    void clear();
void fill_var();
	void set_tree_branches(TTree *outTree);
    static Full_trajectory_info_data &GetInstance()
    {
        static Full_trajectory_info_data trajectories;
        return trajectories;
    };
    std::vector < FullTrajectoryInfo> fAllTrajectoryInfo;
    std::vector < FullTrajectoryInfo> fAllConvElectrons;
	std::vector <int>	particle_to_track;
	std::vector <fastjet::PseudoJet> jets_objects;
    void SetParticleDepEnergy( const std::vector<float> &_particle_dep_energies);
    int  DeltaR_iso(float px, float py, float pz ,size_t idx_m,int particle_loop_pdgid);
private:
	std::vector <int>	particle_pdgid;
        std::vector <int>   particleisIso;
	std::vector <float> particle_pt;
	std::vector <float> particle_eta;
	std::vector <float> particle_phi;
	std::vector <float> particle_e;
	std::vector <float> particle_prod_x;
	std::vector <float> particle_prod_y;
	std::vector <float> particle_prod_z;
    std::vector<float>  particle_dep_energy;

    std::vector<int>    conv_el_fPrimaryPhotonIndex;
    std::vector<float>  conv_el_q;
    std::vector<float>  conv_el_px;
    std::vector<float>  conv_el_py;
    std::vector<float>  conv_el_pz;
    std::vector<float>  conv_el_prod_x;
    std::vector<float>  conv_el_prod_y;
    std::vector<float>  conv_el_prod_z;
    
    std::vector<float> caloExtrapolEta;     // pseudorapidity of the same-pdgid daughter track of maximum kinetic energy in the calorimeter
    std::vector<float> caloExtrapolPhi;     // phi of the same-pdgid daughter track of maximum kinetic energy in the calorimeter
    std::vector<float> idExtrapolEta;     // pseudorapidity of the same-pdgid daughter track of maximum kinetic energy in the outermost id layer
    std::vector<float> idExtrapolPhi;     // phi of the same-pdgid daughter track of maximum kinetic energy in the outermost id layer

};

#endif // __FULL_TRAJECTORY_INFO_VAR_H__
