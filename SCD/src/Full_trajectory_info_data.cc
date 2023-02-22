#include "Full_trajectory_info_data.hh"
#include "DetectorConstruction.hh"
#include <algorithm>
#include <cassert>

Full_trajectory_info_data::Full_trajectory_info_data(/* args */)
{;}


void Full_trajectory_info_data::make_pseudo_jet_particles()
{
    int size_trajectories_list = fAllTrajectoryInfo.size();
    for (int iparticle = 0; iparticle < size_trajectories_list; iparticle++)
    {
        fastjet::PseudoJet particle(fAllTrajectoryInfo.at(iparticle).fMomentum.x(),
                                                  fAllTrajectoryInfo.at(iparticle).fMomentum.y(),
                                                  fAllTrajectoryInfo.at(iparticle).fMomentum.z(),
                                                  fAllTrajectoryInfo.at(iparticle).fEnergy);
        particle.set_user_index(iparticle);
        jets_objects.push_back(particle);
    }
}

int Full_trajectory_info_data::DeltaR_iso(float px, float py, float pz ,size_t idx_m,int pdgid) {

    size_t n_particles = fAllTrajectoryInfo.size();
    
    float ptot_origin = sqrt(px*px+py*py+pz*pz);
    float theta_origin = acos(pz/ptot_origin);
    float eta_origin = -1*log(tan(theta_origin/2.));
    float phi_origin = atan(py/px);

    float pT =  sqrt(px*px+py*py) / 1e3;
    float iso_radius = std::min(10./pT, 0.2);

    int isIso = 0;
    //consider all neutral to be non-isolated
    if(fabs(pdgid) == 22 || fabs(pdgid) == 310 || fabs(pdgid) == 3122 || fabs(pdgid) == 3322 || fabs(pdgid) == 2112) return 0;

    float sum_pt = 0;
    for(size_t part_i=0; part_i < n_particles; part_i++){
	if (idx_m == part_i) continue;
	int pdg_id_loop = fAllTrajectoryInfo.at(part_i).fPDGCode;
         //exclude from the computation neutral particles
         if(fabs(pdg_id_loop) == 22 || fabs(pdg_id_loop) == 310 || fabs(pdg_id_loop) == 3122 || fabs(pdg_id_loop) == 3322 || fabs(pdg_id_loop) == 2112) continue;
	float px_loop   = fAllTrajectoryInfo.at(part_i).fMomentum.x();
	float py_loop   = fAllTrajectoryInfo.at(part_i).fMomentum.y();
	float pz_loop   = fAllTrajectoryInfo.at(part_i).fMomentum.z();
	float phi       = atan(py_loop/px_loop);
	float ptot_loop = sqrt(px_loop*px_loop+py_loop*py_loop+pz_loop*pz_loop);
	float theta     = acos(pz_loop/ptot_loop);
	float eta       = -1*log(tan(theta/2.));
	float dphi      = acos(cos(phi_origin-phi));
	float deta      = eta_origin-eta;
	float dr        = sqrt(dphi*dphi+deta*deta);
	if(dr<iso_radius){
	    sum_pt = sqrt(px_loop*px_loop+py_loop*py_loop) /1e3 + sum_pt;
	}
    }
    if(sum_pt/pT < 0.06) isIso = 1;
    return isIso;
    
}

void Full_trajectory_info_data::fill_var(){
    int size_trajectories_list = fAllTrajectoryInfo.size();
    for (int iparticle = 0; iparticle < size_trajectories_list; iparticle++)
    {
        float px    = fAllTrajectoryInfo.at(iparticle).fMomentum.x();
        float py    = fAllTrajectoryInfo.at(iparticle).fMomentum.y();
        float pz    = fAllTrajectoryInfo.at(iparticle).fMomentum.z();
        float p     = sqrt( px * px + py * py + pz * pz ); // Have some 4-momentum object at some point (but would like to circumvent ROOT). PR.
        float theta = acos( pz / p );
        float eta   = -1.0 * log( tan( 0.5 * theta ) );
        float phi   = GetPhi( px, py );
	
        particle_pdgid.push_back(fAllTrajectoryInfo.at(iparticle).fPDGCode);
	    particleisIso.push_back(DeltaR_iso(px, py, pz,iparticle,int(fAllTrajectoryInfo.at(iparticle).fPDGCode)));
        particle_pt.push_back(sqrtf(px*px + py*py));
        particle_eta.push_back(eta);
        particle_phi.push_back(phi);
        particle_e.push_back(fAllTrajectoryInfo.at(iparticle).fEnergy);
        particle_prod_x.push_back(fAllTrajectoryInfo.at(iparticle).fVertexPosition.x());
        particle_prod_y.push_back(fAllTrajectoryInfo.at(iparticle).fVertexPosition.y());
        particle_prod_z.push_back(fAllTrajectoryInfo.at(iparticle).fVertexPosition.z());
	
        caloExtrapolEta.push_back(fAllTrajectoryInfo.at(iparticle).caloExtrapolEta);
        caloExtrapolPhi.push_back(fAllTrajectoryInfo.at(iparticle).caloExtrapolPhi);
	
        idExtrapolEta.push_back(fAllTrajectoryInfo.at(iparticle).idExtrapolEta);
        idExtrapolPhi.push_back(fAllTrajectoryInfo.at(iparticle).idExtrapolPhi);
	
        particle_to_track.push_back(-1);	
    }

    Tracks_data &tracks_obj = Tracks_data::GetLow();
    int size_tracks = tracks_obj.Tracks_list.size();
    for (int itrack = 0; itrack < size_tracks; itrack++)
    {
	    int trackParticleId = tracks_obj.Tracks_list[itrack].nFinal_State_Particles;
        particle_to_track.at(trackParticleId) = itrack;	
    }
    
    for( const FullTrajectoryInfo& conv_el_tr : fAllConvElectrons ) {
        conv_el_fPrimaryPhotonIndex.push_back( conv_el_tr.fParentID );
        conv_el_q.push_back( conv_el_tr.fPDGCharge );
        conv_el_px.push_back( conv_el_tr.fMomentum.x() );
        conv_el_py.push_back( conv_el_tr.fMomentum.y() );
        conv_el_pz.push_back( conv_el_tr.fMomentum.z() );
        conv_el_prod_x.push_back( conv_el_tr.fVertexPosition.x() );
        conv_el_prod_y.push_back( conv_el_tr.fVertexPosition.y() );
        conv_el_prod_z.push_back( conv_el_tr.fVertexPosition.z() );
    }
    
}



void Full_trajectory_info_data::clear()
{
    fAllTrajectoryInfo.clear();
    fAllConvElectrons.clear();
    particle_pdgid.clear();
    particleisIso.clear();
    particle_pt.clear();
    particle_eta.clear();
    particle_phi.clear();
    particle_e.clear();
    particle_prod_x.clear();
    particle_prod_y.clear();
    particle_prod_z.clear();
    conv_el_fPrimaryPhotonIndex.clear();
    conv_el_q.clear();
    conv_el_px.clear();
    conv_el_py.clear();
    conv_el_pz.clear();
    conv_el_prod_x.clear();
    conv_el_prod_y.clear();
    conv_el_prod_z.clear();
    caloExtrapolEta.clear();
    caloExtrapolPhi.clear();
    idExtrapolEta.clear();
    idExtrapolPhi.clear();
    particle_to_track.clear();
    jets_objects.clear();
    particle_dep_energy.clear();
}

void Full_trajectory_info_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch("particle_pdgid"   ,           "vector<int>", &particle_pdgid);
    outTree->Branch("particle_isIso"   ,           "vector<int>", &particleisIso);
    outTree->Branch("particle_pt"  ,               "vector<float>", &particle_pt);
    outTree->Branch("particle_eta"  ,              "vector<float>", &particle_eta);
    outTree->Branch("particle_phi"  ,              "vector<float>", &particle_phi);
    outTree->Branch("particle_e"   ,               "vector<float>", &particle_e);
    outTree->Branch("particle_prod_x"  ,           "vector<float>", &particle_prod_x);
    outTree->Branch("particle_prod_y"  ,           "vector<float>", &particle_prod_y);
    outTree->Branch("particle_prod_z"  ,           "vector<float>", &particle_prod_z);
    outTree->Branch("particle_track_idx",          "vector<int>", &particle_to_track);
    outTree->Branch("particle_eta_extrap_calo",    "vector<float>", &caloExtrapolEta);
    outTree->Branch("particle_phi_extrap_calo",    "vector<float>", &caloExtrapolPhi);
    outTree->Branch("particle_eta_extrap_its",     "vector<float>", &idExtrapolEta);
    outTree->Branch("particle_phi_extrap_its",     "vector<float>", &idExtrapolPhi);
    outTree->Branch("particle_dep_energy"  ,       "vector<float>", &particle_dep_energy);

    outTree->Branch("conv_el_primary_photon_idx",  "vector<int>",   &conv_el_fPrimaryPhotonIndex);
    outTree->Branch("conv_el_q",                   "vector<float>", &conv_el_q);
    outTree->Branch("conv_el_px",                  "vector<float>", &conv_el_px);
    outTree->Branch("conv_el_py",                  "vector<float>", &conv_el_py);
    outTree->Branch("conv_el_pz",                  "vector<float>", &conv_el_pz);
    outTree->Branch("conv_el_prod_x",              "vector<float>", &conv_el_prod_x);
    outTree->Branch("conv_el_prod_y",              "vector<float>", &conv_el_prod_y);
    outTree->Branch("conv_el_prod_z",              "vector<float>", &conv_el_prod_z);
}

void Full_trajectory_info_data::SetParticleDepEnergy( const std::vector<float> &_particle_dep_energies) {

    particle_dep_energy = _particle_dep_energies;
    
}

