#include "Tracking_func.hh"

Tracking::Tracking(std::vector<FullTrajectoryInfo> FSPs, bool Smearing, std::vector<Track_struct> &Tracks_list, Geometry_definition Geometry)
{
    geometry = Geometry;
    smearing = Smearing;
    int num_FSPs = FSPs.size();
    theta_end_barrel = 2 * atan(exp(-1 * config_json_var.max_eta_barrel));
    theta_end_endcap = 2 * atan(exp(-1 * config_json_var.max_eta_endcap));
    NLayers = geometry.kNLayers;
    R_mid = geometry.layer_mid_radius_flatten;
    LayersPix = geometry.number_of_pixels_flatten;
    d_eta = geometry.layer_deta_flatten;
    d_phi = geometry.layer_dphi_flatten;
    for (int nParticle = 0; nParticle < num_FSPs; nParticle++)
    {
	Trajectory_finder(FSPs.at(nParticle), nParticle, Tracks_list);
    }
    sort(Tracks_list.begin(), Tracks_list.end(), sort_by_pt);
    int size_tracks_list = Tracks_list.size();
    for (int itrack = 0; itrack < size_tracks_list; itrack++)
    {	
        Tracks_list.at(itrack).position_in_list = itrack;
    }
}

// Extrapolated track for charged particle
void Tracking::Trajectory_finder(FullTrajectoryInfo FSP, int nParticle, std::vector<Track_struct> &Tracks_list)
{
    double Hmagnetic = config_json_var.fieldValue;
    double charge = FSP.fPDGCharge; //* 1 elementary charge
    if ((charge != 0))              //&& (abs(DefParticle->GetPDGEncoding())!=11)
    {
        int pdg = (int)FSP.fPDGCode;
        double px = FSP.fMomentum.x(); //* Mev
        double py = FSP.fMomentum.y(); //* Mev
        double pz = FSP.fMomentum.z(); //* Mev
        double energy = FSP.fEnergy;
        double verx = FSP.fVertexPosition.x();
        double very = FSP.fVertexPosition.y();
        double verz = FSP.fVertexPosition.z();
        Track_struct track(pdg, nParticle, energy,
                           FSP.fMass, FSP.fPDGCharge,
                           px, py, pz,
                           verx, very, verz, 
                           FSP.is_conversion_track);
        for (int lay = 0; lay < NLayers; lay++)
        {
            track.x_mid_layer.push_back(-1000000);
            track.y_mid_layer.push_back(-1000000);
            track.z_mid_layer.push_back(-1000000);
            track.eta.push_back(-1000000);
            track.phi.push_back(-1000000);
            track.ind_phi.push_back(-1000000);
            track.ind_eta.push_back(-1000000);
        }
	// Compute the Helix for charged particle in Magnetic field
	// Helix computed until first layer of ECAL1 (config_json_var.r_inn_calo)
        if (Hmagnetic != 0)
        {
            track.rho = fabs(track.pt / (track.charge * Hmagnetic * ligh_speed)); // R of curvature rho = pT/0.3qB
            float x0 = track.initX + track.charge * track.rho * (py / track.pt);
            float y0 = track.initY - track.charge * track.rho * (px / track.pt);
            float r0 = sqrt(sqr(x0) + sqr(y0));
            float a0sqr = sqr(r0 - track.rho);
            float R_of_MF_sqr = sqr(config_json_var.r_inn_calo);
            track.alpha = sqrt(((R_of_MF_sqr - a0sqr) / (4 * track.rho * r0)));

            //* Perigee parametrs
            track.a0 = track.charge * (sqrt(a0sqr));
            track.phiHelix = atan2(y0, x0) + charge * M_PI_2;
            track.phiHelix = track.phiHelix > M_PI ? track.phiHelix - 2 * M_PI : track.phiHelix < -M_PI ? track.phiHelix + 2 * M_PI : track.phiHelix;
            track.theta = acos(pz / sqrt(sqr(track.pt) + sqr(pz)));
            track.q_p = charge / sqrt(sqr(track.pt) + sqr(track.pz));
            //* numerical stability
            if ((sqr(track.initX) + sqr(track.initY) - a0sqr) < 0)
            {
                a0sqr = 0;
            }
            track.z0 = track.initZ - 2 * track.rho * (1 / tan(track.theta)) * asin(sqrt((sqr(track.initX) + sqr(track.initY) - a0sqr) / (4 * sqr(track.rho) + 4 * track.rho * charge * track.a0)));

            if (smearing)
                track.smearing();
            track.p_perigee = charge * (1 / track.q_p);
            track.pt_perigee = track.p_perigee * sin(track.theta);
            track.px_perigee = track.p_perigee * sin(track.theta) * cos(track.phiHelix);
            track.py_perigee = track.p_perigee * sin(track.theta) * sin(track.phiHelix);
            track.pz_perigee = track.p_perigee * cos(track.theta);
            track.mass_perigee = G4ParticleTable::GetParticleTable()->FindParticle("pi+")->GetPDGMass();
            track.energy_perigee = sqrt(sqr(track.p_perigee) + sqr(track.mass_perigee));
            track.rho_perigee = fabs(track.pt_perigee / (charge * Hmagnetic * ligh_speed));
            float param_time = 2 / track.pt_perigee * asin(sqrt((sqr(config_json_var.r_inn_calo) - a0sqr) / (4 * sqr(track.rho_perigee) + 4 * charge * track.rho_perigee * track.a0)));

	    // Helix until detector if not Endcap (equivalent of finding phi at detector position)
            track.x_end_MF = (charge * track.a0 + track.rho_perigee) * cos(track.phiHelix - charge * M_PI_2) + charge * track.rho_perigee * sin(charge * param_time * track.pt_perigee - track.phiHelix);
            track.y_end_MF = (charge * track.a0 + track.rho_perigee) * sin(track.phiHelix - charge * M_PI_2) + charge * track.rho_perigee * cos(charge * param_time * track.pt_perigee - track.phiHelix);

            track.z_end_MF = track.z0 + track.rho_perigee * track.pz_perigee * param_time;

            track.px_end_MF = track.pt_perigee * cos(charge * param_time * track.pt_perigee - track.phiHelix);
            track.py_end_MF = -track.pt_perigee * sin(charge * param_time * track.pt_perigee - track.phiHelix);

            float length = config_json_var.r_inn_calo / tan(theta_end_barrel);
            //* coeff of quadratic equation for traject extrapol
            track.a_coeff = sqr(track.pt_perigee);
            track.b_coeff = 2 * (track.px_end_MF * track.x_end_MF + track.py_end_MF * track.y_end_MF);
	    // Helix in Endcap
            if (std::isnan(track.x_end_MF) || fabs(track.z_end_MF) > length)
            {
                int sign = track.z_end_MF / fabs(track.z_end_MF);
                track.z_end_MF = sign * length;
                float R_min_ECAL = length * tan(theta_end_endcap);
                float time_to_endcap = (track.z_end_MF - track.z0) / (track.rho_perigee * track.pz_perigee);
                float parameter = charge * time_to_endcap * track.pt_perigee - track.phiHelix;
                track.x_end_MF = (charge * track.a0 + track.rho_perigee) * cos(track.phiHelix - charge * M_PI_2) + charge * track.rho_perigee * sin(parameter);
                track.y_end_MF = (charge * track.a0 + track.rho_perigee) * sin(track.phiHelix - charge * M_PI_2) + charge * track.rho_perigee * cos(parameter);
                track.px_end_MF = track.pt_perigee * cos(parameter);
                track.py_end_MF = -track.pt_perigee * sin(parameter);

                if ((sqr(track.x_end_MF) + sqr(track.y_end_MF) >= sqr(R_min_ECAL))) //&& ((parameter)<=M_PI_2))//
                {

                    track.is_in_endcap = true;
                }
                else
                {
                    track.theta = -10000000;
                    track.a0 = -10000000;
                    track.q_p = -10000000;
                    track.z0 = -10000000;
                    track.phiHelix = -10000000;
                    track.Is_reach_calorimeter = false;
                }
            }
	    // continue with Linear Extrapolation in Calorimeter if Magnetic field
            for (int lay = 0; lay < NLayers; lay++)
                TrajectoryInMF(track, lay);
        }
	// no Magnetic field
        else 
        {
	    // compute linear extrapolation
            track.a0 = -10000000;
            track.z0 = -10000000;
            track.phiHelix = atan2(track.py, track.px);;
            track.theta = acos(pz / sqrt(sqr(track.pt) + sqr(pz)));
            track.q_p = charge / sqrt(sqr(track.pt) + sqr(track.pz));
            for (int lay = 0; lay < NLayers; lay++)
            {
                float time_to_layer = R_mid.at(lay) / track.pt;
                float z_f = track.initZ + pz * time_to_layer;
                float x_f = track.initX + px * time_to_layer;
                float y_f = track.initY + py * time_to_layer;
                track.x_mid_layer.at(lay) = x_f;
                track.y_mid_layer.at(lay) = y_f;
                track.z_mid_layer.at(lay) = z_f;
                track.eta.at(lay) = -1 * log(tan(0.5 * acos(z_f / sqrt(sqr(R_mid.at(lay)) + sqr(z_f)))));
                track.phi.at(lay) = atan2(y_f, x_f);
                float deta = d_eta.at(lay); //!
                float dphi = d_phi.at(lay); //!
                track.ind_phi.at(lay) = (int)floor(track.phi.at(lay) / dphi);
                track.ind_eta.at(lay) = (int)floor((config_json_var.max_eta_endcap + track.eta.back()) / deta);
            }
		
        }
	// Apply track inefficiency as function of pT
        if (smearing)
            track.IsTrackReconstructed();

        if (config_json_var.Skip_unuseable_tracks) //Do not store bad tracks
        {
            if(track.Is_Track_Useable())
                Tracks_list.push_back(track);
        }
        else
            Tracks_list.push_back(track);

    }
}

void Tracking::TrajectoryInMF(Track_struct &track, int lay)
{
    float length = config_json_var.r_inn_calo / tan(theta_end_barrel);
    if (track.is_in_endcap)
    {
        int sign = track.z_end_MF / fabs(track.z_end_MF);
       
        float time_to_layer = (R_mid.at(lay) / tan(theta_end_barrel) * sign - track.z_end_MF) / (track.rho_perigee * track.pz_perigee);
        float X_in_mid_layer = track.x_end_MF + track.rho_perigee * track.px_end_MF * time_to_layer;
        float Y_in_mid_layer = track.y_end_MF + track.rho_perigee * track.py_end_MF * time_to_layer;
        float Z_in_mid_layer = track.z_end_MF + track.rho_perigee * track.pz_perigee * time_to_layer;
        float Phi = atan2(Y_in_mid_layer, X_in_mid_layer);
        float Eta = -1 * log(tan(0.5 * acos(Z_in_mid_layer / sqrt(sqr(X_in_mid_layer) + sqr(Y_in_mid_layer) + sqr(Z_in_mid_layer)))));

        float deta = d_eta.at(lay);
        float dphi = d_phi.at(lay);
        track.x_mid_layer.at(lay) = X_in_mid_layer;
        track.y_mid_layer.at(lay) = Y_in_mid_layer;
        track.z_mid_layer.at(lay) = Z_in_mid_layer;

        track.eta.at(lay) = Eta;
        track.phi.at(lay) = Phi;
        track.ind_eta.at(lay) = floor((config_json_var.max_eta_endcap + Eta) / deta);
        track.ind_phi.at(lay) = floor(Phi / dphi);
    }
    else if (!std::isnan(track.x_end_MF) && track.z_end_MF < length)
    {
        float c_coeff = sqr(config_json_var.r_inn_calo) - sqr(R_mid.at(lay));
        float time_to_layer = (-track.b_coeff + sqrt(sqr(track.b_coeff) - 4 * track.a_coeff * c_coeff)) / (2 * track.a_coeff);
        float X_in_mid_layer = track.x_end_MF + track.px_end_MF * time_to_layer;
        float Y_in_mid_layer = track.y_end_MF + track.py_end_MF * time_to_layer;
        float Z_in_mid_layer = track.z_end_MF + track.pz_perigee * time_to_layer;
        float Phi = atan2(Y_in_mid_layer, X_in_mid_layer);
        float Eta = -1 * log(tan(0.5 * acos(Z_in_mid_layer / sqrt( sqr(R_mid.at(lay)) + sqr(Z_in_mid_layer) ) )));
	
        float deta = d_eta.at(lay);
        float dphi = d_phi.at(lay);
        track.x_mid_layer.at(lay) = X_in_mid_layer;
        track.y_mid_layer.at(lay) = Y_in_mid_layer;
        track.z_mid_layer.at(lay) = Z_in_mid_layer;
        track.eta.at(lay) = Eta;
        track.phi.at(lay) = Phi;
        track.ind_eta.at(lay) = floor((config_json_var.max_eta_endcap + Eta) / deta);
        track.ind_phi.at(lay) = floor(Phi / dphi);
    }
    else
    {
        track.Is_reach_calorimeter = false;
    }
}
