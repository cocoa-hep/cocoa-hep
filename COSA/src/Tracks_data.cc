#include "Tracks_data.hh"
#include "Config_reader_var.hh"
#include "Particle_flow_data.hh"

Tracks_data::Tracks_data()
{
    nLayers = 0;
    Tracks_list.clear();
}

void Tracks_data::Fill_perigee_var()
{
    Particle_flow_data &pflow_pbj = Particle_flow_data::GetInstance();
    int size_pflow = pflow_pbj.pflow_list.size();
    if (Tracks_list.size()!=0)
    {
        for (int ilay = 0; ilay < nLayers; ilay++)
        {
            track_extrap_branches["track_eta_layer_" + std::to_string(ilay)]->clear();
            track_extrap_branches["track_phi_layer_" + std::to_string(ilay)]->clear();
            track_extrap_branches["track_x_layer_" + std::to_string(ilay)]->clear();
            track_extrap_branches["track_y_layer_" + std::to_string(ilay)]->clear();
            track_extrap_branches["track_z_layer_" + std::to_string(ilay)]->clear();
        }
        int nTracks = Tracks_list.size();
        for (int track = 0; track < nTracks; track++)
        {
            TrckPDGID.push_back(Tracks_list.at(track).pdgcode);
            TrckPosInRealList.push_back(Tracks_list.at(track).nFinal_State_Particles);
            PerigeeA0.push_back(Tracks_list.at(track).a0);
            PerigeeZ0.push_back(Tracks_list.at(track).z0);
            PerigeeTheta.push_back(Tracks_list.at(track).theta);
            PerigeePhi.push_back(Tracks_list.at(track).phiHelix);
            PerigeeQ_P.push_back(Tracks_list.at(track).q_p);
            track_pflow_object_idx.push_back(-1);
            NotRegistrated.push_back(1*!Tracks_list.at(track).Is_reach_calorimeter+2*!Tracks_list.at(track).Is_inside_R+4*!Tracks_list.at(track).Is_track_reconstracted);
	    track_LHED.push_back( Tracks_list.at(track).GetLHED() );
            for (int ilay = 0; ilay < nLayers; ilay++)
            {
                track_extrap_branches["track_eta_layer_" + std::to_string(ilay)]->push_back(Tracks_list.at(track).eta.at(ilay));
                track_extrap_branches["track_phi_layer_" + std::to_string(ilay)]->push_back(Tracks_list.at(track).phi.at(ilay));
                track_extrap_branches["track_x_layer_" + std::to_string(ilay)]->push_back(Tracks_list.at(track).x_mid_layer.at(ilay));
                track_extrap_branches["track_y_layer_" + std::to_string(ilay)]->push_back(Tracks_list.at(track).y_mid_layer.at(ilay));
                track_extrap_branches["track_z_layer_" + std::to_string(ilay)]->push_back(Tracks_list.at(track).z_mid_layer.at(ilay));
            }
        }
        for (int ipflow = 0; ipflow < size_pflow; ipflow++)
        {
            if ((pflow_pbj.pflow_list.at(ipflow).isTrack))
            {
                track_pflow_object_idx.at(pflow_pbj.pflow_list.at(ipflow).label) = ipflow;
            }
        }
    }
}

void Tracks_data::set_tree_branches(TTree *outTree, int NLayers, std::string Type_of_running)
{
    nLayers = NLayers;
    if (Type_of_running != "PFlow_debug_E_p_template")
    {
	Config_reader_var &config_var = Config_reader_var::GetInstance();
        if ( config_var.doPFlow )
	    outTree->Branch("track_pflow_object_idx", "vector<int>", &track_pflow_object_idx);
        outTree->Branch("track_pdgid"           , "vector<int>", &TrckPDGID        );
        outTree->Branch("track_parent_idx"      , "vector<int>", &TrckPosInRealList   );
        outTree->Branch("track_d0"              , "vector<float>", &PerigeeA0   );
        outTree->Branch("track_z0"              , "vector<float>", &PerigeeZ0   );
        outTree->Branch("track_theta"           , "vector<float>", &PerigeeTheta);
        outTree->Branch("track_phi"             , "vector<float>", &PerigeePhi  );
        outTree->Branch("track_qoverp"          , "vector<float>", &PerigeeQ_P  );
        outTree->Branch("track_not_reg"         , "vector<int>", &NotRegistrated  );
        outTree->Branch("track_lhed"            , "vector<int>", &track_LHED );
    }
    for (int i = 0; i < nLayers; i++)
    {
        std::vector<float>* temp_eta_vec = new std::vector<float>();
        std::vector<float>* temp_phi_vec = new std::vector<float>();
        std::vector<float>* temp_x_vec = new std::vector<float>();
        std::vector<float>* temp_y_vec = new std::vector<float>();
        std::vector<float>* temp_z_vec = new std::vector<float>();
        track_extrap_branches["track_eta_layer_" + std::to_string(i)] = temp_eta_vec;
        if (Type_of_running != "PFlow_debug_E_p_template")
        {
            track_extrap_branches["track_phi_layer_" + std::to_string(i)] = temp_phi_vec;
            track_extrap_branches["track_x_layer_" + std::to_string(i)] = temp_x_vec;
            track_extrap_branches["track_y_layer_" + std::to_string(i)] = temp_y_vec;
            track_extrap_branches["track_z_layer_" + std::to_string(i)] = temp_z_vec;
            outTree->Branch(TString("track_eta_layer_" + std::to_string(i)), "vector<float>", temp_eta_vec);
            outTree->Branch(TString("track_phi_layer_" + std::to_string(i)), "vector<float>", temp_phi_vec);
            outTree->Branch(TString("track_x_layer_" + std::to_string(i)), "vector<float>", temp_x_vec);
            outTree->Branch(TString("track_y_layer_" + std::to_string(i)), "vector<float>", temp_y_vec);
            outTree->Branch(TString("track_z_layer_" + std::to_string(i)), "vector<float>", temp_z_vec);
        }

    }

}

void Tracks_data::Clear() 
{
    track_pflow_object_idx.clear();
    Tracks_list.clear();
    PerigeeA0.clear();
    PerigeeZ0.clear();
    PerigeeTheta.clear();
    PerigeePhi.clear();
    PerigeeQ_P.clear();
    NotRegistrated.clear();
    TrckPDGID.clear();
    TrckPosInRealList.clear();
    track_LHED.clear();

    //  track_extrap_branches.clear();
}


