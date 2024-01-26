//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file eventgenerator/HepMC/HepMCEx02/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
#include "OutputRunAction.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4Trajectory.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "HepMCG4Pythia8Interface.hh"
// #include "DataStorage.hh"
#include "Pythia8/Pythia.h"
#include "ReduceResolution.hh"
#include "Debug_Particle_Flow_func.hh"
#include "Superclustering.hh"
#include "GraphConstructor.hh"
#include "Graph_construction_data.hh"
#include "Jet_Builder_func.hh"
#include "Jet_Builder_data.hh"

#include <memory>

using namespace std;

EventAction::EventAction() : G4UserEventAction()
{;}

EventAction::~EventAction()
{
	;
}

void EventAction::BeginOfEventAction(const G4Event *anEvent) //const G4Event* anEvent
{
	tracks_list_low.Clear();
	cells_data_low.clear();
	cells_data_high.clear();
	topo_clusts.clear();
	//super_clusts.clear();
	pflow_obj.clear();
	superclustering_data.clear();
	trajectories.clear();
	graph_obj.clear();
	graph_obj_high.clear();
	true_jets_obj.clear();
	pflow_jets_obj.clear();
	topo_jets_obj.clear();
	pion_info.clear();
	det_ana.clear();
	// const G4Event* ev = anEvent;

#ifdef DEBUG_HEPMC
	// # This is not active
	// printout primary information

	G4int nVtx = anEvent->GetNumberOfPrimaryVertex();
	G4int i;
	for (i = 0; i < nVtx; i++)
	{
		const G4PrimaryVertex *primaryVertex = anEvent->GetPrimaryVertex(i);
		primaryVertex->Print();
	}
#endif
	(void)anEvent;
}

void EventAction::EndOfEventAction(const G4Event *evt)
{
	(void)evt;

	Config_reader_var &config_var = Config_reader_var::GetInstance();
	auto runAction = (OutputRunAction *)G4RunManager::GetRunManager()->GetUserRunAction();
	if (config_var.Type_of_running == "Detector_Response")
	{
		det_ana.transverse_energy_calculation(cells_data_low.fCell_array);
	}
	else if (config_var.Type_of_running == "TopoClustering_debug")
	{
		ReduceResolution noise_apply;
		noise_apply.apply_noise(cells_data_low.fCell_array);

		std::vector<std::vector<std::vector<Cell>>> neutral = cells_data_low.fCell_array;
		Topo_clust_func clustering_neutral(neutral, config_var.low_resolution, config_var.topological_clustering, "neutral");
		Topo_clusts_data &topo_neutral = Topo_clusts_data::GetInstance_neutral();
		topo_neutral.clear();
		clustering_neutral.topoclustering(topo_neutral.topo_clusts_list);
		topo_neutral.fill_topo_var();

		std::vector<std::vector<std::vector<Cell>>> charge = cells_data_low.fCell_array;
		Topo_clust_func clustering_charge(charge, config_var.low_resolution, config_var.topological_clustering, "charge");
		Topo_clusts_data &topo_charge = Topo_clusts_data::GetInstance_charge();
		topo_charge.clear();
		clustering_charge.topoclustering(topo_charge.topo_clusts_list);
		topo_charge.fill_topo_var();

		std::vector<std::vector<std::vector<Cell>>> noise = cells_data_low.fCell_array;
		Topo_clust_func clustering_noise(noise, config_var.low_resolution, config_var.topological_clustering, "noise");
		Topo_clusts_data &topo_noise = Topo_clusts_data::GetInstance_noise();
		topo_noise.clear();
		clustering_noise.topoclustering(topo_noise.topo_clusts_list);
		topo_noise.fill_topo_var();
	}
	else if (config_var.Type_of_running.find("PFlow_debug") != string::npos)
	{
		Tracking tracking_low(trajectories.fAllTrajectoryInfo, true, tracks_list_low.Tracks_list, config_var.low_resolution);
		ReduceResolution noise_apply;
		noise_apply.apply_noise(cells_data_low.fCell_array);
		Topo_clust_func clustering(cells_data_low.fCell_array, config_var.low_resolution, config_var.topological_clustering, "Standard");
		clustering.topoclustering(topo_clusts.topo_clusts_list);
		cells_data_low.fill_cells_in_topoclusters();
		Debug_Particle_Flow_func debug_pflow(tracks_list_low.Tracks_list, topo_clusts.topo_clusts_list,
										cells_data_low.Cells_in_topoclusters, pion_info, config_var.low_resolution,
										config_var.particle_flow, config_var.Type_of_running);
		pion_info.fill_var();
		tracks_list_low.Fill_perigee_var();
	}
	else if (config_var.Type_of_running == "Standard")
	{
		if (config_var.Use_high_granularity)
		{
			Tracking tracking_low(trajectories.fAllTrajectoryInfo, false, tracks_list_low.Tracks_list, config_var.low_resolution);

			ReduceResolution Down_ptr(cells_data_high.fCell_array, cells_data_low.fCell_array);
			Topo_clust_func clustering(cells_data_low.fCell_array, config_var.low_resolution, config_var.topological_clustering, "Standard");
			clustering.topoclustering(topo_clusts.topo_clusts_list);
			Down_ptr.link_apply(cells_data_high.fCell_array, cells_data_low.fCell_array);
			cells_data_low.fill_cells_in_topoclusters();
			cells_data_high.fill_cells_in_topoclusters();

			if ( config_var.doPFlow )
			    Particle_flow_func pflow(tracks_list_low.Tracks_list, topo_clusts.topo_clusts_list, cells_data_low.Cells_in_topoclusters, pflow_obj.pflow_list, config_var.low_resolution, config_var.particle_flow);

			tracks_list_low.Fill_perigee_var();
			cells_data_low.fill_cell_var();
			cells_data_high.fill_cell_var();
			if ( config_var.doPFlow )
			    pflow_obj.fill_cell_var();
			trajectories.fill_var();

			GraphConstructor graph_construct_high(cells_data_high.Cells_in_topoclusters, tracks_list_low.Tracks_list, trajectories.particle_to_track, graph_obj_high);

			std::vector<float> _particle_dep_energies;
			GraphConstructor graph_construct(cells_data_low.Cells_in_topoclusters, cells_data_high.Cells_in_topoclusters, tracks_list_low.Tracks_list, trajectories.particle_to_track, graph_obj, &_particle_dep_energies);

			trajectories.SetParticleDepEnergy( _particle_dep_energies );
			
			Jet_Builder_func jets_build;
			if ( config_var.doPFlow ) {
			    pflow_obj.make_pseudo_jet_particles();
			    jets_build.build_jets(pflow_obj.jets_objects, pflow_jets_obj, config_var.jet_parameter);
			    pflow_jets_obj.fill_cell_var();
			    jets_build.reset();
			}
			trajectories.make_pseudo_jet_particles();
			jets_build.build_jets(trajectories.jets_objects, true_jets_obj, config_var.jet_parameter);
			true_jets_obj.fill_cell_var(config_var.jet_parameter.radius);
			jets_build.reset();
			topo_clusts.make_pseudo_jet_particles();
			jets_build.build_jets(topo_clusts.jets_objects, topo_jets_obj, config_var.jet_parameter);
			topo_jets_obj.fill_cell_var();
			jets_build.reset();

			runAction->outTree_high->Fill();
		}
		else
		{
			if ( config_var.doSuperclustering )
			{
				//Append conversion tracks to primary tracks list so we can do 
				//superclustering of both electrons and photon conversions in one go
				for (int iconv=0; iconv < trajectories.fAllConvElectrons.size(); iconv++)
					trajectories.fAllTrajectoryInfo.push_back(trajectories.fAllConvElectrons.at(iconv));
			}

			Tracking tracking_low(trajectories.fAllTrajectoryInfo, false, tracks_list_low.Tracks_list, config_var.low_resolution);

			ReduceResolution noise_apply;
			noise_apply.apply_noise(cells_data_low.fCell_array);
			Topo_clust_func clustering(cells_data_low.fCell_array, config_var.low_resolution, config_var.topological_clustering, "Standard");
			clustering.topoclustering(topo_clusts.topo_clusts_list);
			cells_data_low.fill_cells_in_topoclusters();
			topo_clusts.fill_topo_var();

			if ( config_var.doPFlow )
			    Particle_flow_func pflow(tracks_list_low.Tracks_list, topo_clusts.topo_clusts_list, cells_data_low.Cells_in_topoclusters, pflow_obj.pflow_list, config_var.low_resolution, config_var.particle_flow);
			tracks_list_low.Fill_perigee_var();
			cells_data_low.fill_cell_var();
			if ( config_var.doPFlow )
			    pflow_obj.fill_cell_var();
			trajectories.fill_var();
			if ( config_var.doSuperclustering )
			{
				Superclustering superclusters(tracks_list_low.Tracks_list, topo_clusts.topo_clusts_list, cells_data_low.Cells_in_topoclusters, superclustering_data.super_list);
				superclustering_data.fill_supercluster_data();
			}

			std::vector<float> _particle_dep_energies;
			GraphConstructor graph_construct(cells_data_low.Cells_in_topoclusters, tracks_list_low.Tracks_list, trajectories.particle_to_track, graph_obj, &_particle_dep_energies);

			trajectories.SetParticleDepEnergy( _particle_dep_energies );

			Jet_Builder_func jets_build;
			if ( config_var.doPFlow ) {
			    pflow_obj.make_pseudo_jet_particles();
			    jets_build.build_jets(pflow_obj.jets_objects, pflow_jets_obj, config_var.jet_parameter);
			    pflow_jets_obj.fill_cell_var();
			    jets_build.reset();
			}
			trajectories.make_pseudo_jet_particles();
			jets_build.build_jets(trajectories.jets_objects, true_jets_obj, config_var.jet_parameter);
			true_jets_obj.fill_cell_var(config_var.jet_parameter.radius);
			jets_build.reset();
			topo_clusts.make_pseudo_jet_particles();
			jets_build.build_jets(topo_clusts.jets_objects, topo_jets_obj, config_var.jet_parameter);
			topo_jets_obj.fill_cell_var();
			jets_build.reset();
		}
	}

	runAction->outTree_low->Fill();
	truth_record_graph.clear();

}
