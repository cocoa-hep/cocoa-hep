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
/// \file eventgenerator/HepMC/HepMCEx02/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
#ifndef H02_EVENT_ACTION_H
#define H02_EVENT_ACTION_H
#include "G4UserEventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Cells_data.hh"
#include "Tracking_func.hh"
#include "Tracks_data.hh"
#include "Track_var.hh"
#include "Topo_clust_func.hh"
#include "Topo_clust_var.hh"
#include "Topo_clusts_data.hh"
#include "Full_trajectory_info_data.hh"
#include "TruthRecordGraph.hh"
#include "Particle_flow_func.hh"
#include "Particle_flow_data.hh"
#include "Detector_analysis_var.hh"

class EventAction : public G4UserEventAction
{
public:
	EventAction();
	~EventAction();

	virtual void BeginOfEventAction(const G4Event *anEvent); //
	virtual void EndOfEventAction(const G4Event *anEvent);

private:
	Tracks_data &tracks_list_low = Tracks_data::GetLow();
	Cells_data &cells_data_high = Cells_data::GetHigh();
	Cells_data &cells_data_low = Cells_data::GetLow();
	TruthRecordGraph &truth_record_graph = TruthRecordGraph::GetInstance();
	Topo_clusts_data &topo_clusts = Topo_clusts_data::GetInstance();
	Particle_flow_data &pflow_obj = Particle_flow_data::GetInstance();
	Superclustering_data &superclustering_data = Superclustering_data::GetInstance();
	Full_trajectory_info_data &trajectories = Full_trajectory_info_data::GetInstance();
	Graph_construction_data &graph_obj = Graph_construction_data::GetLow();
	Graph_construction_data &graph_obj_high = Graph_construction_data::GetHigh();
	Jet_Builder_data &pflow_jets_obj = Jet_Builder_data::Get_instance_pflow();
	Jet_Builder_data &true_jets_obj = Jet_Builder_data::Get_instance_true();
	Jet_Builder_data &topo_jets_obj = Jet_Builder_data::Get_instance_topo();
	Debug_Particle_Flow_data &pion_info = Debug_Particle_Flow_data::GetInstance();
	Detector_analysis_var &det_ana = Detector_analysis_var::GetInstance();
};

#endif
