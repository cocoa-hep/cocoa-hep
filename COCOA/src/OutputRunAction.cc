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
//
/// \file OutputRunAction.cc
/// \brief Implementation of the OutputRunAction class

#include "OutputRunAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4ProcessManager.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"

using namespace std;


std::string ConvertToStringForVector(std::string Name, int Int_Sufix)
{
	std::string Sufix = std::to_string(Int_Sufix);
	std::string cell_TopoEnergyStr = Name + "_" + Sufix;

	return cell_TopoEnergyStr;
}

OutputRunAction::OutputRunAction(std::string outputfilename, bool save_truth_graph) : G4UserRunAction()
{
	// set printing event number per each event
	G4RunManager::GetRunManager()->SetPrintProgress(1);
	m_outputfilename = outputfilename;
	m_save_truth_graph = save_truth_graph;

	// Create analysis manager
	// The choice of analysis technology is done via selectin of a namespace
	// in B4Analysis.hh
	auto analysisManager = G4AnalysisManager::Instance();

	// Create directories
	//analysisManager->SetHistoDirectoryName("histograms");
	//analysisManager->SetNtupleDirectoryName("ntuple");
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetNtupleMerging(true);
}

OutputRunAction::~OutputRunAction()
{
	;
}

// G4Run *OutputRunAction::GenerateRun()
// {;
// }

void OutputRunAction::BeginOfRunAction(const G4Run *run)
{
	(void)run;

	ChangeLooperParameters(G4Electron::Definition());
	// Get analysis manager
	// auto analysisManager = G4AnalysisManager::Instance();
	CSVReader &InitCSV = InitCSV.GetInstance();
	Config_reader_var &config_var = Config_reader_var::GetInstance();
	
	Tracks_data &track_list_low = Tracks_data::GetLow();
	Cells_data &cells_low = Cells_data::GetLow();
	Topo_clusts_data &topo_clusts = Topo_clusts_data::GetInstance();
	
	TruthRecordGraph &truth_record_graph = TruthRecordGraph::GetInstance();
	Full_trajectory_info_data &final_state_particle = Full_trajectory_info_data::GetInstance();
	Superclustering_data &superclustering_data = Superclustering_data::GetInstance();
	Graph_construction_data &graph_obj = Graph_construction_data::GetLow();
	Jet_Builder_data &pflow_jets_obj = Jet_Builder_data::Get_instance_pflow();
	Jet_Builder_data &topo_jets_obj  = Jet_Builder_data::Get_instance_topo();
	Jet_Builder_data &true_jets_obj  = Jet_Builder_data::Get_instance_true();

	// Open an output file
	//
	G4String fileName = "PFlowNtuple";
	//analysisManager->OpenFile(fileName);
	outf = new TFile(TString(m_outputfilename), "RECREATE");
	TString outTree_name =  config_var.Use_high_granularity ? "Low_Tree" : "Out_Tree";
	outTree_low = new TTree(outTree_name, outTree_name);


	if (config_var.Type_of_running == "Standard")
	{
		Particle_flow_data &pflow_obj = Particle_flow_data::GetInstance();
		cells_low.set_tree_branches(outTree_low);
		track_list_low.set_tree_branches(outTree_low, config_var.low_resolution.kNLayers);
		truth_record_graph.set_tree_branches(outTree_low);
		final_state_particle.set_tree_branches(outTree_low);
		topo_clusts.set_tree_branches(outTree_low);
		if ( config_var.doPFlow ) {
		    pflow_obj.set_tree_branches(outTree_low);
		    pflow_jets_obj.set_tree_branches(outTree_low);
		}
		if ( config_var.doSuperclustering) {
		    superclustering_data.set_tree_branches(outTree_low);
		}
		graph_obj.set_tree_branches(outTree_low);
		true_jets_obj.set_tree_branches(outTree_low);
		topo_jets_obj.set_tree_branches(outTree_low);
		if (config_var.Use_high_granularity)
		{
			outTree_high = new TTree("High_Tree", "High_Tree");
			Cells_data &cells_high = Cells_data::GetHigh();
			Tracks_data &track_list_high = Tracks_data::GetHigh();
			Graph_construction_data &graph_obj_high = Graph_construction_data::GetHigh();

			cells_high.set_tree_branches(outTree_high);
			track_list_high.set_tree_branches(outTree_high, config_var.high_resolution.kNLayers);
			graph_obj_high.set_tree_branches(outTree_high);
		}
	}
	else if (config_var.Type_of_running.find("Detector_") != string::npos)
	{
		Detector_analysis_var &det_ana_obj = Detector_analysis_var::GetInstance();
		det_ana_obj.set_tree_branches(outTree_low, config_var.Type_of_running);
	}
	else if (config_var.Type_of_running == "TopoClustering_debug")
	{
		Topo_clusts_data &topo_neutral = Topo_clusts_data::GetInstance_neutral();
		Topo_clusts_data &topo_charge = Topo_clusts_data::GetInstance_charge();
		Topo_clusts_data &topo_noise = Topo_clusts_data::GetInstance_noise();

		topo_neutral.set_tree_branches(outTree_low);
		topo_charge.set_tree_branches(outTree_low);
		topo_noise.set_tree_branches(outTree_low);
	}
	else if (config_var.Type_of_running.find("PFlow_debug") != string::npos)
	{
		Debug_Particle_Flow_data &pion_info = Debug_Particle_Flow_data::GetInstance();
		pion_info.set_tree_branches(outTree_low, config_var.Type_of_running);
		if (config_var.Type_of_running == "PFlow_debug_E_p_template")
			track_list_low.set_tree_branches(outTree_low, config_var.low_resolution.kNLayers);//Todo add option to save less for template 
	}

}

void OutputRunAction::EndOfRunAction(const G4Run * /*aRun*/)
{
	outf->cd();
	outTree_low->Write();
	Config_reader_var &config_var = Config_reader_var::GetInstance();
	if (config_var.Use_high_granularity) outTree_high->Write();
	outf->Close();
}

std::pair<G4Transportation *, G4CoupledTransportation *>
OutputRunAction::findTransportation(const G4ParticleDefinition *particleDef, bool reportError)
{
	const auto *partPM = particleDef->GetProcessManager();

	G4VProcess *partTransport = partPM->GetProcess("Transportation");
	auto transport = dynamic_cast<G4Transportation *>(partTransport);

	partTransport = partPM->GetProcess("CoupledTransportation");
	auto coupledTransport =
		dynamic_cast<G4CoupledTransportation *>(partTransport);

	if (reportError && !transport && !coupledTransport)
	{
		G4cerr << "Unable to find Transportation process for particle type "
			   << particleDef->GetParticleName()
			   << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
			   << G4endl;
	}

	return std::make_pair(transport, coupledTransport);
	// <G4Transportation*, G4CoupledTransportation*>
}

void OutputRunAction::ChangeLooperParameters(const G4ParticleDefinition *particleDef)
{
	if (particleDef == nullptr)
		particleDef = G4Electron::Definition();
	auto transportPair = findTransportation(particleDef);
	auto transport = transportPair.first;
	auto coupledTransport = transportPair.second;

	if (transport != nullptr)
	{
		// Change the values of the looping particle parameters of Transportation
		if (theWarningEnergy >= 0.0)
			transport->SetThresholdWarningEnergy(theWarningEnergy);
		if (theImportantEnergy >= 0.0)
			transport->SetThresholdImportantEnergy(theImportantEnergy);
		if (theNumberOfTrials > 0)
			transport->SetThresholdTrials(theNumberOfTrials);
	}
	else if (coupledTransport != nullptr)
	{
		// Change the values for Coupled Transport
		if (theWarningEnergy >= 0.0)
			coupledTransport->SetThresholdWarningEnergy(theWarningEnergy);
		if (theImportantEnergy >= 0.0)
			coupledTransport->SetThresholdImportantEnergy(theImportantEnergy);
		if (theNumberOfTrials > 0)
			coupledTransport->SetThresholdTrials(theNumberOfTrials);
	}
}
