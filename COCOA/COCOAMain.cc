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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string> 
#include "OutputRunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#include "TrackingAction.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4Types.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Config_reader_func.hh"
#include "Config_reader_var.hh"

#include "RunTest.hh"

static void show_usage(std::string name)
{
	std::cerr << "Usage: \n" << name << " <option(s)> "
			  << "Options:\n"
			  << "\t--config (-c) <str>\t path to json configuration file\n"
			  << "\t--macro (-m) <str>\t path to Geant4, Pythia8, or HepMC macro file for event generation (can be set in json configuration file)\n"
			  << "\t--output (-o) <str>\t path (incl. name) of output ROOT file to be written (can be set in json configuration file)\n"
			  << "\t--input (-i) <str>\t path to HepMC (.hmc) input file (overrides the default path set in the HepMC macro file)\n"
			  << "\t--seed (-s) <int>\t set random seed\n"
			  << "\t--nevents (-n) <int>\t number of events to generate (default is taken from macro).\n"
			  << "\t--help (-h)\t show this message\n"
			  << "no <option(s)> will call UI interactive command submission\n" 
			  << std::endl;
}


int main(int argc, char **argv)
{
	std::string path_to_config = "./config/config_lowres.json";
	time_t systime = time(NULL);
	G4long seed = (long)systime;
	G4UIExecutive *ui = nullptr;
	std::string root_file_path = "";
	std::string macro_file_path = "";
	std::string input_file_path = "";
	int nEvents = -1;

	if (argc == 1)
	{
		ui = new G4UIExecutive(argc, argv);
	}
	else if (argc == 2)
	{
		std::string arg = argv[1];
		if (arg == "--help" || arg == "-h")
		{
			show_usage(argv[0]);
			return 0;
		}
		else
		{
			std::cerr << "Option "<< arg <<" does not exist!" << std::endl;
			show_usage(argv[0]);
		}
		return 0;
	}
	else
	{
		for (int i = 1; i < argc; ++i)
		{
			std::string arg = argv[i];
			if (arg == "--output" || arg == "-o")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					root_file_path = argv[i];
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--macro" || arg == "-m")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					macro_file_path = argv[i];
				}
				else
				{
					std::cerr << "--macro option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--seed" || arg == "-s")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					seed = stoi(argv[i]);
				}
				else
				{
					std::cerr << "--seed option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--nevents" || arg == "-n")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					nEvents = stoi(argv[i]);
				}
				else
				{
					std::cerr << "--nevents option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--config" || arg == "-c")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					path_to_config = argv[i];
				}
				else
				{
					std::cerr << "--config option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--input" || arg == "-i")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					input_file_path = argv[i];
				}
				else
				{
					std::cerr << "--input option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--help" || arg == "-h")
			{
				show_usage(argv[0]);
				return 0;
			}
			else
			{
				std::cerr << "Option "<< arg <<" does not exist!" << std::endl;
				show_usage(argv[0]);
				return 1;
			}
		}
	}
	Config_reader_var &config_var = Config_reader_var::GetInstance();
	Config_reader_func config_json_func(path_to_config, config_var);
	
	//* choose the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
	G4Random::setTheSeed(seed);

	CLHEP::HepRandom::setTheSeed(seed);
	// }
	Geometry_definition geometry = config_var.low_resolution;
	if (config_var.Use_high_granularity)
	{
		geometry = config_var.high_resolution;
	}

	G4RunManager *runManager = new G4RunManager;

	// User Initialization classes (mandatory)
	//
	G4VUserDetectorConstruction *detector = new DetectorConstruction(geometry);
	runManager->SetUserInitialization(detector);

	//
	G4VUserPhysicsList *physics = new FTFP_BERT;
	runManager->SetUserInitialization(physics);

	runManager->Initialize();

	G4VUserPrimaryGeneratorAction *gen_action = new PrimaryGeneratorAction;
	runManager->SetUserAction(gen_action);

	if ((root_file_path == ""))
	{
		if (config_var.Output_file_path!="")
			root_file_path = config_var.Output_file_path;
		else
		{
			G4cout<<"root_file_path is not given!"<<G4endl;
			return 1;
		}
	}
	
	OutputRunAction *outputrunaction = new OutputRunAction(root_file_path, config_var.Save_truth_particle_graph);
	runManager->SetUserAction(outputrunaction);
	//
	G4UserEventAction *event_action = new EventAction();
	runManager->SetUserAction(event_action);
	//

	G4UserTrackingAction *track_action = new TrackingAction;
	runManager->SetUserAction(track_action);

	G4UserSteppingAction *stepping_action = new SteppingAction(geometry);
	runManager->SetUserAction(stepping_action);

	G4VisManager *visManager = new G4VisExecutive;
	visManager->Initialize();

	//get the pointer to the User Interface manager
	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	if (!ui)
	{ // batch mode

		if ((macro_file_path == ""))
		{
			if (config_var.Macro_file_path!="")
				macro_file_path = config_var.Macro_file_path;
			else 
			{
				G4cout<<"macro_file_path is not given!"<<G4endl;
				return 1;
			}
		}
		
		visManager->SetVerboseLevel("quiet");
		UImanager->ApplyCommand("/generator/pythia8/setSeed " + std::to_string(seed));
		//UImanager->ApplyCommand(G4String("/control/execute ") + macro_file_path);

		ifstream filestream(macro_file_path);
		string line;

		//Parse the file line by line in order to overwrite with user input (e.g. nEvents)
		while ( std::getline(filestream, line) )
		{
			if      (line.empty()) continue;
			else if ( (nEvents > 0) && (line.find( "/run/beamOn" ) != std::string::npos) )
			{
				runManager->BeamOn(nEvents);
			}
			else if ( (input_file_path != "") && (line.find("/generator/hepmcAscii/open") != std::string::npos) )
			{
				if (input_file_path.find(".hmc") != std::string::npos)
				{
					UImanager->ApplyCommand(G4String("/generator/hepmcAscii/open " + input_file_path));
				}
				else
				{
					G4cout << "Input file is not a HepMC file (.hmc)!" << G4endl;
					return 1;
				}
			}
			else if (line.at(0)!='#')
			{
				G4cout << "Applying command: " <<  line << G4endl;
				UImanager->ApplyCommand(G4String(line));
			}
		}
	}
	else
	{ // interactive mode : define UI session
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		// if (ui->IsGUI())
		// {
		// 	UImanager->ApplyCommand("/control/execute gui.mac");
		// }

		ui->SessionStart();
		delete ui;
	}

	// Free the store: user actions, physics_list and detector_description are
	//                 owned and deleted by the run manager, so they should not
	//                 be deleted in the main() program !

	////////////////////////////
	// Test output if desired //
	////////////////////////////

	int exit_code = 0;

	test_t test_type = test_t::kNONE;
	if ( config_var.run_hadron_test )
	    test_type = kHADRON;
	if ( config_var.run_piZero_test )
	    test_type = kPI_ZERO;

	if ( test_type )
	    exit_code = RunTest( test_type );
	
	
	delete visManager;
	delete runManager;
	//  delete UImanager;

	return exit_code;
}
