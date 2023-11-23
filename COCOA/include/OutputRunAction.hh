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
/// \file OutputRunAction.hh
/// \brief Definition of the OutputRunAction class

#ifndef H02RunAction_h
#define H02RunAction_h 1

#include <string>

#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"

#include "CSVReader.hh"

#include "globals.hh"
// #include "DetectorGeometryDefinitions.hh"
//#include "g4root.hh"
// #include "DataStorage.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Version.hh"
#if G4VERSION_NUMBER >= 1100
#include "G4AnalysisManager.hh"
#endif
#include "G4UserRunAction.hh"
// #include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Track_var.hh"
#include "Tracks_data.hh"
#include "Cells_data.hh"
#include "Topo_clusts_data.hh"
#include "Particle_flow_data.hh"
#include "Superclustering_data.hh"
#include "TruthRecordGraph.hh"
#include "Full_trajectory_info_data.hh"
#include "Graph_construction_data.hh"
#include "Debug_Particle_Flow_data.hh"
#include "Jet_Builder_data.hh"
#include "Detector_analysis_var.hh"

using namespace std;

class G4Run;
class G4ParticleDefinition;
class G4Transportation;
class G4CoupledTransportation;

//#ifdef __MAKECINT__
//#pragma link C++ class vector<vector<vector<float> > >+;
//#endif

/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following
/// physics quantities:
/// The same values are also saved in the ntuple.
///
/// In EndOfRunAction(), the accumulated statistic and computed
/// dispersion is printed.
///


class OutputRunAction : public G4UserRunAction
{
public:
	OutputRunAction(std::string outputfilename, bool save_truth_graph = false);
	virtual ~OutputRunAction();

	// virtual G4Run *GenerateRun();

	virtual void BeginOfRunAction(const G4Run *);
	virtual void EndOfRunAction(const G4Run *);

	std::pair<G4Transportation *, G4CoupledTransportation *>
	findTransportation(const G4ParticleDefinition *particleDef,
					   bool reportError = true);
	void ChangeLooperParameters(const G4ParticleDefinition *particleDef);

	TFile *outf;
	TTree *outTree_high;
	TTree *outTree_low;

	std::string m_outputfilename;
	bool m_save_truth_graph;
	int typeofrun;


private:
	// Values for initialising 'loopers' parameters of Transport process
	G4int theNumberOfTrials = 10;			 // Default will not overwrite
	G4double theWarningEnergy = 300 * MeV;	 // Default values - non operational
	G4double theImportantEnergy = 300 * MeV; // Default - will not overwrite

	// int    theVerboseLevel = 0;
};

#endif
