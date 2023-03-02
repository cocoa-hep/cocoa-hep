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
/// \file eventgenerator/HepMC/HepMCEx02/include/SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
//
#ifndef H02_STEPPING_ACTION_H
#define H02_STEPPING_ACTION_H

#include "G4UserSteppingAction.hh"
#include "TRandom.h"
#include "Config_reader_var.hh"
#include "Cells_data.hh"
#include "Full_trajectory_info_data.hh"
#include "Detector_analysis_var.hh"


class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(Geometry_definition geometry);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* astep);

        int *CellIndex(const char* cellName, double RhoPos, double EtaPos, double PhiPos) ;
    
	private:
		Config_reader_var& config_json_var = Config_reader_var::GetInstance();
		std::vector< std::vector <long double> > High_cone_min_length_ECAL;
		std::vector< std::vector <long double> > High_cone_max_length_ECAL;
		std::vector< std::vector <long double> > High_cone_min_length_HCAL;
		std::vector< std::vector <long double> > High_cone_max_length_HCAL;
		std::vector <long double> cone_min_length_flatten;
		std::vector <long double> cone_max_length_flatten;
		long double theta_min;
		Geometry_definition geometry;
		char* Name_creation(char *name, int low_layer, int high_layer);
                    
};

#endif

