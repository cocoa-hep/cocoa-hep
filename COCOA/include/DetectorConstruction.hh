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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef H02DetectorConstruction_h
#define H02DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Types.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "Config_reader_var.hh"
#include "CalorimeterConstruction.hh"
#include "InnerConstruction.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

long double EtaToTheta(long double eta);
long double ThetaToEta(long double theta);
long double GetPhi(long double px, long double py);

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
	DetectorConstruction(Geometry_definition Geometry);
	~DetectorConstruction();
	virtual G4VPhysicalVolume *Construct();

	Config_reader_var &config_json_var = Config_reader_var::GetInstance();
	G4Material *defaultMaterial;
	G4Material *lead;
	G4Material *iron;
	G4Material *plastic;
	G4Material *elSi;
	void Build_Iron_Gap_legacy(G4LogicalVolume *expHallLV);
	long double theta_min;

private:
	// G4VPhysicalVolume*   fGapPV;      // the gap physical volume
	Geometry_definition geometry;
	G4bool fCheckOverlaps; // option to activate checking of volumes overlaps
};

#endif
