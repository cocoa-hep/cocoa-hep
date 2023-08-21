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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
// #include "OutputRunAction.hh"
#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4ChordFinder.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVParameterised.hh"
#include "G4ReflectionFactory.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include <json/json.h>
#include "G4RunManager.hh"
#include "G4Cons.hh"

// G4ThreadLocal
// G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

long double ThetaToEta(long double theta) {
    return - log( tan( 0.5 * theta ) );
}

long double EtaToTheta(long double eta) {
    return 2.0 * atan( exp( - eta ) );
}

long double GetPhi(long double px, long double py) {
	long double phi = atan2( py, px );
	return phi;    
}

DetectorConstruction::DetectorConstruction(Geometry_definition Geometry): G4VUserDetectorConstruction()
{
	fCheckOverlaps = config_json_var.check_geometry_overlap;
	geometry = Geometry;
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
	// ==============================================================
	// Materials
	// ==============================================================

	G4NistManager *nistManager = G4NistManager::Instance();
	iron = nistManager->FindOrBuildMaterial("G4_Fe");
	elSi = nistManager->FindOrBuildMaterial("G4_Si");

	// Argon gas
	long double a, z, density;
	density = 1.782e-03 * g / cm3;

	defaultMaterial = new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = CLHEP::universe_mean_density,
												 kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

	// Print materials


	theta_min = 2 * atan(exp(-1 * config_json_var.max_eta_barrel));
	// long double half_lengthZ = r_inn*(1.0/(tan(theta_min)));
	long double length = geometry.layer_out_radius_flatten.back() / tan(theta_min);

	long double worldSizeXY = 4. * geometry.layer_out_radius_flatten.back(); //#1.2 * Total_Calo_Length;
	long double worldSizeZ = 4. * length;											 //1.5 * Total_Calo_Length;

	// ==============================================================
	// Experimental Hall (world)
	// ==============================================================
	G4Box *expHallSolid =
		new G4Box("EXP_HALL", worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2);

	G4LogicalVolume *expHallLV =
		new G4LogicalVolume(expHallSolid, defaultMaterial, "EXP_HALL_LV");

	// visualization attributes
	G4VisAttributes *expHallVisAtt =
		new G4VisAttributes(true, G4Colour(1., 1., 1.));
	//expHallVisAtt-> SetForceWireframe(TRUE);
	expHallLV->SetVisAttributes(expHallVisAtt);

	G4PVPlacement *expHall = new G4PVPlacement(
		0,				 // no rotation
		G4ThreeVector(), // at (0,0,0)
		expHallLV,		 // its logical volume
		"Exp_HALL",		 // its name
		0,				 // its mother  volume
		false,			 // no boolean operation
		0,				 // copy number
		fCheckOverlaps);

	CalorimeterConstruction Calorimeter(expHallLV, fCheckOverlaps, geometry);
	InnerConstruction       InnerDetector(expHallLV, defaultMaterial, iron, elSi, fCheckOverlaps);
	return expHall;
}

void DetectorConstruction::Build_Iron_Gap_legacy(G4LogicalVolume *expHallLV)
{
	long double l_gap = (geometry.layer_out_radius_ECAL.back().back()) / tan(theta_min);
	//* Make gap between ECAL3 and HCAL1

	G4Tubs *Gap = new G4Tubs("Inner_part", (geometry.layer_out_radius_ECAL.back().back()),
							 geometry.layer_inn_radius_HCAL.front().front(), l_gap, 0, config_json_var.max_phi);
	G4LogicalVolume *Gap_LV = new G4LogicalVolume(Gap, iron, "Gap_LV");

	new G4PVPlacement(
		0, // no rotation
		G4ThreeVector(0., 0., 0.),
		Gap_LV,	   // its logical volume
		"Gap_PL",  // its name
		expHallLV, // its mother  volume
		false,	   // no boolean operation
		0,		   // copy number
		fCheckOverlaps);

	double length_cone_min = ((geometry.layer_out_radius_ECAL.back().back())) / tan(theta_min);
	double length_cone_max = length_cone_min + (geometry.layer_inn_radius_HCAL.front().front() - geometry.layer_out_radius_ECAL.back().back()) / tan(theta_min);
	for (int direction = 1; direction > -2; direction = direction - 2)
	{
		if (direction == -1)
		{
			double buf = length_cone_min;
			length_cone_min = length_cone_max;
			length_cone_max = buf;
		}
		long double d_theta_next = 2 * atan(exp(-1 * (config_json_var.max_eta_endcap)));
		G4Tubs *Cone_Gap = new G4Tubs("Cone_Gap", length_cone_max*tan(d_theta_next), geometry.layer_inn_radius_HCAL.front().front(),
									  fabs((geometry.layer_inn_radius_HCAL.front().front() - (geometry.layer_out_radius_ECAL.back().back() + 1e-3)) / (2 * tan(theta_min))), 0, 2 * M_PI);
									  
		// G4Cons *Cone_Gap = new G4Cons("Cone_Gap", (length_cone_max)*tan(d_theta_next), geometry.layer_inn_radius_HCAL.front().front(),
		// 							  (length_cone_min)*tan(d_theta_next), geometry.layer_out_radius_ECAL.back().back(),
		// 							  fabs((geometry.layer_inn_radius_HCAL.front().front() - (geometry.layer_out_radius_ECAL.back().back() + 1e-3)) / (2 * tan(theta_min))), 0, 2 * M_PI);

		G4LogicalVolume *Iron_LV_posdir = new G4LogicalVolume(Cone_Gap, iron, "Cone_Gap_LV");

		new G4PVPlacement(
			0,																																																											//* rotation
			G4ThreeVector(0, 0, -1 * direction * ((geometry.layer_out_radius_ECAL.back().back()) + 0.5 * (geometry.layer_inn_radius_HCAL.front().front() - (geometry.layer_out_radius_ECAL.back().back()))) / tan(theta_min)), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
			Iron_LV_posdir,																																																								//* its logical volume
			"Iron_PL",																																																									//* its name
			expHallLV,																																																									//* its mother  volume
			false,																																																										//* no boolean operation
			0,																																																											//* copy number
			fCheckOverlaps);
	}
}
