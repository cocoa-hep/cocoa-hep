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

// #include "OutputRunAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
// #include "DataStorage.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "DetectorGeometryDefinitions.hh"
//#include "FullTrajectoryInfo.hh"
#include <math.h>
#include <string>
#include <algorithm>

char* SteppingAction::Name_creation(char *name, int low_layer, int high_layer)
{
	name[4] = (low_layer+49);
	name[6] = (high_layer+49);
	return name;
}
SteppingAction::SteppingAction(Geometry_definition Geometry) : G4UserSteppingAction()
{

	geometry = Geometry;
	theta_min = 2 * atan(exp(-1 * config_json_var.max_eta_barrel));
	cone_min_length_flatten = geometry.layer_inn_radius_flatten;
	cone_max_length_flatten = geometry.layer_out_radius_flatten;
	int nLow_Layers = geometry.layer_out_radius_flatten.size();
	for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
	{
		cone_min_length_flatten.at(ilow_layer) = geometry.layer_inn_radius_flatten.at(ilow_layer) / tan(theta_min);
		cone_max_length_flatten.at(ilow_layer) = geometry.layer_out_radius_flatten.at(ilow_layer) / tan(theta_min);
	}
}

SteppingAction::~SteppingAction()
{
	;
}

int *SteppingAction::CellIndex(const char* cellName, double XPos, double YPos, double ZPos)
{
    //
    // Determine layer, eta and phi indices among the set of low resolution cells.
    //
    // For the layer index the cell name is used. It is assumed that this name starts with
    // ECAL<layer index + 1>_ or HCAL<layer index + 1>_ if the passed named actually belongs to
    // a calorimeter cell (more precisely, the 'CAL' subtring counts).
    //
	int R_Bin(-1), Eta_Bin(-1), Phi_Bin(-1);

	static int ZXYBin[3] = {-1};
	std::string name = cellName;
	if ( name.substr( 1, 3 ) != "CAL" ) {
	    ZXYBin[0] = -1;
	    ZXYBin[1] = -1;
	    ZXYBin[2] = -1;
	    return ZXYBin;
	}
	
	int         kFirstUnderscore      = name.find_first_of( "_" );
	int         nDigitsMainLayerIndex = kFirstUnderscore - 4;
	std::string nameSubLayerPart      = name.substr( kFirstUnderscore + 1, name.length() - ( kFirstUnderscore + 1 ) );
	int         nDigistSubLayerIndex  = nameSubLayerPart.find_first_of( "_" );

	int mainLayerIndex = atoi( name.substr( 4, nDigitsMainLayerIndex ).c_str() ) - 1;
	int subLayerIndex  = atoi( name.substr( 4 + nDigitsMainLayerIndex + 1, nDigistSubLayerIndex ).c_str() ) - 1;
	R_Bin              = subLayerIndex;
	//
	// Now count all previous sublayers and add the result to R_Bin
	//
	bool isECAL = true;
	if ( name.substr( 0, 1 ) == "H" )
	    isECAL = false;

	
	for( size_t iMainLayer = 0; iMainLayer < geometry.layer_inn_radius_ECAL.size(); ++iMainLayer ) {
	    if ( isECAL && iMainLayer == mainLayerIndex )
		break;
	    R_Bin += geometry.layer_inn_radius_ECAL[iMainLayer].size();
	}
	if ( !isECAL ) {
	    for( size_t iMainLayer = 0; iMainLayer < geometry.layer_inn_radius_HCAL.size(); ++iMainLayer ) {
		if ( iMainLayer == mainLayerIndex )
		    break;
		R_Bin += geometry.layer_inn_radius_HCAL[iMainLayer].size();
	    }
	}
	
	double r_sqr = pow(XPos, 2) + pow(YPos, 2);

	double PhiPos = atan2(YPos, XPos);
	if (PhiPos < 0)
		PhiPos += config_json_var.max_phi;
	// Phi_Bin = (int) floor(PhiPos/divided_tube_dPhi);

	double EtaPos = -1 * log(tan(0.5 * acos(ZPos / pow(r_sqr + pow(ZPos, 2), 0.5))));
	// Eta_Bin = (int) floor((config_json_var.max_eta_endcap+EtaPos)/d_eta_old);
	if (EtaPos >= -1 * config_json_var.max_eta_barrel && EtaPos <= config_json_var.max_eta_barrel)
	    //
	    // Phi definition used here does not perfectly match the one used in the cell definition. A rotation in phi is considered Ok though.
	    //
	{
		int nLow_Layers = geometry.number_of_pixels_flatten.size();
		for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
		{
		    if ( r_sqr >= pow(geometry.layer_inn_radius_flatten.at(ilow_layer), 2) &&
			 r_sqr < pow(geometry.layer_out_radius_flatten.at(ilow_layer), 2) ) {
			
				Phi_Bin = (int)floor(PhiPos / geometry.layer_dphi_flatten.at(ilow_layer));
				Eta_Bin = (int)floor((config_json_var.max_eta_endcap + EtaPos) / (geometry.layer_deta_flatten.at(ilow_layer)));
				break;

		    }
		}
	}
	else if (EtaPos >= -1 * config_json_var.max_eta_endcap && EtaPos < -1 * config_json_var.max_eta_barrel)
	{
		int nLow_Layers = geometry.number_of_pixels_flatten.size();
		for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
		{
			if (ZPos >= -1 * cone_max_length_flatten.at(ilow_layer) &&
				ZPos < -1 * cone_min_length_flatten.at(ilow_layer))
			{
				Phi_Bin = (int)floor(PhiPos / geometry.layer_dphi_flatten.at(ilow_layer));
				Eta_Bin = (int)floor((config_json_var.max_eta_endcap + EtaPos) / (geometry.layer_deta_flatten.at(ilow_layer)));
				break;
			}
		}
	}
	else if (EtaPos > config_json_var.max_eta_barrel && EtaPos <= config_json_var.max_eta_endcap)
	{
		int nLow_Layers = geometry.number_of_pixels_flatten.size();
		for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
		{
			if (ZPos >= cone_min_length_flatten.at(ilow_layer) &&
				ZPos < cone_max_length_flatten.at(ilow_layer))
			{
				Phi_Bin = (int)floor(PhiPos / geometry.layer_dphi_flatten.at(ilow_layer));
				Eta_Bin = (int)floor((config_json_var.max_eta_endcap + EtaPos) / (geometry.layer_deta_flatten.at(ilow_layer)));
				break;
			}
		}
	}

	ZXYBin[0] = R_Bin;
	ZXYBin[1] = Eta_Bin;
	ZXYBin[2] = Phi_Bin;
	return ZXYBin;
}


void SteppingAction::UserSteppingAction(const G4Step *astep)
{
	G4Track *aTrack = astep->GetTrack();
	
	// prevent infinite loops
	if ( aTrack->GetCurrentStepNumber() > 1e5 )
	    aTrack->SetTrackStatus( fStopAndKill );
	
	auto edep = astep->GetTotalEnergyDeposit();
	Full_trajectory_info_data &trajectories = Full_trajectory_info_data::GetInstance();
	
	G4ThreeVector PreStepPoint = astep->GetPreStepPoint()->GetPosition();
	G4ThreeVector PostStepPoint = astep->GetPostStepPoint()->GetPosition();
	G4ThreeVector momentum      = aTrack->GetDynamicParticle()->GetMomentum();
	G4TouchableHandle touch1 = astep->GetPreStepPoint()->GetTouchableHandle();
	G4LogicalVolume *lvol    = touch1->GetVolume()->GetLogicalVolume();
	G4int ParentID           = aTrack->GetParentID();
	// G4int originID           = -1;
	// if ( aTrack->GetDynamicParticle()->GetPrimaryParticle() )
	//     originID = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTrackID();
	G4int trackID            = aTrack->GetTrackID();
	G4int trackPdgId         = aTrack->GetDefinition()->GetPDGEncoding();
	std::string volume_name  = touch1->GetVolume()->GetName();
	G4double eKin            = aTrack->GetKineticEnergy();
	
	if ( volume_name.substr(1, 3) == "CAL" ) {
	    
	    for ( size_t iPrimaryParticle = 0; iPrimaryParticle < trajectories.fAllTrajectoryInfo.size(); ++iPrimaryParticle ) {
		
		FullTrajectoryInfo &trajectory = trajectories.fAllTrajectoryInfo[iPrimaryParticle];
		
		if ( trackPdgId == trajectory.fPDGCode &&
		     std::find( trajectory.vTrackID.begin(), trajectory.vTrackID.end(), trackID ) != trajectory.vTrackID.end() ) {
		    
		    if ( eKin > trajectory.caloExtrapolMaxEkin ) {
			
			trajectory.caloExtrapolMaxEkin = eKin;
			trajectory.caloExtrapolEta     = ThetaToEta( acos( PreStepPoint.z() / sqrt( sqr( PreStepPoint.x() ) + sqr( PreStepPoint.y() ) + sqr( PreStepPoint.z() ) ) ) );
			trajectory.caloExtrapolPhi     = GetPhi( PreStepPoint.x(),
								 PreStepPoint.y() );
			break;
			
		    }
		    
		}
		
	    }
	    
	} else if ( std::strstr( volume_name.c_str(), "outermostInner" ) ) {

	    for ( size_t iPrimaryParticle = 0; iPrimaryParticle < trajectories.fAllTrajectoryInfo.size(); ++iPrimaryParticle ) {
		
		FullTrajectoryInfo &trajectory = trajectories.fAllTrajectoryInfo[iPrimaryParticle];
		
		if ( trackPdgId == trajectory.fPDGCode &&
		     std::find( trajectory.vTrackID.begin(), trajectory.vTrackID.end(), trackID ) != trajectory.vTrackID.end() ) {
		    
		    if ( eKin > trajectory.idExtrapolMaxEkin ) {
			
			trajectory.idExtrapolMaxEkin = eKin;
			trajectory.idExtrapolEta     = ThetaToEta( acos( PreStepPoint.z() / sqrt( sqr( PreStepPoint.x() ) + sqr( PreStepPoint.y() ) + sqr( PreStepPoint.z() ) ) ) );
			trajectory.idExtrapolPhi     = GetPhi( PreStepPoint.x(),
							       PreStepPoint.y() );
			break;
			
		    }
		    
		}
		
	    }
	    
	}
	
	//*Geantino

	G4ParticleDefinition* particle0 = aTrack->GetDefinition();
	if (particle0 == G4Geantino::Geantino()&&lvol->GetName()!="EXP_HALL_LV")//G4ChargedGeantino::ChargedGeantino()
	{
		Detector_analysis_var &det_ana_obj = Detector_analysis_var::GetInstance();
		G4double int_l  = lvol->GetMaterial()->GetNuclearInterLength();//GetNuclearInterLength()/cm;
		G4double rad_l  = lvol->GetMaterial()->GetRadlen();//GetNuclearInterLength()/cm;
		G4double step_l = astep->GetStepLength();
		int lay_count = 0;
		bool if_track = true;
		int num_ecal_layers = geometry.number_of_pixels_ECAL.size();
		for (int iecal_low = 0; iecal_low < num_ecal_layers; iecal_low++)
		{
			int num_sub_ecal_layrs = geometry.number_of_pixels_ECAL.at(iecal_low).size();
			for (int iecal_high = 0; iecal_high < num_sub_ecal_layrs; iecal_high++)
			{
				if (lvol->GetName()==Name_creation(strdup("ECALN_N_forward_LV"),iecal_low,iecal_high)||
 				lvol->GetName()==Name_creation(strdup("ECALN_N_back_LV"),iecal_low,iecal_high)||
				lvol->GetName()==Name_creation(strdup("ECALN_N_Endcap_forward_LV"),iecal_low,iecal_high)||
				lvol->GetName()==Name_creation(strdup("ECALN_N_Endcap_back_LV"),iecal_low,iecal_high))
				{
					det_ana_obj.add_lengths(step_l/rad_l, step_l/int_l,lay_count + 1);
					if_track = false;
				}
				lay_count++;
			}
		}
		lay_count++;
		int num_hcal_layers = geometry.number_of_pixels_HCAL.size();
		for (int ihcal_low = 0; ihcal_low < num_hcal_layers; ihcal_low++)
		{
			int num_sub_hcal_layers = geometry.number_of_pixels_HCAL.at(ihcal_low).size();
			for (int ihcal_high = 0; ihcal_high < num_sub_hcal_layers; ihcal_high++)
			{
				if (lvol->GetName()==Name_creation(strdup("HCALN_N_forward_LV"),ihcal_low,ihcal_high)||
				lvol->GetName()==Name_creation(strdup("HCALN_N_back_LV"),ihcal_low,ihcal_high)||
				lvol->GetName()==Name_creation(strdup("HCALN_N_Endcap_forward_LV"),ihcal_low,ihcal_high)||
				lvol->GetName()==Name_creation(strdup("HCALN_N_Endcap_back_LV"),ihcal_low,ihcal_high))
				{
					det_ana_obj.add_lengths(step_l/rad_l, step_l/int_l,lay_count + 1);
					if_track = false;
				}
				lay_count++;
			}
		}
		if (lvol->GetName().substr(3, 7) == "ironGap" )
		{
			lay_count = num_ecal_layers;
			det_ana_obj.add_lengths(step_l/rad_l, step_l/int_l,lay_count + 1);
		}
		else if (if_track)
		{
			det_ana_obj.add_lengths(step_l/rad_l, step_l/int_l,0);
		}
	}

	//* Geantino end
	if (aTrack->GetTrackStatus() != fAlive && edep == 0.0)
		return;

	if ((int)trajectories.fAllTrajectoryInfo.size() == 0)
		return;
	// ----- assigning the trajectory to step ---- //
	bool foundTraj(false);
	int mTraj(-1);
	// int mParent(-1);
	// ----- search over all the trajectories ---------- //
	for (int iTraj = (int)trajectories.fAllTrajectoryInfo.size() - 1; iTraj >= 0; iTraj--)
	{
		// ----- search for all the existing trackID in a given trajectory ---- //
		for (int iParent = (int)trajectories.fAllTrajectoryInfo.at(iTraj).vTrackID.size() - 1; iParent >= 0; iParent--)
		{

			if (ParentID == trajectories.fAllTrajectoryInfo.at(iTraj).vTrackID.at(iParent))
			{
				foundTraj = true;
				mTraj = iTraj;
				break;
			} // if( trajectories.fAllTrajectoryInfo.at(iTraj).vParentID[iParent] == ParentID )

		} // for(int iParent = 0; iParent < (int)trajectories.fAllTrajectoryInfo.at(iTraj).vParentID.size(); iParent++  )

		if (foundTraj)
			break;
	} // for(int iTraj = 0; iTraj < (int)trajectories.fAllTrajectoryInfo.size(); iTraj++)

	// ==================================================
	// Emulation of Energy Sampling : Use the user-defined
	// sampling fraction to randomly decide whether an energy
	// deposition is taken into account.
	// ==================================================
	float samplingFraction = ( volume_name.substr( 0, 1 ) == "E" ? config_json_var.samplingFraction_ECAL : config_json_var.samplingFraction_HCAL );
	if ( gRandom->Uniform() > samplingFraction )
	    edep = 0.0;
	
	if (foundTraj && edep > 0.)
	{

	        std::string volume_name = touch1->GetVolume()->GetName();
	        int *Bin                = CellIndex( volume_name.c_str(),
						     PreStepPoint.x(),
						     PreStepPoint.y(),
						     PreStepPoint.z() );

		G4double Charge = trajectories.fAllTrajectoryInfo.at(mTraj).fPDGCharge;
		G4double Etot(0), Ech(0.), Enu(0.); //, EHadCh(0.), EHadNu(0.), EEM(0.);
		if (Charge == 0)					//|| (abs(trajectories.fAllTrajectoryInfo[mTraj].fPDGCode)==11)
		{
			Enu = edep;
		}
		else
		{
			Ech = edep;
		}

		Particle_dep_in_cell *conv_el = nullptr;
		for ( size_t iConvEl = 0; iConvEl < trajectories.fAllConvElectrons.size(); ++iConvEl ) {
		    if ( std::find( trajectories.fAllConvElectrons[iConvEl].vTrackID.begin(),
				    trajectories.fAllConvElectrons[iConvEl].vTrackID.end(),
				    trackID ) != trajectories.fAllConvElectrons[iConvEl].vTrackID.end() ) {
			conv_el = new Particle_dep_in_cell;
			conv_el->PDG_ID = trajectories.fAllConvElectrons[iConvEl].fPDGCode;
			conv_el->particle_pos_in_true_list = iConvEl;
			conv_el->Energy = edep / samplingFraction;
			break;
		    }
		}
		
		
		if ((*Bin >= 0 && *Bin < geometry.kNLayers))
		{
			if ((*(Bin + 1) >= 0 && *(Bin + 1) < geometry.number_of_pixels_flatten.at(*Bin)) && (*(Bin + 2) >= 0 && *(Bin + 2) < geometry.number_of_pixels_flatten.at(*Bin)))
			{
				Etot = Ech + Enu;
				// runData->AddTotalEnergy(Ech, Enu);
				Particle_dep_in_cell ptrc;

				ptrc.Energy = Etot / samplingFraction;
				ptrc.PDG_ID = trajectories.fAllTrajectoryInfo.at(mTraj).fPDGCode;

				ptrc.particle_pos_in_true_list = mTraj;
				
				if (config_json_var.Use_high_granularity)
				{
					Cells_data &cells_data = Cells_data::GetHigh();
					cells_data.add_cell_info(*Bin, *(Bin + 1), *(Bin + 2), Ech / samplingFraction, Enu / samplingFraction, ptrc, conv_el);
				}
				else
				{
					Cells_data &cells_data = Cells_data::GetLow();
					cells_data.add_cell_info(*Bin, *(Bin + 1), *(Bin + 2), Ech / samplingFraction, Enu / samplingFraction, ptrc, conv_el);
				}
			}
			// runData->AddCaloCell(*Bin, *(Bin+1), *(Bin+2), Ech, Enu, ptrc);
		} // if(*Bin >= 0 && *(Bin+1) >= 0 && *(Bin+2) >= 0 )
		delete conv_el;
	}	  // if(foundTraj && edep > 0.)
	G4ParticleDefinition *particleType = aTrack->GetDefinition();
	if ((particleType == G4MuonPlus::MuonPlusDefinition()) || (particleType == G4MuonMinus::MuonMinusDefinition()))
		return;

}

