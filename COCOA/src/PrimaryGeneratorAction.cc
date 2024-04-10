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

#include "OutputRunAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "HepMCG4Pythia8Interface.hh"
#include "HepMCG4Reader.hh"
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4GeneralParticleSource.hh"

// #include "DataStorage.hh"

#include "G4RunManager.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction()
{
	// default generator is particle gun.
	fCurrentGenerator = fParticleGun = new G4ParticleGun();
	//fCurrentGenerator= fParticleGun= new G4GeneralParticleSource();
	fCurrentGeneratorName = "particleGun";
	// fHepmc = new HepMCG4Reader();
	//#ifdef G4LIB_USE_PYTHIA
	fPythiaGen = new HepMCG4Pythia8Interface();
	fHepMCGen  = new HepMCG4Reader();
	//#else
	//  fPythiaGen= 0;
	//#endif

	fParticleGun_ = new G4ParticleGun(1);

	fGentypeMap["particleGun"] = fParticleGun;
	fGentypeMap["hepmc"]  = fHepMCGen;
	fGentypeMap["pythia8"]     = fPythiaGen;

	fMessenger = new PrimaryGeneratorMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
	// auto runAction = (OutputRunAction *)G4RunManager::GetRunManager()->GetUserRunAction();
	if (fCurrentGenerator)
	{
		//* particleGun
		if (fCurrentGeneratorName == "particleGun")
		{
			Detector_analysis_var &det_ana_obj = Detector_analysis_var::GetInstance();
			G4double Phi0 = 0; //phi_fr + (phi_to - phi_fr) * G4UniformRand();
			G4double eta1 = det_ana_obj.get_eta();//runAction->eta_step; //a + (b - a) * G4UniformRand();

			double PtPi_P = 15 * GeV; //(10 + 5 * G4UniformRand())  * GeV;
			double pxPi_P_init = PtPi_P * cos(Phi0);
			double pyPi_P_init = PtPi_P * sin(Phi0);
			double pzPi_P_init = PtPi_P * sinh(eta1);

			fParticleGun_ = nullptr;
			fParticleGun_ = new G4ParticleGun(1);
			auto particleDefinition1 = G4ParticleTable::GetParticleTable()->FindParticle("geantino"); //"pi0"geantino#chargedgeantino
			fParticleGun_->SetParticleDefinition(particleDefinition1);
			fParticleGun_->SetParticleMomentum(G4ThreeVector(pxPi_P_init, pyPi_P_init, pzPi_P_init));
			fParticleGun_->GeneratePrimaryVertex(anEvent);
		}
		//* particleGun end
		//* pythia8
		else if (fCurrentGeneratorName == "pythia8" or fCurrentGeneratorName == "hepmc")
		{
			fGentypeMap[fCurrentGeneratorName]->GeneratePrimaryVertex(anEvent);
		}
		//* Pythia8 end
	}
	else
		G4Exception("PrimaryGeneratorAction::GeneratePrimaries",
					"InvalidSetup", FatalException,
					"Generator is not instanciated.");
}
