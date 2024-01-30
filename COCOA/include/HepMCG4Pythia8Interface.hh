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
/// \file eventgenerator/HepMC/HepMCEx03/include/HepMCG4Pythia8Interface.hh
/// \brief Definition of the HepMCG4Pythia8Interface class
//
//

#ifndef HEPMC_G4_PYTHIA8_INTERFACE_H
#define HEPMC_G4_PYTHIA8_INTERFACE_H

#include "HepMCG4Interface.hh"
#include "G4UImessenger.hh"

#include "Pythia8/Pythia8ToHepMC3.h"
#include "Pythia8/Pythia.h"
#include "Config_reader_var.hh"

class HepMCG4Pythia8Messenger;


class HepMCG4Pythia8Interface : public HepMCG4Interface
{
protected:
	Config_reader_var &config_var = Config_reader_var::GetInstance(); 
	G4int verbose;
	HepMC3::Pythia8ToHepMC3 ToHepMC;
	Pythia8::Pythia pythia;
	Pythia8::Event sum_events;

	HepMCG4Pythia8Messenger *messenger;

	virtual HepMC3::GenEvent *GenerateHepMCEvent();


public:
	HepMCG4Pythia8Interface();
	~HepMCG4Pythia8Interface();
	void fillParticle(int pdgid, double ee, double eta_init, double phi_init, int status, Pythia8::Event &event, Pythia8::ParticleData &pdt, bool atRest = false);
	// set/get methods
	void SetVerboseLevel(G4int i);
	G4int GetVerboseLevel() const;

	// call pyxxx
	void CallPythiaInit(); //
	void CallPythiaStat();
	void CallPythiaReadString(G4String par);

	Pythia8::Event GetPythiaObject();

	// random numbers operations
	void SetRandomSeed(G4int iseed);
	void PrintRandomStatus(std::ostream &ostr = G4cout) const; //
	void PrintRandomStatus() const;

	// setup user parameters (empty in default).
	// Implement your parameters in a delived class if you want.
	virtual void SetUserParameters();

	virtual void Print() const;
	int poisson(double nAvg, Pythia8::Rndm &rndm);
};

inline void HepMCG4Pythia8Interface::SetVerboseLevel(G4int i)
{
	verbose = i;
}

inline G4int HepMCG4Pythia8Interface::GetVerboseLevel() const
{
	return verbose;
}

class MyUserHooks : public Pythia8::UserHooks
//class MyUserHooks : public Pythia8::UserHooksPtr
{

public:
	MyUserHooks();
	~MyUserHooks();
	//* Allow a veto for the process level, to gain access to decays.
	bool canVetoProcessLevel();
	virtual bool canModifySigma();

	virtual double multiplySigmaBy(const Pythia8::SigmaProcess *sigmaProcessPtr, const Pythia8::PhaseSpace *phaseSpacePtr, bool /*inEvent*/);

	//* this is a veto at the event process and it is suited for user-veto definition within fiducial regions of the detector
	//* it is run during the main evolution scheme and it is this slower but more information is available
	bool doVetoProcessLevel(Pythia8::Event &process);

private:
	Pythia8::SlowJet *slowJet;
	double eff;
	double p;
	double f;
};

#endif
