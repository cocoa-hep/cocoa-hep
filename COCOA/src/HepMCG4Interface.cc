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
/// \file eventgenerator/HepMC/HepMCEx01/src/HepMCG4Interface.cc
/// \brief Implementation of the HepMCG4Interface class
//
//
// #include "OutputRunAction.hh"

#include "HepMCG4Interface.hh"
#include "OutputRunAction.hh"

#include "G4RunManager.hh"
#include "G4LorentzVector.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include <vector>


HepMCG4Interface::HepMCG4Interface()
	: hepmcEvent(0)
{
	;
}

HepMCG4Interface::~HepMCG4Interface()
{
	delete hepmcEvent;
}


G4bool HepMCG4Interface::CheckVertexInsideWorld(const G4ThreeVector &pos) const
{
	G4Navigator *navigator = G4TransportationManager::GetTransportationManager()
								 ->GetNavigatorForTracking();

	G4VPhysicalVolume *world = navigator->GetWorldVolume();
	G4VSolid *solid = world->GetLogicalVolume()->GetSolid();
	EInside qinside = solid->Inside(pos);

	if (qinside != kInside)
		return false;
	else
		return true;
}

void HepMCG4Interface::HepMC2G4(const HepMC3::GenEvent *hepmcevt,
				G4Event *g4event)
{
	float eta_primary=-100;
	float phi_primary=-100;

	//* loop for vertex
	for ( HepMC3::ConstGenVertexPtr vitr : hepmcevt->vertices() ) {

		//* is vertex real ?
		G4bool qvtx = false;
		for ( HepMC3::ConstGenParticlePtr prt : vitr->particles_out() )
		{
		    if (!(prt->end_vertex()) && prt->status() == 1)
			{
				qvtx = true;
				break;
			}
		}
		if (!qvtx)
		{
		    for ( HepMC3::ConstGenParticlePtr prt : vitr->particles_out() )
			{ 
			    if(prt->status()  > 20 && prt->status() < 30)
				{
				    eta_primary =  ( prt->momentum().eta() );
				    phi_primary =  ( prt->momentum().phi() );
				}
			}
		    continue;
		}

		//* check world boundary
		HepMC3::FourVector pos = vitr->position();
		G4LorentzVector xvtx(pos.x(), pos.y(), pos.z(), pos.t());
		if (!CheckVertexInsideWorld(xvtx.vect() * mm))
		    continue;

		//* create G4PrimaryVertex and associated G4PrimaryParticles
		G4PrimaryVertex *g4vtx =
			new G4PrimaryVertex(xvtx.x() * mm, xvtx.y() * mm, xvtx.z() * mm,
								xvtx.t() * mm / c_light);

		std::vector<HepMC3::ConstGenParticlePtr> selected_particles;

		for ( HepMC3::ConstGenParticlePtr prt : vitr->particles_out() )
		{

			if (prt->status() != 1)
				continue;

			if ( config_json_var.fiducial_cuts.min_pT > 0.0 &&
			     prt->momentum().perp() < config_json_var.fiducial_cuts.min_pT )
			{
				continue;
			}
			if (fabs(prt->momentum().eta()) > config_json_var.max_eta_endcap )
			    continue;
	
			// Cut away everything separated from the primary eta,phi by more than dR (used for single-jet data)
			if(config_json_var.fiducial_cuts.dR_cut > 0)
			{
				float eta_cut = prt->momentum().eta();
				float phi_cut = prt->momentum().phi();
				float dphi = acos(cos(phi_cut-phi_primary));
				float deta = eta_cut-eta_primary;
				float dR = sqrt(dphi*dphi+deta*deta);
				if(dR > config_json_var.fiducial_cuts.dR_cut)
					continue;
			}
			
			// if the particle decayed "too far" into the detector, replace it with its parent. otherwise this function just returns the original particle.
			HepMC3::ConstGenParticlePtr pptr = m_truthrecordgraph.check_prod_location(prt);

			m_truthrecordgraph.add_to_vector(pptr, m_truthrecordgraph.m_interesting_particles);

			m_truthrecordgraph.add_all_moving_parents(pptr, m_truthrecordgraph.m_interesting_particles);

			G4int pdgcode = (pptr)->pdg_id();
			// skip neutrinos
			if (abs(pdgcode) == 12 || abs(pdgcode) == 14 || abs(pdgcode) == 16)
			{
				continue;
			}

			//add to selected particles (selected to be passed to GEANT), avoid duplicates.
			m_truthrecordgraph.add_to_vector(pptr, selected_particles);
		}

		for ( HepMC3::ConstGenParticlePtr prt : selected_particles )
		{
			HepMC3::FourVector mom = prt->momentum();
			int pdgcode = prt->pdg_id();
			G4LorentzVector p(mom.px(), mom.py(), mom.pz(), mom.e());
			G4PrimaryParticle *g4prim = new G4PrimaryParticle(pdgcode, p.x() * MeV, p.y() * MeV, p.z() * MeV);
			g4vtx->SetPrimary(g4prim);
			m_truthrecordgraph.m_final_state_particles.push_back(prt);
		}

		g4event->AddPrimaryVertex(g4vtx);
	}
	m_truthrecordgraph.fill_truth_graph();
}

HepMC3::GenEvent *HepMCG4Interface::GenerateHepMCEvent()
{
	HepMC3::GenEvent *aevent = new HepMC3::GenEvent();
	return aevent;
}

void HepMCG4Interface::GeneratePrimaryVertex(G4Event *anEvent)
{
	// delete previous event object
	delete hepmcEvent;

	// generate next event
	hepmcEvent = GenerateHepMCEvent();
	if (!hepmcEvent)
	{

		G4RunManager::GetRunManager()->AbortRun();
		return;
	}
	
	HepMC2G4(hepmcEvent, anEvent);
}
