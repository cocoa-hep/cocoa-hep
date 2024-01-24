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
/// \file eventgenerator/HepMC/HepMCEx01/src/HepMCG4Reader.cc
/// \brief Implementation of the HepMCG4Reader class
//
//

#include "HepMCG4Reader.hh"
#include "HepMCG4ReaderMessenger.hh"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/Print.h"

#include <iostream>
#include <fstream>

using namespace HepMC3;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Reader::HepMCG4Reader()
    :  filename(""), verbose(0), i_first_event(0), hepmc3_reader(nullptr)
{
    // if ( filename != "" )
    // 	hepmc3_reader = new HepMC3::ReaderRootTree(filename.c_str());

  messenger= new HepMCG4ReaderMessenger(this);
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Reader::~HepMCG4Reader()
{
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Reader::Initialize()
{

  if ( filename == "" )
      return;
  
  // delete hepmc3_reader;
  
  hepmc3_reader = HepMC3::deduce_reader( filename );
  hepmc3_reader->skip( i_first_event );
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC3::GenEvent* HepMCG4Reader::GenerateHepMCEvent()
{

  GenEvent* evt = new GenEvent(Units::MEV,Units::MM);
  hepmc3_reader->read_event(*evt);
  evt->set_units(HepMC3::Units::MEV, HepMC3::Units::MM);
  	
  if(!evt) return 0; // no more event
  if(verbose>0) HepMC3::Print::content(*evt);

  return evt;
}
