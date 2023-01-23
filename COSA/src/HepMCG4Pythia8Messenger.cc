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


#ifdef G4LIB_USE_PYTHIA8

#include "HepMCG4Pythia8Messenger.hh"




#include <sstream>
#include <fstream>
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"


HepMCG4Pythia8Messenger::HepMCG4Pythia8Messenger(HepMCG4Pythia8Interface* agen)
  : gen(agen)
{
  dir= new G4UIdirectory("/generator/pythia8/");
  dir-> SetGuidance("Commands for Pythia 8 event generation");

  verbose= new G4UIcmdWithAnInteger("/generator/pythia8/verbose",this);
  verbose-> SetGuidance("set verbose level");
  verbose-> SetParameterName("verboseLevel", false, false);
  verbose-> SetRange("verboseLevel>=0 && verboseLevel<=2");

  print= new G4UIcmdWithoutParameter("/generator/pythia8/print", this);
  print-> SetGuidance("print user information.");

  cpythiainit= new G4UIcommand("/generator/pythia8/init", this);
  cpythiainit-> SetGuidance("call INIT");
  G4UIparameter* beam= new G4UIparameter("beam particle (pdg code)", 's', false);
  cpythiainit-> SetParameter(beam);
  G4UIparameter* target= new G4UIparameter("target particle (pdg code)", 's', false);
  cpythiainit-> SetParameter(target);
  G4UIparameter* eCM= new G4UIparameter("energy of system in CM frame (GeV)", 'd', false);
  cpythiainit-> SetParameter(eCM);

  cpythiastat= new G4UIcmdWithoutParameter("/generator/pythia8/stat", this);
  cpythiastat-> SetGuidance("print statistics");

  cpythiaread= new G4UIcommand("/generator/pythia8/read",this);
  cpythiaread-> SetGuidance("call PYTHIAREAD");
  G4UIparameter* parameter= new G4UIparameter ("Parameter", 's', false);
  cpythiaread-> SetParameter(parameter);

  setUserParameters=
    new G4UIcmdWithoutParameter("/generator/pythia8/setUserParameters",this);
  setUserParameters->
    SetGuidance("Set user parameters in the Pythia common blocks");

  setSeed= new G4UIcmdWithAnInteger("/generator/pythia8/setSeed", this);
  setSeed-> SetGuidance("set initial seed.");

  printRandomStatus=
    new G4UIcmdWithAString("/generator/pythia8/printRandomStatus", this);
  printRandomStatus-> SetGuidance("print random number status.");
  printRandomStatus-> SetParameterName("filename", true, false);
  printRandomStatus-> SetDefaultValue("std::cout");

  quarkgluon = new G4UIcmdWithAnInteger("/generator/pythia8/QuarkGluon", this);
  quarkgluon->SetGuidance("set quark or gluon.");
  minEnergy = new G4UIcmdWithADouble("/generator/pythia8/minEnergy", this);
  minEnergy->SetGuidance("set min Energy.");
  maxEnergy = new G4UIcmdWithADouble("/generator/pythia8/maxEnergy", this);
  maxEnergy->SetGuidance("set max Energy.");
  minEta = new G4UIcmdWithADouble("/generator/pythia8/minEta", this);
  minEta->SetGuidance("set min Eta.");
  maxEta = new G4UIcmdWithADouble("/generator/pythia8/maxEta", this);
  maxEta->SetGuidance("set max Eta.");
}

HepMCG4Pythia8Messenger::~HepMCG4Pythia8Messenger()
{
  delete verbose;
  delete print;
  delete cpythiainit;
  delete cpythiastat;
  delete cpythiaread;
  delete setUserParameters;
  delete setSeed;
  delete printRandomStatus;
  delete quarkgluon;
  delete minEnergy;
  delete maxEnergy;
  delete minEta;
  delete maxEta;

  delete dir;
}

void HepMCG4Pythia8Messenger::SetNewValue(G4UIcommand* command,
                                          G4String newValues)
{

  if(command == verbose) {  // /verbose ...
    G4int level= verbose-> GetNewIntValue(newValues);
    gen-> SetVerboseLevel(level);

  } else if (command == print) { // /print ...
    gen-> Print();
    
  } else if (command == cpythiainit) { // /pythiainit ...
    gen-> CallPythiaInit();

  } else if (command == cpythiastat) { // /pythiastat ...
    gen-> CallPythiaStat();

  } else if (command == cpythiaread) { // /pythiaread ...
    G4String s= newValues;
    gen-> CallPythiaReadString(s);

  } else if (command == setUserParameters) { // /setUserParameters ...
    gen-> SetUserParameters();

  } else if (command == setSeed) { // /setSeed ...
    G4int iseed= setSeed-> GetNewIntValue(newValues);
    gen-> SetRandomSeed(iseed);

  } else if (command == printRandomStatus) { // /printRandomStatus ...
    G4String s= newValues;
    if (newValues == "std::cout") {
      // gen->PrintRandomStatus();
    } 
    else {
      // to a file (overwrite mode)
      std::ofstream ofs;
      ofs.open(s.c_str(), std::ios::out);
      //ofs.open(randomStatusFileName.c_str(), std::ios::out|std::ios::app);
      ofs.setf(std::ios::fixed | std::ios::showpoint);
      gen->PrintRandomStatus(ofs);
      ofs.close();
    }
  }
  else if (command==quarkgluon)
  {
    GQparticle = quarkgluon-> GetNewIntValue(newValues);
  }
  else if (command==minEnergy)
  {
    MinEnergy = minEnergy->GetNewDoubleValue(newValues);
  }
  else if (command==maxEnergy)
  {
    MaxEnergy = maxEnergy->GetNewDoubleValue(newValues);
  }
  else if (command==minEta)
  {
    MinEta = minEta->GetNewDoubleValue(newValues);
  }
  else if (command==maxEta)
  {
    MaxEta = maxEta->GetNewDoubleValue(newValues);
  }
}

G4String HepMCG4Pythia8Messenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  if (command == verbose) {
    cv= verbose-> ConvertToString(gen->GetVerboseLevel());
  }
  return cv;
}

#endif
