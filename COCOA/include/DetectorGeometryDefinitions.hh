#ifndef __DETECTORGEOMETRYDEFINITIONS_H__
#define __DETECTORGEOMETRYDEFINITIONS_H__

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
/// \file COCOA/include/DetectorGeometryDefinitions.hh
/// \brief Definition of detector constants used in COCOA project.

#include "G4SystemOfUnits.hh"

#include "globals.hh" // used for G4int, long double

enum class CaloIdx
{
    ECAL1,
    ECAL2,
    ECAL3,
    HCAL1,
    HCAL2,
    HCAL3
};

// Geometry parameters cell size of dR=0.1 in mm
constexpr long double dR_01 = 125 * cm;

constexpr long double Z_Source = 150 * cm;

// constexpr G4int CELL_NUMBER[kNLayers] = {kNoEM1_Xcell, kNoEM2_Xcell, kNoEM3_Xcell, kNoHCAL1_Xcell, kNoHCAL2_Xcell, kNoHCAL3_Xcell};

// Detector geometry

constexpr long double widthiron_add = 350 * um;
constexpr long double thiknesPix = 150 * um;
constexpr long double r_inn_trkPix0 = 39 * mm;
constexpr long double r_out_trkPix0 = r_inn_trkPix0 + thiknesPix;
constexpr long double r_inn_trkPix1 = 75 * mm;
constexpr long double r_out_trkPix1 = r_inn_trkPix1 + thiknesPix;
constexpr long double r_inn_trkPix2 = 155 * mm;
constexpr long double r_out_trkPix2 = r_inn_trkPix2 + thiknesPix;
constexpr long double r_inn_trkPix3 = 213 * mm;
constexpr long double r_out_trkPix3 = r_inn_trkPix3 + thiknesPix;
constexpr long double r_inn_trkPix4 = 271 * mm;
constexpr long double r_out_trkPix4 = r_inn_trkPix4 + thiknesPix;
constexpr long double thiknesStr = 320 * um;
constexpr long double r_inn_trkStr0 = 405 * mm;
constexpr long double r_out_trkStr0 = r_inn_trkStr0 + thiknesStr;
constexpr long double r_inn_trkStr1 = 562 * mm;
constexpr long double r_out_trkStr1 = r_inn_trkStr1 + thiknesStr;
constexpr long double r_inn_trkStr2 = 762 * mm;
constexpr long double r_out_trkStr2 = r_inn_trkStr2 + thiknesStr;
constexpr long double r_inn_trkStr3 = 1000 * mm;
constexpr long double r_out_trkStr3 = r_inn_trkStr3 + thiknesStr;
constexpr long double pos_EndCap_trkPix0 = 350 * mm;
constexpr long double pos_EndCap_trkPix1 = 420 * mm;
constexpr long double pos_EndCap_trkPix2 = 530 * mm;
constexpr long double pos_EndCap_trkPix3 = 670 * mm;
constexpr long double pos_EndCap_trkPix4 = 870 * mm;
constexpr long double pos_EndCap_trkPix5 = 1100 * mm;
constexpr long double pos_EndCap_trkPix6 = 1400 * mm;
constexpr long double pos_EndCap_trkPix7 = 2000 * mm;
constexpr long double pos_EndCap_trkPix8 = 2300 * mm;
constexpr long double pos_EndCap_trkPix9 = 2650 * mm;

constexpr long double pos_EndCap_trkStr0 = 1300 * mm;
constexpr long double pos_EndCap_trkStr1 = 1600 * mm;
constexpr long double pos_EndCap_trkStr2 = 1900 * mm;
constexpr long double pos_EndCap_trkStr3 = 2250 * mm;
constexpr long double pos_EndCap_trkStr4 = 2650 * mm;
constexpr int NumberOfIronLayers = 2;

#endif // __DETECTORGEOMETRYDEFINITIONS_H__