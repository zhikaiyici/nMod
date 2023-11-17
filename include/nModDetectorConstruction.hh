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
// $Id: nModDetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file nModDetectorConstruction.hh
/// \brief Definition of the nModDetectorConstruction class

#ifndef nModDetectorConstruction_h
#define nModDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "UserDataInput.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class nModDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    nModDetectorConstruction();
    virtual ~nModDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  protected:
    G4LogicalVolume*  fScoringVolume;

  private:
      // Option to switch on/off checking of volumes overlaps
      //
      G4bool checkOverlaps; // = true;
      UserDataInput* userData;
      /*G4double dtctrx = userData.GetDectorDimensionX();
      G4double dtctry = userData.GetDectorDimensionY();
      G4double dtctrz = userData.GetDectorDimensionZ();
      G4double sourceDetectorDistance = userData.GetDistanceOfSD();*/
      G4double moderatorThickness;// = userData.GetModeratorThickness();

      void ConstructModerator(G4LogicalVolume* logicWorld);
      void ConstructDetector(G4LogicalVolume* logicWorld);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

