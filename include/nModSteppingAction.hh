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
// $Id: nModSteppingAction.hh 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file nModSteppingAction.hh
/// \brief Definition of the nModSteppingAction class

#ifndef nModSteppingAction_h
#define nModSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include "nModEventAction.hh"
#include "nModRunAction.hh"

#include "UserDataInput.hh"
#include <vector>

using namespace std;

class nModEventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class nModSteppingAction : public G4UserSteppingAction
{
  public:
    nModSteppingAction(nModEventAction* eventAction, nModRunAction* runAction);
    virtual ~nModSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);
private:
    nModEventAction* fEventAction;
    nModRunAction* fRunAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
