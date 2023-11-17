﻿//
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
// $Id: nModRunAction.hh 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file nModRunAction.hh
/// \brief Definition of the nModRunAction class

#ifndef nModRunAction_h
#define nModRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

#include <vector>
#include <list>

using namespace std;

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class nModRunAction : public G4UserRunAction
{
public:
	nModRunAction();
	virtual ~nModRunAction();

	// virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void   EndOfRunAction(const G4Run*);

	inline void AddNeutronCount() { fNeutronCount += 1; /* G4int neutronCount = fNeutronCount.GetValue(); G4cout << neutronCount << G4endl; */ }
	inline void AddNeutronTrack(G4double stepLength) { fNeutronTrack += stepLength; }

private:
	//static G4int neutronCount;
	//static G4double neutronTrack;
	G4Accumulable<G4int> fNeutronCount;
	G4Accumulable<G4double> fNeutronTrack;
};

#endif
