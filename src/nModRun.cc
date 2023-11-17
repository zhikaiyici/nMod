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
// $Id: B1Run.cc 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file B1Run.cc
/// \brief Implementation of the B1Run class

#include "nModRun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nModRun::nModRun()
	: G4Run(),
	//momentaDirection(0),
	outPosition(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nModRun::~nModRun()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nModRun::Merge(const G4Run* run)
{
	G4cout << "---------------------------Merge-----------------------" << G4endl;
	// getchar();

	const nModRun* localRun = static_cast<const nModRun*>(run);

	//std::list<G4double> localMomentaDirection = localRun->momentaDirection;
	//momentaDirection.merge(localMomentaDirection);

	std::list<G4double> localOutPosition = localRun->outPosition;
	outPosition.merge(localOutPosition);

	//G4cout << "size of eneryDeposit before merge: " << outPosition.size() << G4endl;
	//G4cout << "size of localEneryDeposit: " << localOutPosition.size() << G4endl;
	//for (std::list<G4double>::iterator itr = outPosition.begin(); itr != outPosition.end(); ++itr)
	//{
	//	if (*itr != 0)
	//	{
	//		G4cout << *itr << G4endl;
	//	}
	//}

	G4Run::Merge(run);

	//G4cout << "size of eneryDeposit after merge: " << outPosition.size() << G4endl;
	//G4cout << "size of localEneryDeposit: " << localOutPosition.size() << G4endl;
	//for (std::list<G4double>::iterator itr = outPosition.begin(); itr != outPosition.end(); ++itr)
	//{
	//	if (*itr != 0)
	//	{
	//		G4cout << *itr << G4endl;
	//	}
	//}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void nModRun::PushBackMomentumDirection(G4double mdir)
{
	momentaDirection.push_back(mdir); 
}*/

void nModRun::PushBackOutPosition(G4double outpstn)
{
	outPosition.push_back(outpstn);
	//for (std::list<G4double>::iterator itr = outPosition.begin(); itr != outPosition.end(); ++itr)
	//{
	//	G4cout << *itr << G4endl;
	//}
	//G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


