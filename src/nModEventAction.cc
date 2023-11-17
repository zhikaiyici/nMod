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
// $Id: nModEventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file nModEventAction.cc
/// \brief Implementation of the nModEventAction class

#include "nModRun.hh"
#include "nModEventAction.hh"
#include "nModRunAction.hh"
#include "nModSteppingAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

////vector<vector<G4double>> nModEventAction::EnergyDeposit;// = *(new vector<vector<G4double>>);
//list<vector<G4double>> nModEventAction::EnergyDeposit = *(new list<vector<G4double>>);

nModEventAction::nModEventAction(nModRunAction* /*runAction*/)
	: G4UserEventAction(),
	momentumDirection(G4ThreeVector())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nModEventAction::~nModEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nModEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
	momentumDirection = G4ThreeVector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nModEventAction::EndOfEventAction(const G4Event* event)
{
	//nModRun* fnModRun
	//	= static_cast<nModRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
	//if (momentumDirection != G4ThreeVector())
	//{
	//	G4double directionX = momentumDirection.getX();
	//	G4double directionY = momentumDirection.getY();
	//	G4double directionZ = momentumDirection.getZ();

	//	//G4cout << "event momentumDirection: " << momentumDirection
	//	//	<< "xyz: (" << directionX << ", " << directionY << ", " << directionZ << ")" << G4endl;
	//	//getchar();
	//	//fnModRun->PushBackMomentumDirection(directionX);
	//	//fnModRun->PushBackMomentumDirection(directionY);
	//	//fnModRun->PushBackMomentumDirection(directionZ);
	//}

	G4int eventID = event->GetEventID();
	G4int numberofevents = UserDataInput::GetNumberOfEvents();
	if (eventID == 0 || ((eventID + 1) % (numberofevents / 10) == 0))
	{
		if (eventID == 0)
		{
			G4int numofevents = UserDataInput::GetNumberOfEvents();
			G4double doubNumOfEvents = numofevents;
			G4double moderatorThickness = UserDataInput::GetModeratorThickness() / cm;
			G4cout << G4endl << " " << doubNumOfEvents << " event(s) will be simulated."
				   << G4endl
				   << " The thickness of moderator is " << moderatorThickness << " cm"
				   << G4endl
				   << G4endl;
		}
		G4int per =(int)( (1. * eventID + 1) / (numberofevents * 0.01));
		//G4cout << " numberofevents: "<< numberofevents << " eventID: "<< eventID <<G4endl;
		G4long seconds = time(NULL); // 格林威治时间
		seconds = seconds + 8 * 3600; // 北京时间
		G4int secondnow = seconds % 60;
		G4int minutes = (seconds - secondnow) / 60;
		G4int minutenow = minutes % 60;
		G4int hours = (minutes - minutenow) / 60;
		G4int hournow = hours % 24;
		G4cout << " Time now: " << setw(2) << hournow << ":" << setw(2) << minutenow << ":" << setw(2) << secondnow
			   << ". " << setw(3) << per << "% of simulation completed."
			   << G4endl;
		//getchar();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
