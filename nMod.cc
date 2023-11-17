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
// $Id: exampleB1.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "nModDetectorConstruction.hh"
//#include "nModPhysicsList.hh"
#include "nModActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4ScoringManager.hh"
#include "G4UImanager.hh"
#include "QGSP_BIC_AllHP.hh"
#include "Shielding.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "UserDataInput.hh"
#include <cmath>
#include <list>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
	//std::list<G4ThreeVector> test(0);
	//G4ThreeVector avector = G4ThreeVector(1, 2, 3);
	//G4cout << "avector: " << avector << G4endl;
	//test.push_back(avector);
	//for (list<G4ThreeVector>::iterator itr = test.begin(); itr != test.end(); ++itr)
	//{
	//	G4cout << *itr << G4endl;
	//}
	//getchar();

	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4long seed = time(NULL);
	G4Random::setTheSeed(seed);

	UserDataInput userdata;
	userdata.ReadInputData();

	// Construct the default run manager
	//
#ifdef G4MULTITHREADED
	G4MTRunManager* runManager = new G4MTRunManager;
	G4int numberOfThreads = userdata.GetNumberOfThreads();
	//cout << "GetNumberOfThreads" << numberOfThreads << endl;
	runManager->SetNumberOfThreads(numberOfThreads);
#else
	G4RunManager* runManager = new G4RunManager;
#endif

	// Activate UI-command base scorer
	G4ScoringManager* scManager = G4ScoringManager::GetScoringManager();
	scManager->SetVerboseLevel(1);

	// Set mandatory initialization classes
	//
	// Detector construction
	runManager->SetUserInitialization(new nModDetectorConstruction());

	//// Physics list
	//runManager->SetUserInitialization(new nModPhysicsList());
	//G4cout << "nModPhysicsList" << G4endl;

	// Physics list
	G4VModularPhysicsList* physicsList = new QGSP_BIC_AllHP;
	//G4VModularPhysicsList* physicsList = new Shielding;
	G4double cutsValue = 10. * um;
	physicsList->SetDefaultCutValue(cutsValue);
	physicsList->SetVerboseLevel(1);
	runManager->SetUserInitialization(physicsList);

	// User action initialization
	runManager->SetUserInitialization(new nModActionInitialization());

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	G4bool uiStatus = userdata.GetStatusOfUI();
	if (uiStatus != FALSE)
	{
		// Detect interactive mode (if no arguments) and define UI session
		//
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);

		// Initialize visualization
		//
		G4VisManager* visManager = new G4VisExecutive;
		// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
		// G4VisManager* visManager = new G4VisExecutive("Quiet");
		visManager->Initialize();

		UImanager->ApplyCommand("/run/verbose 2");
		UImanager->ApplyCommand("/event/verbose 0");
		UImanager->ApplyCommand("/tracking/verbose 0");
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;

		// Job termination
		// Free the store: user actions, physics_list and detector_description are
		// owned and deleted by the run manager, so they should not be deleted 
		// in the main() program !
		delete visManager;
		delete runManager;
	}
	else
	{
		if (argc > 1)
		{
			// batch mode
			G4String command = "/control/execute ";
			G4String fileName = argv[1];
			UImanager->ApplyCommand(command + fileName);
		}
		else
		{
			// Get the pointer to the User Interface manager
			G4UImanager* UImanager = G4UImanager::GetUIpointer();
			UImanager->ApplyCommand("/run/verbose 2");
			UImanager->ApplyCommand("/event/verbose 0");
			UImanager->ApplyCommand("/tracking/verbose 0");

			// initialize G4 kernel
			runManager->Initialize();
		}
		G4int numberOfEvents = userdata.GetNumberOfEvents();
		runManager->BeamOn(numberOfEvents);

		delete runManager;

		G4cout << G4endl << "Press Enter to exit." << G4endl;
		getchar();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
