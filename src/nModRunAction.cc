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
// $Id: nModRunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file nModRunAction.cc
/// \brief Implementation of the nModRunAction class

#include "nModRunAction.hh"
#include "nModPrimaryGeneratorAction.hh"
#include "nModDetectorConstruction.hh"
#include "nModRun.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicalConstants.hh"
#include "nModSteppingAction.hh"
#include "nModEventAction.hh"
#include "UserDataInput.hh"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
/*#include <direct.h>
#include <io.h>*/

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*G4Accumulable<G4int> nModRunAction::fNeutronCount = 0;
G4Accumulable<G4double> nModRunAction::fNeutronTrack = 0.;*/
//G4int nModRunAction::neutronCount = 0;
//G4double nModRunAction::neutronTrack = 0.;

nModRunAction::nModRunAction()
	: G4UserRunAction(),
	fNeutronCount(0),
	fNeutronTrack(0.)
{
	// Register accumulable to the accumulable manager
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->RegisterAccumulable(fNeutronCount);
	accumulableManager->RegisterAccumulable(fNeutronTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nModRunAction::~nModRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nModRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	/*// inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);*/

	// reset accumulables to their initial values
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Reset();

	//neutronCount = 0;
	//neutronTrack = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nModRunAction::EndOfRunAction(const G4Run* run)
{
	G4int nofEvents = run->GetNumberOfEvent();
	if (nofEvents == 0) return;

	// Merge accumulables
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Merge();

	/*// 判断输出目标文件夹是否存在，不存在则创建
	if (_access("output", 4) != 0)
	{
		_mkdir("output");
		//G4cout << "已创建output" << G4endl;
	}*/

	// Print
	//  
	if (IsMaster()) {
		G4double gammaPercentage = UserDataInput::GetGammaPercentage();
		//G4int numberOfNeutron = (1. - gammaPercentage) * nofEvents;

		G4double numofevent = nofEvents;
		ostringstream ostrsEventNumber, ostrsGammaPercentage;
		ostrsEventNumber << setprecision(1) << numofevent;
		ostrsGammaPercentage << gammaPercentage;
		G4String strEventNumber = ostrsEventNumber.str();
		G4String strGammaPercentage = ostrsGammaPercentage.str();

		G4double moderatorThickness = UserDataInput::GetModeratorThickness() / cm;
		ostringstream ostrsModeratorThickness; ostrsModeratorThickness << moderatorThickness;
		G4String strModeratorThickness = ostrsModeratorThickness.str();

		G4double neutronTrack = fNeutronTrack.GetValue();
		G4double neutronFluence;
		G4String sourceType = UserDataInput::GetSourceType();
		if (sourceType == "Point")
		{
			G4double volume = 4 * pi * (pow(moderatorThickness + 1. * mm / cm, 3) - pow(moderatorThickness, 3)) / 3.;
			neutronFluence = neutronTrack / volume;
		}
		else if (sourceType == "Parallel")
		{
			G4double volume = 4 * pi * (pow(moderatorThickness + 1. * mm / cm, 3) - pow(moderatorThickness, 3)) / 3.;
			neutronFluence = neutronTrack / volume;
			//neutronFluence = neutronTrack / (1 * cm3);
		}

		G4int neutronCount = fNeutronCount.GetValue();
		ofstream numOfThermalNeutron, fileMomentumDirection;
		G4String fileName =
			"output/nmod_" + strEventNumber + "(gamma_" + strGammaPercentage + "_" + sourceType + ").txt";
		numOfThermalNeutron.open(fileName, ios::app);
		numOfThermalNeutron << moderatorThickness << setw(20) << neutronFluence << setw(20) << neutronCount << G4endl;
		numOfThermalNeutron.close();

		//const nModRun* fnModRun = static_cast<const nModRun*>(run);
		//list<G4double> momentaDirection = fnModRun->GetMomentaDirection();
		//G4String directionFileName = 
		//	"output/direction_" + strEventNumber + "(gamma_" + strGammaPercentage + "_" + sourceType + ").txt";
		//fileMomentumDirection.open(directionFileName, ios::out);
		//for (list<G4double>::iterator itr = momentaDirection.begin(); itr != momentaDirection.end(); ++itr)
		//{
		//	fileMomentumDirection << *itr << G4endl;
		//}
		//fileMomentumDirection.close();

		G4cout
			<< G4endl
			<< "--------------------End of Global Run-----------------------";
		G4cout
			<< G4endl
			<< " The number of thermal neutron is " << neutronCount
			<< G4endl
			<< " The neutron fluence is " << neutronFluence << " cm-2"
			<< G4endl
			<< "------------------------------------------------------------"
			<< G4endl
			<< G4endl;
	}
	else {
		G4cout
			<< G4endl
			<< "--------------------End of Local Run------------------------"
			<< G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
