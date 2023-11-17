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
// $Id: nModSteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file nModSteppingAction.cc
/// \brief Implementation of the nModSteppingAction class

#include "nModSteppingAction.hh"
#include "nModDetectorConstruction.hh"
#include "nModRun.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4SystemOfUnits.hh"

#include "UserDataInput.hh"

#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

/*
vector<G4double> nModSteppingAction::EnergyDepositOf2nd;// = *(new vector<G4double>);
vector<G4int> nModSteppingAction::trackIDs;// = *(new vector<G4int>);
*/

nModSteppingAction::nModSteppingAction(nModEventAction* eventAction, nModRunAction* runAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nModSteppingAction::~nModSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nModSteppingAction::UserSteppingAction(const G4Step* step)
{

  G4Track* theTrack = step->GetTrack();
  G4int parentID = theTrack->GetParentID();
  //if (parentID != 0) return;

  //G4float stepLength = step->GetStepLength();		//cm
  //G4StepPoint* prestepPoint = step->GetPreStepPoint();
  //G4StepPoint* poststepPoint = step->GetPostStepPoint();
  //G4StepStatus poststepStatus = poststepPoint->GetStepStatus();
  //G4float PreKineticEnergy = prestepPoint->GetKineticEnergy(); 	//MeV
  //G4float PostKineticEnergy = poststepPoint->GetKineticEnergy(); 	//MeV
  //G4cout << "PreKineticEnergy: " << PreKineticEnergy / keV << G4endl;
  //G4cout << "PostKineticEnergy: " << PostKineticEnergy / keV << G4endl;
  //G4cout << "stepLength: " << stepLength /um << G4endl;
  //G4ThreeVector prepoint = prestepPoint->GetPosition();
  //G4ThreeVector postpoint = poststepPoint->GetPosition();
  //G4cout << "prepoint: " << prepoint / um << G4endl;
  //G4cout << "postpoint: " << postpoint / um << G4endl << G4endl;

  //通过慢化剂后的中子能量
  G4StepPoint* poststepPoint = step->GetPostStepPoint();
  G4StepStatus poststepStatus = poststepPoint->GetStepStatus();
  G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
  G4String PhyVolName = touch->GetVolume()->GetName();
  G4String strPrtclName = theTrack->GetDefinition()->GetParticleName();
  G4float PostKineticEnergy = poststepPoint->GetKineticEnergy(); 	//MeV
  if (poststepStatus == fGeomBoundary)
  {
      if (PhyVolName == "Moderator" && strPrtclName == "neutron")
      {
          if (PostKineticEnergy <= 0.0253 * eV)
          {
              fRunAction->AddNeutronCount();
              //G4ThreeVector postPosition = poststepPoint->GetPosition() / cm;
              //G4double R
              //    = postPosition.getX() * postPosition.getX() + postPosition.getY() * postPosition.getY() + postPosition.getZ() * postPosition.getZ();
              //cout << "PostKineticEnergy: " << PostKineticEnergy / eV << "  R: " << R << endl;
              //getchar();
              //ofstream moderatedneutron;
              //G4String MDRTEDNTRN = "moderatedneutron.txt";
              //moderatedneutron.open(MDRTEDNTRN, ios::app);
              //moderatedneutron << PostKineticEnergy / eV << G4endl;
              //moderatedneutron.close();

              //G4ThreeVector outPosition = poststepPoint->GetPosition();
              ////G4double positionX = outPosition.getX();
              ////G4double positionY = outPosition.getY();
              ////G4double positionZ = outPosition.getZ();
              //////static nModRun* fnModRun = new nModRun();
              ////nModRun* fnModRun
              ////    = static_cast<nModRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
              ////fnModRun->PushBackOutPosition(positionX);
              ////fnModRun->PushBackOutPosition(positionY);
              ////fnModRun->PushBackOutPosition(positionZ);
              //
              //ofstream outPositionFile;
              //G4String outPositionFileNmae = "outPosition.txt";
              //outPositionFile.open(outPositionFileNmae, ios::app);
              //outPositionFile << outPosition << G4endl;
              //outPositionFile.close();

              //G4ThreeVector momentumDirection = poststepPoint->GetMomentumDirection();
              //ofstream outDirectionFile;
              //G4String outDirectionFileName = "outDirection.txt";
              //outDirectionFile.open(outDirectionFileName, ios::app);
              //outDirectionFile << momentumDirection << G4endl;
              //outDirectionFile.close();

              //G4double directionX = momentumDirection.getX();
              //G4double directionY = momentumDirection.getY();
              //G4double directionZ = momentumDirection.getZ();

              //G4cout << "step momentumDirection: " << momentumDirection
              //    << "  xyz: (" << directionX << ", " << directionY << ", " << directionZ << ")" << G4endl;
              //getchar();

              //nModRun* fnModRun
              //    = static_cast<nModRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
              //fnModRun->PushBackMomentumDirection(directionX);
              //fnModRun->PushBackMomentumDirection(directionY);
              //fnModRun->PushBackMomentumDirection(directionZ);

              //fEventAction->SetMyMomentumDirection(momentumDirection);
          }
      }
  }

  // 统计热中子的径迹总长度
  if (PhyVolName == "CalFlux" && strPrtclName == "neutron")
  {
      if (PostKineticEnergy <= 0.0253 * eV)
      {
          G4double stepLength = step->GetStepLength(); // cm
          fRunAction->AddNeutronTrack(stepLength);
          //cout << " PostKineticEnergy:" << PostKineticEnergy << endl << "stepLength: " << stepLength / um << endl;
          //getchar();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

