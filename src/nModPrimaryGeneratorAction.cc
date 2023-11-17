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
// $Id: nModPrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file nModPrimaryGeneratorAction.cc
/// \brief Implementation of the nModPrimaryGeneratorAction class

#include "nModPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "UserDataInput.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

nModPrimaryGeneratorAction::nModPrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	fGamma(nullptr), fNeutron(nullptr),
	primaryParticleEnergy(0. * MeV),
	fParticleGun(0)
{
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);

	// default particle kinematic
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	fGamma = particleTable->FindParticle(particleName = "gamma");
	fNeutron = particleTable->FindParticle(particleName = "neutron");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nModPrimaryGeneratorAction::~nModPrimaryGeneratorAction()
{
	delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nModPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//this function is called at the begining of ecah event
	//

	// In order to avoid dependence of PrimaryGeneratorAction
	// on DetectorConstruction class we get Envelope volume
	// from G4LogicalVolumeStore.

	UserDataInput userData;
	G4double gammaPercentage = userData.GetGammaPercentage();
	vector<G4double> gammaEnergy = userData.GetGammaEnergy();
	vector<G4double> gammaNormSpectrum = userData.GetGammaNormSpectrum();
	vector<G4double> neutronEnergy = userData.GetNeutronEnergy();
	vector<G4double> neutronNormSpectrum = userData.GetNeutronNormSpectrum();

	G4ParticleDefinition* particle;
	G4double random = G4UniformRand();
	if (random < gammaPercentage)
	{
		particle = fGamma;
		primaryParticleEnergy = EnergySampling(gammaEnergy, gammaNormSpectrum);
	}
	else
	{
		particle = fNeutron;
		primaryParticleEnergy = EnergySampling(neutronEnergy, neutronNormSpectrum);
	}

	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleEnergy(primaryParticleEnergy);

	G4ThreeVector pstnvctr, drctnvctr;
	Sampling(pstnvctr, drctnvctr);
	fParticleGun->SetParticlePosition(pstnvctr);
	fParticleGun->SetParticleMomentumDirection(drctnvctr);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}


G4double nModPrimaryGeneratorAction::EnergySampling(vector<G4double> energy, vector<G4double> normSpectrum)
{
	G4double random = G4UniformRand();
	G4double sum = 0;
	G4int i = 0;
	while (sum < random)
	{
		sum += normSpectrum[i];
		i++;
	}
	primaryParticleEnergy = energy[i - 1] * MeV;
	return primaryParticleEnergy;
}

void nModPrimaryGeneratorAction::Sampling(G4ThreeVector& pstnvctr, G4ThreeVector& drctnvctr)
{
	G4String choice = UserDataInput::GetSourceType();
	if (choice == "Point")
	{
		G4double phi = twopi * G4UniformRand();
		G4double costheta = 2. * G4UniformRand() - 1.;
		G4double sintheta = sqrt(1 - pow(costheta, 2));
		drctnvctr = G4ThreeVector(sintheta * cos(phi), sintheta * sin(phi), costheta);

		pstnvctr = G4ThreeVector();
	}
	else if (choice = "Parallel") //平行束源
	{
		/*G4double pstnx = -0.5 * 1 * cm + 1 * cm * G4UniformRand();
		G4double pstny = -0.5 * 1 * cm + 1 * cm * G4UniformRand();
		pstnvctr = G4ThreeVector(pstnx, pstny, 0.);*/
		pstnvctr = G4ThreeVector();

		drctnvctr = G4ThreeVector(0., 0., -1.);
	}

	/*UserDataInput userData;

	G4double dtctrx = userData.GetDectorDimensionX();
	G4double dtctry = userData.GetDectorDimensionY();
	// G4double dtctrz = userdata.GetDectorDimensionZ();
	G4double sourceDetectorDistance = userData.GetDistanceOfSD();
	G4double pstns = 0.5 * sourceDetectorDistance; // + dtctrz / 2.;
	G4double pstnx, pstny;

	G4String choice = userData.GetSourceType();
	if (choice == "Point") // 点源
	{
		pstnx = 0.;
		pstny = 0.;
		pstnvctr = G4ThreeVector(pstnx, pstny, pstns);

		G4double theta;
		G4double phi = 2. * pi * G4UniformRand();
		G4double phi0 = atan(dtctry / dtctrx);
		G4double phi1 = 0.5 * pi + atan(dtctrx / dtctry);
		G4double maxalpha;
		if ((phi > 0 && phi <= phi0) || (phi > phi1 && phi <= phi0 + pi) ||(phi > phi1 + pi && phi <= 2. * pi))
			maxalpha = atan(0.5 * dtctrx * sqrt(1 + tan(phi) * tan(phi)) / sourceDetectorDistance);
		else
			maxalpha = atan(0.5 * dtctry * sqrt(1 + tan(pi * 0.5 - phi) * tan(pi * 0.5 - phi)) / sourceDetectorDistance);

		G4double mintheta = pi - maxalpha;
		theta = mintheta + (pi - mintheta) * G4UniformRand();

		drctnvctr = G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
	}
	else if (choice == "Surface") //面源
	{
		// 粒子位置和方向的抽样取决于面源的大小
		// 只适用于面源面积小于探测器面积的情况,面源面积若大于探测器面积有可能陷入死循环
		// 先抽取phi，0到2pi；再抽取theta，pi/2到pi
		// 计算此方向所在直线与探测器所在平面的交点是否在探测器上，不在则舍弃，重新抽取theta
		// 抽样效率较低

		G4double srfcsrsx = 0.5 * dtctrx;
		G4double srfcsrsy = 0.5 * dtctry;
		pstnx = -0.5 * srfcsrsx + srfcsrsx * G4UniformRand();
		pstny = -0.5 * srfcsrsy + srfcsrsy * G4UniformRand();
		pstnvctr = G4ThreeVector(pstnx, pstny, pstns);

		G4double theta;
		G4double phi = 2. * pi * G4UniformRand();
		theta = 0.5 * pi + 0.5 * pi * G4UniformRand(); // pi/2到pi

		G4double xtemp = pstnx - sourceDetectorDistance * tan(theta) * cos(phi);
		G4double ytemp = pstny - sourceDetectorDistance * tan(theta) * sin(phi);
		while (abs(xtemp) > 0.5 * dtctrx || abs(ytemp) > 0.5 * dtctry)
		{
			theta = 0.5 * pi + 0.5 * pi * G4UniformRand();
			xtemp = pstnx - sourceDetectorDistance * tan(theta) * cos(phi);
			ytemp = pstny - sourceDetectorDistance * tan(theta) * sin(phi);
		}

		drctnvctr = G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
	}*/
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
