//*************************************************************
// G4 Physics Class --> NRFCrossSectionData: source file
// for Nuclear Resonance Fluorescence cross sections
// Created: N. Kikuzawa, Japan Atomic Energy Agency, 24-Nov-06
// Modified by Hani Negm to include the angular distribution
// for the polarization effect. 2012-2014
//*************************************************************

#include "../include/NRFCrossSectionData.hh"
#include <vector>
#include <sstream>
#include <algorithm>
#include <fstream>
#include "G4StokesVector.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4NistManager.hh"
#include <iostream>
#include <cmath>

using namespace std;

//-----------------------------------------------//

NRFCrossSectionData::NRFCrossSectionData(G4Isotope *anIsotope)
{

	Init(anIsotope);
}

NRFCrossSectionData::~NRFCrossSectionData()
{
}

G4bool NRFCrossSectionData::IsApplicable(const G4DynamicParticle *, const G4Isotope *anIsotope)
{

	if (anIsotope->GetZ() == theZ && anIsotope->GetN() == theN)
	{
		return true;
	}

	return false;
}

void NRFCrossSectionData::Init(G4Isotope *anIsotope)
{

	// gets isotope data

	theName = anIsotope->GetName();
	theN = anIsotope->GetN();
	theZ = anIsotope->GetZ();

	theLastSpin = 0;
	theLastParity = 0;
	hasCrossSectionData = false;

	// initialize a new level data

	theLevelData = new LevelData;

	// string that contains the name of the input file

	G4String dirName = getenv("PWD");
	dirName = dirName + "/CrossSection/";
	std::ostringstream streamName;
	streamName << dirName << "z" << theZ << ".a" << theN;
	G4String file = streamName.str();
	std::ifstream NRFCrossSectionDataFile(file);

	// warning if the files is not found and the function returns

	if (!NRFCrossSectionDataFile)
	{
		//  		G4cout << "NRFCrossSectionData -­‐ Warning: "
		//  		<< " cannot find cross section data file:" << "z" << theZ
		//		<< ".a" << theN << G4endl;
		return;
	}

	char inputChars[80] = {' '};
	G4String inputLine;

	while (!NRFCrossSectionDataFile.getline(inputChars, 80).eof())
	{

		inputLine = inputChars;
		inputLine = inputLine.strip(G4String::trailing);

		if (inputChars[0] != '#' && inputLine.length() != 0)
		{

			G4double a(0.0), b(0.0), d(0.0), e(0.0), f2(0.0), g1(0.0), h(0.0);
			G4String c(string{' '});
			std::istringstream tmpStream(inputLine);

			tmpStream >> a >> b >> c >> d >> e >> f2 >> g1 >> h;

			MyLevelData *aLevelData = new MyLevelData();

			aLevelData->Excite = a * keV;
			aLevelData->Spin = b;
			aLevelData->Width = d * eV;
			aLevelData->CrossSection = e * microbarn * MeV;
			aLevelData->F2 = f2;
			aLevelData->BranchingRatio = g1;
			aLevelData->ExcitedLevel = h * keV;

			if (c == (G4String)("+"))
			{
				aLevelData->Parity = 1;
			}
			else if (c == (G4String)("-"))
			{
				aLevelData->Parity = -1;
			}
			else
			{
				aLevelData->Parity = 0;
			}

			theLevelData->push_back(aLevelData);
		}

		if (theLevelData->size() > 0)
		{
			hasCrossSectionData = true;
		}
	}

	NRFCrossSectionDataFile.close();
	theGroundStateSpin = theLevelData->at(0)->Spin;
	theGroundStateParity = theLevelData->at(0)->Parity;

	// if (theLevelData->size() >0) PrintLevelData();
}

G4double NRFCrossSectionData::GetCrossSection(const G4DynamicParticle *aDynamicGamma)
{

	G4bool squareCS = false;

	theLastCrossSection = DBL_MIN;
	if (hasCrossSectionData == false)
		return theLastCrossSection;

	theLastMomentumDirection = aDynamicGamma->GetMomentumDirection();

	for (G4int i = 0; i < G4int(theLevelData->size()); i++)
	{

		if (squareCS)
		{ // condition for square cross section

			G4double deltaE = std::abs(aDynamicGamma->GetKineticEnergy() - theLevelData->at(i)->Excite) / eV;
			G4double HWHM = theLevelData->at(i)->Width * 0.5 / eV;

			if (deltaE < HWHM)
			{

				theLastSpin = theLevelData->at(i)->Spin;
				theLastParity = theLevelData->at(i)->Parity;
				theLastBranchingRatio = theLevelData->at(i)->BranchingRatio;
				theFirstExcitedLevel = theLevelData->at(i)->ExcitedLevel;
				theLastCrossSection = theLevelData->at(i)->CrossSection / theLevelData->at(i)->Width;

				return theLastCrossSection;
			}
		}
		else
		{ // else statement for gaussian cross section

			G4double deltaE = std::abs(aDynamicGamma->GetKineticEnergy() - theLevelData->at(i)->Excite) / eV;
			G4double sigma = (theLevelData->at(i)->Width / 2.354) / eV;

			if (deltaE <= 4 * sigma)
			{

				G4double gaus = exp(-(pow(deltaE, 2) / (2.0 * pow(sigma, 2))));

				theLastSpin = theLevelData->at(i)->Spin;
				theLastParity = theLevelData->at(i)->Parity;
				theLastBranchingRatio = theLevelData->at(i)->BranchingRatio;
				theFirstExcitedLevel = theLevelData->at(i)->ExcitedLevel;
				theLastCrossSection = (theLevelData->at(i)->CrossSection / theLevelData->at(i)->Width) * gaus;

				return theLastCrossSection;
			}
		}
	}

	return theLastCrossSection;
}

G4ThreeVector NRFCrossSectionData::GetMomentumDirection(const G4DynamicParticle *aParticle)
{
	return ComputeMomentumDirection(aParticle);
}

G4ThreeVector NRFCrossSectionData::ComputeMomentumDirection(const G4DynamicParticle *aParticle)
{
	//	G4double aSpin = 0;
	G4int aParity = 0;
	G4double F2;

	for (G4int i = 0; i < G4int(theLevelData->size()); i++)
	{

		G4double deltaE = std::abs(aParticle->GetKineticEnergy() - theLevelData->at(i)->Excite) / eV;
		G4double HWHM = theLevelData->at(i)->Width * 0.5 / eV;
		G4double sigma = (theLevelData->at(i)->Width / 2.354) / eV;

		if (deltaE <= 4 * sigma)
		{
			//	aSpin = theLevelData->at(i)->Spin - theGroundStateSpin;
			aParity = theLevelData->at(i)->Parity * theGroundStateParity;
			F2 = theLevelData->at(i)->F2;
		}
	}

	G4double greject(1.5);
	G4double cosTheta;
	G4ThreeVector aLastMomentumDirection = aParticle->GetMomentumDirection();
	G4StokesVector aPolarization(aParticle->GetPolarization());

	// NRF  transition
	G4ThreeVector aMomentumDirection = aParticle->GetMomentumDirection();

	do
	{
		aLastMomentumDirection = G4RandomDirection();
		cosTheta = std::cos(aLastMomentumDirection.getTheta());
		aLastMomentumDirection.rotateUz(aMomentumDirection);

		G4ThreeVector aPolarizationFrame = G4PolarizationHelper::GetFrame(aMomentumDirection, aPolarization);
		G4ThreeVector nInteractionFrame = G4PolarizationHelper::GetFrame(aMomentumDirection, aLastMomentumDirection);
		G4double cosPhi1 = nInteractionFrame * aPolarizationFrame;
		G4double Phi1 = std::acos(cosPhi1);

		G4double P2 = 0.5 * (3 * cosTheta * cosTheta - 1);
		G4double P22 = 0.0;

		if (std::abs(aPolarization.mag() - 1.) < 1.e-6)
		{
			P22 = 3 * (1 - cosTheta * cosTheta);
		}
		else
		{
			P22 = 0.0;
		}

		greject = 1 + F2 * F2 * (P2 + aParity * 0.5 * cos(2 * Phi1) * P22);

		if (greject < 0 || greject > 1.5)
		{
			G4cout << " ERROR ::" << G4endl;
		}
	} while (greject < 1.5 * G4UniformRand());

	return aLastMomentumDirection;
}

void NRFCrossSectionData::DumpPhysicsTable(const G4ParticleDefinition &)
{
	//	PrintLevelData();
}

void NRFCrossSectionData::PrintLevelData()
{

	G4cout << "-­‐-­‐-­‐-­‐-­‐ Material = " << theName << " -­‐-­‐-­‐-­‐-­‐ " << G4endl;

	if (theLevelData->size() == 0)
	{
		G4cout << "*** Warning *** No level Data ";
		return;
	}

	for (G4int i = 0; i < G4int(theLevelData->size()); i++)
	{
		G4String aParity;
		if (theLevelData->at(i)->Parity > 0)
		{
			aParity = "+";
		}
		else if (theLevelData->at(i)->Parity < 0)
		{
			aParity = "‐";
		}
		else
		{
			aParity = " ";
		}

		G4cout << i << ":Level(keV)= " << std::setw(5) << theLevelData->at(i)->Excite / keV
			   << " " << theLevelData->at(i)->Spin << aParity
			   << ", Width(eV)=" << std::setw(5) << theLevelData->at(i)->Width / eV
			   << ", CrossSection(ubarn MeV)= "
			   << theLevelData->at(i)->CrossSection / microbarn / MeV
			   << ", F2=" << theLevelData->at(i)->F2
			   << ", Branching Ratio=" << theLevelData->at(i)->BranchingRatio
			   << ", The excited level=" << theLevelData->at(i)->ExcitedLevel
			   << G4endl;
	}

	G4cout << "-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐-­‐- " << G4endl;
}
