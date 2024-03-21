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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
//

#include <stdlib.h>
#include "AppsAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <G4Threading.hh>
#include "G4Run.hh"
#include "G4RunManager.hh"
#include <iostream>

G4ThreadLocal AppsAnalysisManager* AppsAnalysisManager::fgInstance = nullptr;
namespace {
    // Mutex to lock AppsAnalysisManager constructor
    G4Mutex AppsAnalysisManagerMutex = G4MUTEX_INITIALIZER;
} // namespace

AppsAnalysisManager::AppsAnalysisManager()
{
}

AppsAnalysisManager::~AppsAnalysisManager() 
{
}

AppsAnalysisManager* AppsAnalysisManager::Instance()
{
    if (fgInstance == nullptr) {
        static G4ThreadLocalSingleton<AppsAnalysisManager> inst;
        fgInstance = inst.Instance();
    }
    return fgInstance;
}

void AppsAnalysisManager::Write()
{
	G4String fFilename = "AppsSim_" + std::to_string(G4Threading::G4GetThreadId()) + ".root";
	fFile = new TFile(fFilename.c_str(), "RECREATE");
	ftree->Write();
}

void AppsAnalysisManager::CloseFile()
{
	fFile->Close();
}

void AppsAnalysisManager::CreateELIADETree()
{
	ftree = new TTree("ELIADE_Tree", "ELIADE_Tree");
	ftree->Branch("Mod", &Mod, "Mod/b");
	ftree->Branch("Ch", &Ch, "Ch/b");
	ftree->Branch("TimeStamp", &TimeStamp, "TimeStamp/l");
	ftree->Branch("FineTS", &FineTS, "FineTS/D");
	ftree->Branch("ChargeLong", &ChargeLong, "ChargeLong/s");
	ftree->Branch("ChargeShort", &ChargeShort, "ChargeShort/s");
	ftree->Branch("RecordLength", &RecordLength, "RecordLength/i");
	ftree->Branch("Signal", &Signal, "Signal");
}

void AppsAnalysisManager::FillELIADETreeColumn(G4int id, G4double value)
{
	switch (id)
	{
	
	case 0://Mod
		Mod = UChar_t(value);
		break;
	case 1://Ch
		Ch = UChar_t(value);
		break;
	case 2://TimeStamp
		TimeStamp = ULong64_t(value);
		break;
	case 3://FineTS
		FineTS = Double_t(value);
		break;
	case 4://ChargeLong
		ChargeLong = UShort_t(value);
		break;
	case 5://ChargeShort
		ChargeShort = UShort_t(value);
		break;
	case 6://RecordLength
		RecordLength = UInt_t(value);
		break;
	}
}

void AppsAnalysisManager::AddELIADETreeRow()
{
	ftree->Fill();
	Mod = ' ';
	Ch = ' ';
	TimeStamp = 0;
	FineTS = 0.;
	ChargeLong = 0;
	ChargeShort = 0;
	RecordLength = 0;
}


