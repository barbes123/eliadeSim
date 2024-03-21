#include "AppsRunAction.hh"
#include "AppsPrimaryGeneratorAction.hh"
#include "AppsDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "AppsAnalysisManager.hh"
#include <iostream>
#include "AppsInput.hh"
#include "G4SDManager.hh"
#include "tbb/concurrent_vector.h"
#include "AppsRun.hh"
#include <chrono>
#include <ctime>
#include <iomanip>

using namespace std;

AppsRunAction::AppsRunAction()
	: G4UserRunAction()
{
}

G4Run *AppsRunAction::GenerateRun()
{

	return new AppsRun();
}

AppsRunAction::~AppsRunAction()
{
}

void AppsRunAction::BeginOfRunAction(const G4Run *run)
{

	const AppsRun *localRun = static_cast<const AppsRun *>(run);

	if (IsMaster() && run->GetRunID() == 0)
	{

		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);
		std::ostringstream oss;
		oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");

		outFileName = "AppsSim_" + oss.str();

		string command = "mkdir " + outFileName;
		G4int status = system(command.c_str());

		command = "cp InputObjects.txt " + outFileName + "/";
		status = system(command.c_str());

		command = "cp SimParameters.txt " + outFileName + "/";
		status = system(command.c_str());
	}

	G4RunManager::GetRunManager()->SetRandomNumberStore(false);

	// Get analysis manager

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	// open the output file

	G4String fileName = "AppsSim_" + to_string(run->GetRunID()) + ".root";

	// Get input manager

	AppsInput *FInput = AppsInput::Instance();

	// ELIADE tree

	if (IsMaster() && run->GetRunID() == 0)
	{
		analysisManager->CreateELIADETree();
	}

}

void AppsRunAction::EndOfRunAction(const G4Run *theRun)
{

	const AppsRun *localRun = static_cast<const AppsRun *>(theRun);

	// G4cout << " core at the end " << endl;

	if (IsMaster())
	{

		for (unsigned int i = 0; i < localRun->AnalysisObjects.size(); i++)
		{

			for (unsigned int o = 0; o < localRun->AnalysisObjects[i].size(); o++)
			{

				localRun->AnalysisObjects[i][o]->BuildDigitalEvents(0);
				localRun->AnalysisObjects[i][o]->RecordDigitalEvents();
			}
		}

		AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

		analysisManager->Write();
		analysisManager->CloseFile();
	}

	if (IsMaster())
	{

		string command = "hadd AppsSim_" + to_string(theRun->GetRunID()) + ".root " + "AppsSim_*.root";
		int status = system(command.c_str());

		command = "mv AppsSim_" + to_string(theRun->GetRunID()) + ".root " + outFileName + "/";
		status = system(command.c_str());

		command = "rm AppsSim_*.root ";
		status = system(command.c_str());
	}

	delete AppsAnalysisManager::Instance();
}
