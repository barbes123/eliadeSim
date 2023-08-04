
#include "AppsDetectorConstruction.hh"
#include "AppsActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "AppsInput.hh"
#include "G4Run.hh"
#include "AppsRun.hh"
#include "AppsPhysicsList.hh"

int main(int argc, char **argv)
{

	AppsInput *FInput = AppsInput::Instance();

	G4UIExecutive *ui = 0;

	if (argc == 1)
	{
		ui = new G4UIExecutive(argc, argv);
	}

	G4Random::setTheEngine(new CLHEP::RanecuEngine);

	G4MTRunManager *runManager = new G4MTRunManager;

	runManager->SetNumberOfThreads(FInput->GetNumberOfCores());
	runManager->SetUserInitialization(new AppsDetectorConstruction());
	runManager->SetUserInitialization(new AppsPhysicsList());

	//    	G4VModularPhysicsList* physicsList = new QBBC;
	//	physicsList->SetVerboseLevel(1);
	//	runManager->SetUserInitialization(physicsList);

	//	delete physicsList;

	runManager->SetUserInitialization(new AppsActionInitialization());
	runManager->Initialize();

	G4int isVisual = FInput->GetInterfaceMode(); // visual flag

	if (isVisual)
	{

		G4VisManager *visManager = new G4VisExecutive;
		visManager->Initialize();
		G4UImanager *UImanager = G4UImanager::GetUIpointer();

		if (!ui)
		{

			G4String command = "/control/execute ";
			G4String fileName = argv[1];
			UImanager->ApplyCommand(command + fileName);
		}
		else
		{

			UImanager->ApplyCommand("/control/execute init_vis.mac");
			ui->SessionStart();
			delete ui;
		}
	}
	else
	{

		long long int totalEvents = FInput->GetNumberOfEvents()[0]; // number of events for single source type

		int maxInt = std::numeric_limits<int>::max() * 0.9;

		int numbRuns = totalEvents / maxInt;
		int lastRunEvents = totalEvents - ((numbRuns)*maxInt);

		for (int i = 0; i <= numbRuns; i++)
		{

			G4long seeds[2];
			time_t systime = time(NULL);
			seeds[0] = (long)systime;
			seeds[1] = (long)(systime * G4UniformRand());
			G4Random::setTheSeeds(seeds);

			if (i == numbRuns)
			{
				runManager->BeamOn(lastRunEvents);
			}
			else
			{
				runManager->BeamOn(maxInt);
			}
		}
	}

	/*string command = "hadd $(ls -t | head -n 1)/AppsSim.root $(ls -t | head -n 1)/AppsSim_*";
	int status = system(command.c_str());
	command = "cp $(ls -t | head -n 1)/AppsSim.root .";
	status = system(command.c_str());*/

	delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
