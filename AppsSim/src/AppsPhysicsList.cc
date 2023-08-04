#include "AppsPhysicsList.hh"
#include <iostream>
#include "NeutronHPphysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"

using namespace std;

AppsPhysicsList::AppsPhysicsList()
	: G4VUserPhysicsList()
{

	// Default physics
	RegisterPhysics(new G4DecayPhysics());

	//Radioactive decay
	RegisterPhysics(new G4RadioactiveDecayPhysics());

	// EM physics
	RegisterPhysics(new G4EmStandardPhysics());
	// Neutron Physics
	RegisterPhysics( new NeutronHPphysics("neutronHP"));  

}

AppsPhysicsList::~AppsPhysicsList()
{
}


void AppsPhysicsList::SetCuts()
{
	SetCutsWithDefault();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
