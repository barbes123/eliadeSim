#include "G4PhotoElastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4Elastic.hh"
#include "G4Electron.hh"
#include <iostream>

using namespace std;

G4PhotoElastic::G4PhotoElastic(const G4String &processName,
                               G4ProcessType type) : G4VEmProcess(processName, type),
                                                     isInitialised(false)
{
  SetStartFromNullFlag(false);
  SetBuildTableFlag(true);
  SetMinKinEnergyPrim(0.1 * MeV);
  SetSplineFlag(true);
}

G4PhotoElastic::~G4PhotoElastic()
{
}

G4bool G4PhotoElastic::IsApplicable(const G4ParticleDefinition &p)
{
  return (&p == G4Gamma::Gamma());
}

void G4PhotoElastic::InitialiseProcess(const G4ParticleDefinition *)
{
  if (!isInitialised)
  {
    isInitialised = true;

    if (!EmModel(0))
    {
      SetEmModel(new G4Elastic(), 0);
    }

    EmModel(0)->SetLowEnergyLimit(0.1 * MeV);

    EmModel(0)->SetHighEnergyLimit(20 * MeV);

    AddEmModel(0, EmModel(0));
  }
}

void G4PhotoElastic::PrintInfo()
{
}
