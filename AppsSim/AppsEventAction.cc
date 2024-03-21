#include "AppsEventAction.hh"
#include "AppsRunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include <iostream>
#include "AppsInput.hh"
#include "G4SDManager.hh"
#include "AppsHit.hh"
#include "G4Threading.hh"
#include "Randomize.hh"
#include <chrono>

using namespace std;

AppsEventAction::AppsEventAction(AppsRunAction *)
    : G4UserEventAction()
{
}

AppsEventAction::~AppsEventAction()
{
}

void AppsEventAction::BeginOfEventAction(const G4Event *)
{
}

void AppsEventAction::EndOfEventAction(const G4Event *)
{
}
