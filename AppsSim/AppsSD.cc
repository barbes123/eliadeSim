#include "AppsSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <iostream> 
#include "G4VProcess.hh"
using namespace std;



AppsSD::AppsSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}



AppsSD::~AppsSD() 
{}


void AppsSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new AppsHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}


G4bool AppsSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  // energy deposit

	
  	G4double edep = aStep->GetTotalEnergyDeposit();


	

	if (!(edep==0.0) || aStep->IsFirstStepInVolume()){ 
	
		AppsHit* newHit = new AppsHit();	

		newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  		newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
                                               ->GetCopyNumber());
  		newHit->SetEdep(edep);
  		newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

  		newHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
			
	  	newHit->SetFaceEnergy(aStep->GetPreStepPoint()->GetKineticEnergy());
		newHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());


  		fHitsCollection->insert( newHit );
	}


	//	


	
	

  return true;
}


void AppsSD::EndOfEvent(G4HCofThisEvent*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
