#include "G4VSensitiveDetector.hh"

#include "AppsHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;



class AppsSD : public G4VSensitiveDetector
{
  public:
    AppsSD(const G4String& name, const G4String& hitsCollectionName);
    virtual ~AppsSD();
  

    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    AppsHitsCollection* fHitsCollection;
};


