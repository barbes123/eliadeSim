#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class AppsHit : public G4VHit
{
  public:
    AppsHit();
    AppsHit(const AppsHit&);
    virtual ~AppsHit();

    // operators
    const AppsHit& operator=(const AppsHit&);
    G4int operator==(const AppsHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();


    // Set methods
    void SetTrackID  (G4int track)      { fTrackID = track; }
    void SetChamberNb(G4int chamb)      { fChamberNb = chamb; }
    void SetEdep     (G4double de)      { fEdep = de; }
    void SetPos      (G4ThreeVector xyz){ fPos = xyz; }
   

    void SetTime     (G4double time)   {fGlobalTime = time; } 	// my add
    void SetFaceEnergy (G4double energy) {fFaceEnergy = energy; }	// my add
    void SetPosition (G4ThreeVector position){ fPosition = position; }
	


    // Get methods
    G4int GetTrackID() const     { return fTrackID; }
    G4int GetChamberNb() const   { return fChamberNb; }
    G4double GetEdep() const     { return fEdep; }
    G4ThreeVector GetPos() const { return fPos; }
    G4double GetTime() const     { return fGlobalTime; }
    G4double GetFaceEnergy() const { return fFaceEnergy; }
    G4ThreeVector GetPosition() const { return fPosition; }


  private:
	

	
		G4double fFaceEnergy;	
		G4ThreeVector fPosition;     		

	
      G4int         fTrackID;
      G4int         fChamberNb;
      G4double      fEdep;
      G4ThreeVector fPos; 
      G4double      fGlobalTime;
};



typedef G4THitsCollection<AppsHit> AppsHitsCollection;

extern G4ThreadLocal G4Allocator<AppsHit>* AppsHitAllocator;


inline void* AppsHit::operator new(size_t)
{
  if(!AppsHitAllocator)
      AppsHitAllocator = new G4Allocator<AppsHit>;
  return (void *) AppsHitAllocator->MallocSingle();
}



inline void AppsHit::operator delete(void *hit)
{
  AppsHitAllocator->FreeSingle((AppsHit*) hit);
}



