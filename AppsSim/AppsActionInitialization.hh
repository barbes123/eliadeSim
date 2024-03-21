#ifndef AppsActionInitialization_h
#define AppsActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"
/// Action initialization class.

class AppsActionInitialization : public G4VUserActionInitialization
{
  public:
    AppsActionInitialization();
    virtual ~AppsActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
