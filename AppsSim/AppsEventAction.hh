#ifndef AppsEventAction_h
#define AppsEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>


using namespace std;



class AppsRunAction;

/// Event action class
///

class AppsEventAction : public G4UserEventAction
{
  public:
    AppsEventAction(AppsRunAction* runAction);
    virtual ~AppsEventAction();


    	virtual void BeginOfEventAction(const G4Event* event);
    	virtual void EndOfEventAction(const G4Event* event);


  private:

	




};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
