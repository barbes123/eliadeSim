#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "../include/NRFInteractions.hh"



class AppsPhysicsList: public G4VModularPhysicsList
{
  public:
    AppsPhysicsList();
   ~AppsPhysicsList();

  protected:
  

    	virtual void SetCuts();   
};




