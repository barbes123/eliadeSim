#include "AppsActionInitialization.hh"
#include "AppsPrimaryGeneratorAction.hh"
#include "AppsRunAction.hh"
#include "AppsEventAction.hh"
#include "AppsSteppingAction.hh"

#include "G4Threading.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AppsActionInitialization::AppsActionInitialization()
 : G4VUserActionInitialization()
{


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AppsActionInitialization::~AppsActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AppsActionInitialization::BuildForMaster() const
{


  AppsRunAction* runAction = new AppsRunAction;
  SetUserAction(runAction);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AppsActionInitialization::Build() const
{
	

	SetUserAction(new AppsPrimaryGeneratorAction);

  	AppsRunAction* runAction = new AppsRunAction;
  	SetUserAction(runAction);
  
  	AppsEventAction* eventAction = new AppsEventAction(runAction);
	SetUserAction(eventAction);
  
  //SetUserAction(new AppsSteppingAction(eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
