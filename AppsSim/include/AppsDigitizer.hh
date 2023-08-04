#ifndef AppsDigitizer_h
#define AppsDigitizer_h 1

#include "globals.hh"
#include <vector>
#include "tbb/concurrent_vector.h"

using namespace std;


class AppsDigitizer
{
  public:
      	AppsDigitizer();
      	virtual ~AppsDigitizer();


	tbb::concurrent_vector<tuple<double, double, int>> StoredEvents;

	void AddEvent(G4double time, G4double energy, G4int pileupflag = 1);
	
	tbb::concurrent_vector<tuple<double, double, int>> GetStoredEvents(){ return StoredEvents; }
	
	void BuildDigitalEvent();
	void ClearStoredEvents();

  protected:  
	

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

