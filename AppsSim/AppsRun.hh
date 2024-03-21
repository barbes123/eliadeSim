#ifndef AppsRun_h
#define AppsRun_h 1

#include "G4Run.hh"
#include "AppsAnalysisObject.hh"
#include <vector>
#include "G4Threading.hh"
#include "G4RunManager.hh"

using namespace std;


class AppsRun : public G4Run
{
  public:

  	AppsRun();
  	virtual ~AppsRun();
  	virtual void RecordEvent(const G4Event*);
  	virtual void Merge(const G4Run*);

//	void BuildDigitalEvents(const G4int id);

	void BuildAnalysisObjects();

	unsigned int digitalBufferSize = 1000000;
	unsigned int rawEventsBufferSize = 1000000;
	unsigned int digitalGlobalBufferSize = 10;

//	G4double DetectorResolution(vector<G4double> fwhmPar, G4double energy);

//	void FillHistogram(G4int id);
//	void RecordEvents(G4int i, G4int o, G4int id);

	vector<vector<tuple<double, double, int>>> DigitalEvents;

	vector<vector<AppsAnalysisObject *>> AnalysisObjects;
	
	void SetFinalRunFlag(G4bool flag){ finalRun = flag; }
	G4bool GetFinalRunFlag(){ return finalRun; }

	vector<G4int> HitList;

	G4double beamTime = 0;

	G4int mergeTime = 1;
	G4int mergeIncrement = 1;
	G4bool finalRun = true;

	

	unsigned long long int lastMerge = 0;
	
	G4int nrCores = G4Threading::G4GetNumberOfCores();

	vector<G4int> mergeCounter;



 private:
  
};

#endif

