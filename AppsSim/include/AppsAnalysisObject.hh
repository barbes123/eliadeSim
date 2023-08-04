#ifndef AppsAnalysisObject_h
#define AppsAnalysisObject_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include <vector>
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include <random>
#include "G4VHitsCollection.hh"
#include <tuple>
using namespace std;


class AppsAnalysisObject
{
  public:

	AppsAnalysisObject();
      	virtual ~AppsAnalysisObject();
		

	void Initialize();

	void SetObjectName(G4String name){ objectName = name; }						// set function for the name of the digitizer object
	G4String GetObjectName(){ return objectName; }							// get function for the name of the digitizer object
		
	tuple<double,double,int> ProcessEvent(G4VHitsCollection * hitCol);				// the hit collection from the event are passed to the analysis objects to create the tuple 
	void AddEvent(tuple<double, double, int> event);						// of the event
	


	void RecordSumEvents();
	void RecordDigitalEvents();
	void RecordAddBackEvents();
	
	void BuildDigitalEvents(unsigned long long currentTime);

	void SetDetectorResolution(vector<G4double> res){ detectorRes = res; }
	void SetDetectorResolution(G4double slope, G4double intercept){ detectorRes.clear(); detectorRes.push_back(slope); detectorRes.push_back(intercept); }
	vector<G4double> GetDetectorResolution(){ return detectorRes; }

	void SetFullBufferFlag(G4bool flag){ fullBuffer = flag; }
	G4bool GetFullBufferFlag(){ return fullBuffer; }

	void SetAddbackFlag(G4bool flag){ addbackFlag = flag; }
	G4bool GetAddbackFlag(){ return addbackFlag; }

	void SetDigitalBufferSize(G4int buffsize){ digitalBufferSize = buffsize; }
	G4int GetDigitalBufferSize(){ return digitalBufferSize; }
	
	void SetAddBackBufferSize(G4int buffsize){ addbackBufferSize = buffsize; }
	G4int GetAddBackBufferSize(){ return addbackBufferSize; }


	void FillDigitalHistograms(G4int index);
		
	void SetObjectHistograms(vector<G4String> hist){ objectHistograms = hist; }
	vector<G4String> GetObjectHistograms(){ return objectHistograms; }
	
	void SetFinalMergeFlag(G4bool flag){ finalMerge = flag; }
	G4bool GetFinalMergeFlag(){ return finalMerge; }

	G4double DetectorResolution(G4double energy);

	void ClearDigitalEvents();
	void ClearRawEvents();

	void SetDigitalFlag(G4bool flag){ digitalFlag = flag; }
	G4bool GetDigitalFlag(){ return digitalFlag; }

	G4bool WasHit = false;	

	vector<tuple<double, double, int>> digitalEvents;
	vector<tuple<double, double, int>> rawEvents;
	vector<tuple<double, double, int>> addBackEvents;

	unsigned int addbackBufferSize = 10;		
	unsigned int digitalBufferSize = 10;
	unsigned int digitalGlobalBufferSize = 2;

	void AddAddBackEvent(tuple<double, double, int>  event){ addBackEvents.push_back(event); }

	void AddAddBackEvents(vector<tuple<double, double, int>>  events);



	void SetAddBackEvents(vector<tuple<double, double, int>>  events){ addBackEvents = events; }
	vector<tuple<double, double, int>> GetAddBackEvents(){ return addBackEvents; }
	void ClearAddBackEvents(){ addBackEvents.clear(); }
	
	void SetDigitalEvents(vector<tuple<double, double, int>>  events){ digitalEvents = events; }	
	void AddDigitalEvents(vector<tuple<double, double, int>>  events);	

	void SetObjectId(G4double id){ objectId = id; }
	G4double GetObjectId(){ return objectId; }

	void SetMod(G4double id){ Mod = id; }
	G4double GetMod(){ return Mod; }

	void SetCh(G4double id){ Ch = id; }
	G4double GetCh(){ return Ch; }

	void SetLastTime(G4double time){ lastTime = time; }
	G4double GetLastTime(){ return lastTime; }


	void SetLinearResFlag(G4bool flag){ linearRes = flag; }
	G4bool GetLinearResFlag(){ return linearRes; }

	void AddRate(G4double rate){ theRates.push_back(rate); }
	vector<G4double> GetRates(){ return theRates; }
		
	
	vector<tuple<double, double, int>> sumEvents;

  protected:  
	
	vector<G4double> theRates;
	
	vector<G4String> objectHistograms;
	G4String objectName;
	G4bool primitiveFlag;

	G4double objectId;
	G4double Mod = -1;
	G4double Ch = -1;

	G4bool addbackFlag = false;	
	G4bool fullBuffer = false;	
	G4bool digitalFlag = false;
	G4bool finalMerge = false;
	
	


	G4bool linearRes = true;

	G4double lastTime = 0;

	vector<G4double> detectorRes; 


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

