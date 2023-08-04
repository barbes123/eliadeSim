#ifndef AppsGeometryObject_h
#define AppsGeometryObject_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include <vector>
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include <random>
#include "AppsAnalysisObject.hh"


using namespace std;

// Geometry object class to define all the geometry objects.

class AppsGeometryObject
{
  public:
      AppsGeometryObject(G4String type);
      AppsGeometryObject();
      virtual ~AppsGeometryObject();

	// Object name

	void SetObjectName(G4String name) { ObjectName = name; }
	G4String GetObjectName() { return ObjectName; }
	
	// Object mother volume

	void SetObjectMotherVolume(G4String name) { ObjectMotherVolume = name; }
	G4String GetObjectMotherVolume() { return ObjectMotherVolume; }

	// Object position 
	
	void SetObjectPosition(G4double x, G4double y, G4double z) { ObjectPosition.setX(x); ObjectPosition.setY(y); ObjectPosition.setZ(z); }
	void SetObjectPosition(G4ThreeVector pos){ ObjectPosition = pos; }
	G4ThreeVector GetObjectPosition() { return ObjectPosition; }
	
	// Object rotation

	void SetObjectRotation(G4double x, G4double y, G4double z);
	vector<G4double> GetObjectRotationAngle() { return ObjectRotationAngle; } 
	G4RotationMatrix * GetObjectRotation() { return ObjectRotation; }	

	// Object material 

	void AddObjectMaterial(G4String mat) { ObjectMaterial.push_back(mat); }
	void ClearObjectMaterial() { ObjectMaterial.clear(); }
	vector<G4String> GetObjectMaterial() { return ObjectMaterial; }

	// Object shape parameters

	void AddShapeParameter(G4double par){ ShapeParameters.push_back(par); }
	void ClearShapeParameters() { ShapeParameters.clear(); }
	vector<G4double> GetShapeParameters() { return ShapeParameters; }

	// Object active status

	void SetObjectActiveStatus(G4bool active) { ObjectActiveStatus = active; }
	G4bool GetObjectActiveStatus() { return ObjectActiveStatus; }	

	// Object output options 

	void AddOutputOptions(G4String option) { OutputOptions.push_back(option); }
	vector<G4String> GetOutputOptions() { return OutputOptions; }
	
	// Object class 
	
	void SetObjectType(G4String type) { ObjectType = type; }
	G4String GetObjectType(){ return ObjectType; }
	
	// Object overlap	

	void SetObjectOverlap(G4bool overlap) { ObjectOverlap = overlap; }
	G4bool GetObjectOverlap() { return ObjectOverlap; }	
	void Initialize();
	

	// Object Histogram type

	void AddObjectHistogramType(G4String hist){ ObjectHistogramType.push_back(hist); }
	void SetObjectHistogramType(vector<G4String> hist){ ObjectHistogramType = hist; }
	void ClearObjectHistogramType() { ObjectHistogramType.clear(); }	
	vector<G4String> GetObjectHistogramType(){ return ObjectHistogramType; }
		
	
	// Object type

	void AddObjectFeature(G4String hist){ ObjectFeatures.push_back(hist); }
	void ClearObjectFeatures() { ObjectFeatures.clear(); }	
	vector<G4String> GetObjectFeatures(){ return ObjectFeatures; }


	// Object ActiveVolumes

	void SetObjectActiveVolumes(vector<G4String> active) { ActiveVolumes = active; }
	void AddObjectActiveVolume(G4String active){ ActiveVolumes.push_back(active); }
	vector<G4String> GetObjectActiveVolumes(){ return ActiveVolumes; }	


	// Object detector resolution

	G4double AddDetectorResolution();
	void SetObjectFWHMParam(vector<double> fwhmPar){FWHMParameters = fwhmPar; }
	vector<G4double> GetObjectFWHMParam(){ return FWHMParameters; }
	
	// analysis objects

	//void AddAnalysisObject(AppsAnalysisObject * analysis){ AnalysisObjects.push_back(analysis); }
	vector<vector<AppsAnalysisObject *>> GetAnalysisObjects(){ return AnalysisObjects; }	
	vector<vector<AppsAnalysisObject *>> AnalysisObjects;


	// movable parametrs

	vector<G4double> GetMovingParameters(){ return MovingParameters; }	
	void SetMovingParameters(vector<G4double> mPar){ MovingParameters = mPar; }
	void AddMovingParameters(G4double mPar){ MovingParameters.push_back(mPar); }
		
	


  protected:  

	vector<G4double> FWHMParameters;			// slope [0] and intercept [1] from the detector FWHM 

	vector<G4String> ObjectFeatures;
	vector<G4String> ObjectHistogramType;
	vector<G4String> ActiveVolumes;

	G4ThreeVector ObjectPosition;
  	G4RotationMatrix * ObjectRotation;

	G4String ObjectName;
	G4String ObjectMotherVolume;	
	G4String ObjectType;

	G4bool ObjectOverlap;

	vector<G4String> ObjectMaterial;
	G4bool ObjectActiveStatus;
	vector<G4String> OutputOptions;
	vector<G4String> Extra;
		
	vector<G4double> ShapeParameters;
	vector<G4double> ObjectRotationAngle;

	vector<G4double> MovingParameters;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

