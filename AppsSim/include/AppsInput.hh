#ifndef AppsInput_h
#define AppsInput_h 1

#include "globals.hh"
#include "TH1.h"
#include "G4RunManager.hh"
#include "AppsGeometryObject.hh"
#include <fstream>
#include <vector>
#include <atomic>

using namespace std;

class AppsInput

{
public:
	AppsInput();
	~AppsInput();

	void ReadGeometryObjects(); // In charge of parsing the input geometry file
	void ReadRunParameters();	// In charge of parsing the input parameter file
	void ReadNeutronSpectrum();	// In charge of parsing the input spectrum file
	void DisplayMenu();			// In charge of the used interactive menu
	void PrintSimParameters();	// Function to display the simulation  parameters

	void PrintRunParameters(); // Function to display the individual run parameters

	void PrintGeometryParameters(); // function to print the objects parameters

	void Initialize();			 // Function called at initialization to read the input files, initialize parameters
	void InitializeParameters(); // function to initialize the input parameters

	void CheckSimParameters(); // function to check if the input parameters are correct

	void CheckGeometryObjects();

	// Gamma energy parameters

	void AddGammaEnergy(G4double energy) { GammaEnergy.push_back(energy); }
	vector<G4double> GetGammaEnergy() { return GammaEnergy; }

	void AddGammaBandwidth(G4double bandwidth) { GammaBandwidth.push_back(bandwidth); }
	vector<G4double> GetGammaBandwidth() { return GammaBandwidth; }

	void SetEnergyIncrement(G4double increment) { EnergyIncrement = increment; }
	G4double GetEnergyIncrement() { return EnergyIncrement; }

	// Event and run parameters

	void SetNumberOfRuns(G4double rnumber) { NumberOfRuns = rnumber; }
	G4int GetNumberOfRuns() { return NumberOfRuns; }

	void AddNumberOfEvents(G4double enumber) { NumberOfEvents.push_back(enumber); }
	vector<long long int> GetNumberOfEvents() { return NumberOfEvents; }

	void SetInterfaceMode(G4bool interface) { InterfaceMode = interface; }
	G4bool GetInterfaceMode() { return InterfaceMode; }

	void SetNumberOfCores(G4int nCores) { NumberOfCores = nCores; }
	G4int GetNumberOfCores() { return NumberOfCores; }

	// Beam parameters

	void SetPolarAngle(G4double min, G4double max)
	{
		PolarAngle.push_back(min);
		PolarAngle.push_back(max);
	}
	vector<G4double> GetPolarAngle() { return PolarAngle; }

	void SetAzimuthalAngle(G4double min, G4double max)
	{
		AzimuthalAngle.push_back(min);
		AzimuthalAngle.push_back(max);
	}
	vector<G4double> GetAzimuthalAngle() { return AzimuthalAngle; }

	void SetSourcePosition(G4double x, G4double y, G4double z)
	{
		SourcePosition.setX(x);
		SourcePosition.setY(y);
		SourcePosition.setZ(z);
	}
	G4ThreeVector GetSourcePosition() { return SourcePosition; }

	// Geometry objects

	void AddGeometryObject(AppsGeometryObject object) { GeometryObjects.push_back(object); }
	vector<AppsGeometryObject> GetGeometryObjects() { return GeometryObjects; }

	static AppsInput *Instance();

	vector<AppsGeometryObject> GeometryObjects; // holds the objects from the geometry

	// implemented as GBS

	void SetBeamType(G4String beamtype) { BeamType = beamtype; }
	G4String GetBeamType() { return BeamType; }

	// source parameters

	void AddActivityParameter(G4double activityPar) { ActivityParameters.push_back(activityPar); }
	void ClearActivityParameters() { ActivityParameters.clear(); }
	vector<G4double> GetActivityParameters() { return ActivityParameters; }

	void SetDecayNucleus(G4String nucleus) { DecayNucleus = nucleus; }
	G4String GetDecayNucleus() { return DecayNucleus; }

	void SetBunchParticles(G4double particles) { BunchParticles = particles; }
	G4double GetBunchParticles() { return BunchParticles; }

	TH1D* GetNeutronSpectrum() { return neutron_spectrum; }

	// active volumes

	G4int GetActiveNumber();

	// Time structure

	void AddTimeStructureParameter(G4int timeParameter) { TimeStructure.push_back(timeParameter); }
	void ClearTimeSTructureParameters() { TimeStructure.clear(); }
	vector<G4int> GetTimeStructureParameters() { return TimeStructure; }

	// beam Histogram type

	void AddBeamHistogramType(G4String hist) { BeamHistogramType.push_back(hist); }
	void SetBeamHistogramType(vector<G4String> hist) { BeamHistogramType = hist; }
	void ClearBeamHistogramType() { BeamHistogramType.clear(); }
	vector<G4String> GetBeamHistogramType() { return BeamHistogramType; }

	// Beam shape parameters

	void AddShapeParameter(G4double par) { ShapeParameters.push_back(par); }
	void ClearShapeParameters() { ShapeParameters.clear(); }
	vector<G4double> GetShapeParameters() { return ShapeParameters; }

	// Beam polarization

	void SetPolarization(G4double x, G4double y, G4double z)
	{
		Polarization.setX(x);
		Polarization.setY(y);
		Polarization.setZ(z);
	}
	G4ThreeVector GetPolarization() { return Polarization; }

private:
	G4int NumberOfCores = 0; // the number of cores

	G4int ActiveVolumes = 0;			 // the number of active volumes from the Geometry objects
	vector<G4double> ActivityParameters; // activity parameters of the source [0] activity [1] uncertainity
	G4String DecayNucleus;				 // the nucleus of the calibration source
	G4double BunchParticles;
	G4String BeamType; // string that holds the type of time structure -> GBS or calibration source

	static AppsInput *fgInstance;
	vector<G4double> GammaEnergy;		  // energy of the gamma beam
	vector<G4double> GammaBandwidth;	  // FWHM of the gamma beam
	G4double EnergyIncrement;			  // increment of the energy in multiple energy runs
	G4int NumberOfRuns;					  // keeps the number of runs
	vector<long long int> NumberOfEvents; // keeps the number of events
	G4bool InterfaceMode;				  // keeps display info
	vector<G4double> PolarAngle;		  // [0] min angle, [1] max angle
	vector<G4double> AzimuthalAngle;	  // [0] min angle, [1] max angle
	G4ThreeVector SourcePosition;		  // Threevector for source position
	G4int DefaultNameCount = 0;			  // number of box objects with undefined name

	vector<G4String> BeamHistogramType;
	vector<double> ShapeParameters;
	G4ThreeVector Polarization; // Threevector for source position
	vector<int> TimeStructure;
	TH1D *neutron_spectrum;

	// the list of available input parameter options
	vector<G4String> InputParametersOption = {"display interface", "number of runs", "number of events", "gamma energy", "bandwidth", "energy increment",
											  "polar angle", "azimuthal angle", "source position", "beam type", "bunch particles",
											  "decay activity", "decay nucleus", "number of cores", "time structure", "output", "shape parameters", "polarization"};

	// the list of available input objects
	vector<G4String> InputGeometryObjects = {"ELIADE", "Monitor", "CloverDetector", "Tube", "Box", "DiagDetector", "ScintDetector", "World", "Sphere", "Collimator", "PIXEL", "Mask", "CeBrDetector", "VEGA"};
	vector<G4String> InputObjectParameters = {"position", "rotation", "material", "mother volume", "name", "shape parameters", "active", "overlap", "histogram", "features", "movable"};
};

//, "detector distance","shield", "output", "extra", };

#endif
