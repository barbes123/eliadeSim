#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

#include "globals.hh"
#include "TF2.h"
#include "Randomize.hh"
#include "TF1.h"
#include "TH1.h"
#include <random>
#include "GammaGenerator.hh"

class G4ParticleGun;
class G4Event;
class G4Box;

class AppsPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	AppsPrimaryGeneratorAction();
	virtual ~AppsPrimaryGeneratorAction();
	virtual void GeneratePrimaries(G4Event *);
	const G4ParticleGun *GetParticleGun() const { return fParticleGun; }
	static AppsPrimaryGeneratorAction *Instance();
	G4double GetParticleTime();			  // sets up the particle time
	G4double GetParticleEnergy();		  // sets up the particle energy
	G4ThreeVector GetParticlePosition();  // sets up the particle position
	G4ThreeVector GetParticleDirection(); // sets up the particle direction

	std::vector<G4double> AngularDistrCo;
	std::vector<G4double> AngularIntCo;

	G4double GetRandom(std::vector<G4double> array);

	void GenerateGBSEvent(G4Event *anEvent);	// generates an event for the gbs type
	void GenerateSourceEvent(G4Event *anEvent); // generate a calibration source event
	void GeneratePuBeEvent(G4Event *anEvent); // generate a calibration PuBe source event
	void GenerateSimpleEvent(G4Event *anEvent); // generate a simple event
	void GenerateSurfEvent(G4Event *anEvent);	// generate particle on a flat surface
	void GenerateSpeckEvent(G4Event *anEvent);	// generate speckled particle on a flat surface
	void GenerateFunctionEvent(G4Event *anEvent);
	void GenerateVEGAEvent(G4Event *anEvent);

	static long long Macropulse;
	static long long Micropulse;
	static long long Time;

private:
	G4ParticleGun *fParticleGun;
	G4String BeamType;
	G4double BunchParticles;
	G4double MacroBunchParticles;
	G4double GammaEnergy;
	G4double GammaBandwidth;
	G4ThreeVector sourcePos;
	G4String DecayNucleus;
	G4double partEnergy;			// particle energy
	G4ThreeVector beamPolarization; // beam polarization

	std::vector<G4double> PolarAngle;								   // the polar angles of the beam - [0] min angle, [1] max angle
	std::vector<G4double> AzimuthalAngle;							   // the azimuthal angles of the beam - [0] min angle, [1] max angle
	std::vector<G4double> sourceActivity;
	std::vector<G4int> TimeStructure;
	std::vector<G4double> shapeParameters;

	TF1 *functionDistribution;
	GammaGenerator *vegaSource;

	std::piecewise_constant_distribution<> bremDistribution;
	std::mt19937 bremGenerator;

	double SurfR1 = 0;
	double SurfR2 = 0;
	// int totalSurf = G4RandFlat::shoot(50,100);
	int totalSurf = 1000;
	int refSurf = 5;
	int currentRun = 0;
};
