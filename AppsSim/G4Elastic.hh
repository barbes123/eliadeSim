#include "G4VEmModel.hh"
#include <vector>

using namespace std;

struct Amplitude;
struct TAmplitude;

class G4ParticleChangeForGamma;

class G4Elastic : public G4VEmModel
{

public:
	G4Elastic(const G4ParticleDefinition *p = 0, const G4String &nam = "Elastic Scattering");

	virtual ~G4Elastic();

	virtual void Initialise(const G4ParticleDefinition *, const G4DataVector &);

	G4double ComputeCrossSectionPerAtom(double ekin, G4Isotope *isotope);

	virtual void SampleSecondaries(std::vector<G4DynamicParticle *> *,
								   const G4MaterialCutsCouple *,
								   const G4DynamicParticle *,
								   G4double tmin,
								   G4double maxEnergy);

	G4double CrossSectionPerVolume(const G4Material *material, const G4ParticleDefinition *p, G4double ekin,
								   G4double emin,
								   G4double emax);

	const G4Isotope *SelectRandomIsotope(const G4DynamicParticle *aDynamicGamma, const G4Material *aMaterial);

	int lastM;
	int lastE;

	void ConstructTable();
	void CalculateIntegrals();

	vector<vector<vector<double>>> ElasticData;
	vector<vector<vector<double>>> IntegralData;

	Amplitude Thomson(double Z, double A, double E, double thet);
	Amplitude Delbruck(double Z, double A, double E, double thet);
	Amplitude Rayleigh(double Z, double A, double E, double thet);
	Amplitude GDR(double Z, double A, double E, double thet);

	double InchGDR(double Z, double A, double E, double thet);

	double ComputeCrossS(double Z, double A, double E, double thet);
	double GetRandom();
	double Interpolate(double x1, double y1, double x2, double y2, double value);

	vector<double> *GetData(const G4Isotope *isotope, double E);

	vector<double> *GetIntegrals();

	double Integrate(vector<double> *array);

	void ReadData();

	vector<vector<double>> DelAmp;
	vector<vector<double>> GDRPar;

	vector<vector<double>> xfact;	 // vector that contains the momentum transfer
	vector<vector<double>> formfact; // vector that contains the form factors

	vector<vector<vector<Amplitude>>> SMAmp;

	vector<double> SMEnergy;
	vector<double> SMTheta;

	vector<vector<double>> MFormFactors;
	vector<vector<double>> MXValues;

	vector<vector<double>> ANMRealFactors;
	vector<vector<double>> ANMRealEValues;

	vector<vector<double>> ANMImagFactors;
	vector<vector<double>> ANMImagEValues;

	vector<double> CCTheta;
	vector<vector<double>> CCAmp;

protected:
	G4ParticleDefinition *theGamma;
	G4ParticleDefinition *theElectron;
	G4ParticleChangeForGamma *fParticleChange;

private:
	vector<G4double> test;

	G4Elastic &operator=(const G4Elastic &right);
	G4Elastic(const G4Elastic &);

	vector<G4Isotope *> IsotopeVector;

	int GranP;
	int GranT;
	int GranE;
	double MaxEnergy;
	double MinEnergy;
	double EnInterval;
	std::vector<G4double> xsec;

	int stokesXi1 = 0;
};
