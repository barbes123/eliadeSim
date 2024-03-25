#include "AppsPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iostream>
#include "AppsInput.hh"
#include <random>
#include "G4RandomDirection.hh"
#include <chrono>
#include "AppsAnalysisManager.hh"
#include <atomic>
#include "AppsRun.hh"

G4Mutex myLowEPrimGenMutex = G4MUTEX_INITIALIZER;

using namespace std;
/*
static std::atomic<long long> Macropulse(0);
static std::atomic<long long> Micropulse(0);
static std::atomic<long long> Time(0);
*/

long long AppsPrimaryGeneratorAction::Macropulse = 0;
long long AppsPrimaryGeneratorAction::Micropulse = 0;
long long AppsPrimaryGeneratorAction::Time = 0;

AppsPrimaryGeneratorAction::AppsPrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	  fParticleGun(0)
{

	AppsInput *FInput = AppsInput::Instance();

	fParticleGun = new G4ParticleGun(1);

	// default particle kinematic
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

	// defined type of particle

	G4ParticleDefinition *particle = particleTable->FindParticle("gamma");
	fParticleGun->SetParticleDefinition(particle);

	// set source position (point like, fixed position)
	sourcePos = FInput->GetSourcePosition();

	// retrive beam time structure type
	BeamType = FInput->GetBeamType();

	// retrive and set beam polarization
	beamPolarization = FInput->GetPolarization();
	fParticleGun->SetParticlePolarization(beamPolarization);

	if (BeamType == "speck")
	{

		// retrive angular data for the source

		PolarAngle = FInput->GetPolarAngle();
		AzimuthalAngle = FInput->GetAzimuthalAngle();

		// retrive the beam time structure
		TimeStructure = FInput->GetTimeStructureParameters();

		// retrive the number of beam particles
		BunchParticles = FInput->GetBunchParticles();
		MacroBunchParticles = BunchParticles * TimeStructure[0]; // number of macro bunch particles

		// retrive the beam energy vector, currently implemented for single energy
		vector<G4double> GammaEnergyVector = FInput->GetGammaEnergy();
		GammaEnergy = GammaEnergyVector[0] * MeV;

		// retrive the beam bandwidth vector, currently implemented for single bandwidth
		vector<G4double> GammaBandwidthVector = FInput->GetGammaBandwidth();
		GammaBandwidth = GammaBandwidthVector[0] * MeV;

		// retrive shape parameters
		shapeParameters = FInput->GetShapeParameters();

		std::random_device rd;	// Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(-shapeParameters[0], shapeParameters[0]);

		SurfR1 = distribution(gen);
		SurfR2 = distribution(gen);
	}

	if (BeamType == "surf")
	{

		// retrive angular data for the source

		PolarAngle = FInput->GetPolarAngle();
		AzimuthalAngle = FInput->GetAzimuthalAngle();

		// retrive the beam time structure
		TimeStructure = FInput->GetTimeStructureParameters();

		// retrive the number of beam particles
		BunchParticles = FInput->GetBunchParticles();
		MacroBunchParticles = BunchParticles * TimeStructure[0]; // number of macro bunch particles

		// retrive the beam energy vector, currently implemented for single energy
		vector<G4double> GammaEnergyVector = FInput->GetGammaEnergy();
		GammaEnergy = GammaEnergyVector[0] * MeV;

		// retrive the beam bandwidth vector, currently implemented for single bandwidth
		vector<G4double> GammaBandwidthVector = FInput->GetGammaBandwidth();
		GammaBandwidth = GammaBandwidthVector[0] * MeV;

		// retrive shape parameters
		shapeParameters = FInput->GetShapeParameters();
	}

	if (BeamType == "GBS")
	{

		// retrive angular data for the source

		PolarAngle = FInput->GetPolarAngle();
		AzimuthalAngle = FInput->GetAzimuthalAngle();

		// retrive the beam time structure
		TimeStructure = FInput->GetTimeStructureParameters();

		// retrive the number of beam particles
		BunchParticles = FInput->GetBunchParticles();
		MacroBunchParticles = BunchParticles * TimeStructure[0]; // number of macro bunch particles

		// retrive the beam energy vector, currently implemented for single energy
		vector<G4double> GammaEnergyVector = FInput->GetGammaEnergy();
		GammaEnergy = GammaEnergyVector[0] * MeV;

		// retrive the beam bandwidth vector, currently implemented for single bandwidth
		vector<G4double> GammaBandwidthVector = FInput->GetGammaBandwidth();
		GammaBandwidth = GammaBandwidthVector[0] * MeV;
	}

	if (BeamType == "VEGA")
	{
		// not sure about the random engine
		vegaSource = new GammaGenerator(new CLHEP::MixMaxRng);

		// retrive the beam time structure
		TimeStructure = FInput->GetTimeStructureParameters();

		// retrive the number of beam particles
		BunchParticles = FInput->GetBunchParticles();
		MacroBunchParticles = BunchParticles * TimeStructure[0]; // number of macro bunch particles
	}

	if (BeamType == "source")
	{

		// schedules the number of events to be processed for one real second source indent

		// retrieve the decay nucleus
		DecayNucleus = FInput->GetDecayNucleus();

		for (int i = 0; i < 1000; i++)
		{

			G4double angle = i * (CLHEP::pi / 1000);
			double value = (1 + (1.0 / 8) * pow(cos(angle), 2) + (1.0 / 24) * pow(cos(angle), 4)) * sin(angle);
			AngularDistrCo.push_back(value);
		}

		sourceActivity = FInput->GetActivityParameters();

		if (Macropulse == 0)
		{
			Macropulse = G4RandGauss::shoot(sourceActivity[0] / 100, sourceActivity[1] / 100);
		}
	}

	if (BeamType == "PaNi")
	{

		// schedules the number of events to be processed for one real second source indent


		for (int i = 0; i < 1000; i++)
		{

			G4double angle = i * (CLHEP::pi / 1000);
			double value = (1 + (1.0 / 8) * pow(cos(angle), 2) + (1.0 / 24) * pow(cos(angle), 4)) * sin(angle);
			AngularDistrCo.push_back(value);
		}
		
		FInput->ClearActivityParameters();
		FInput->AddActivityParameter(222000.);
		FInput->AddActivityParameter(500.);

		sourceActivity = FInput->GetActivityParameters();

		if (Macropulse == 0)
		{
			Macropulse = G4RandGauss::shoot(sourceActivity[0] / 100, sourceActivity[1] / 100);
		}
	}

	if (BeamType == "simple")
	{

		// retrive angular data for the source

		PolarAngle = FInput->GetPolarAngle();
		AzimuthalAngle = FInput->GetAzimuthalAngle();

		// retrive the beam energy vector, currently implemented for single energy
		vector<G4double> GammaEnergyVector = FInput->GetGammaEnergy();
		GammaEnergy = GammaEnergyVector[0] * MeV;

		// retrive the beam bandwidth vector, currently implemented for single bandwidth
		vector<G4double> GammaBandwidthVector = FInput->GetGammaBandwidth();
		GammaBandwidth = GammaBandwidthVector[0] * MeV;
	}

	if (BeamType == "multiple")
	{

		// retrive angular data for the source

		PolarAngle = FInput->GetPolarAngle();
		AzimuthalAngle = FInput->GetAzimuthalAngle();

		// retrive the number of beam particles
		BunchParticles = FInput->GetBunchParticles();

		// retrive the beam energy vector, currently implemented for single energy
		vector<G4double> GammaEnergyVector = FInput->GetGammaEnergy();
		GammaEnergy = GammaEnergyVector[0] * MeV;

		// retrive the beam bandwidth vector, currently implemented for single bandwidth
		vector<G4double> GammaBandwidthVector = FInput->GetGammaBandwidth();
		GammaBandwidth = GammaBandwidthVector[0] * MeV;
	}

	if (BeamType == "function")
	{

		PolarAngle = FInput->GetPolarAngle();
		AzimuthalAngle = FInput->GetAzimuthalAngle();

		shapeParameters = FInput->GetShapeParameters();

		if (shapeParameters[0] == 0)
		{

			functionDistribution = new TF1(
				"func", [&](double *x, double *p)
				{ return p[0] * pow(log(x[0]), -1) +
						 p[1] * pow(log(x[0]), -2) + p[2] * pow(log(x[0]), -3) + p[3] * pow(log(x[0]), -4); },
				4.2, 7.5, 4);
			functionDistribution->FixParameter(0, -409.063);
			functionDistribution->FixParameter(1, 1394.98);
			functionDistribution->FixParameter(2, -1178.07);
			functionDistribution->FixParameter(3, 330.936);
		}

		if (shapeParameters[0] == 10)
		{

			functionDistribution = new TF1(
				"func", [&](double *x, double *p)
				{ return pow(x[0] / p[0] - p[1], 2) + p[2]; },
				4.2, 7.5, 3);
			functionDistribution->SetParameter(0, 2);
			functionDistribution->SetParameter(1, 2.1);
			functionDistribution->SetParameter(2, 1);
		}

		if (shapeParameters[0] == 20)
		{

			// retrive the beam energy vector, currently implemented for single energy
			vector<G4double> GammaEnergyVector = FInput->GetGammaEnergy();
			GammaEnergy = GammaEnergyVector[0] * MeV;

			string fileName = "9MeV.dat";

			if (GammaEnergy == 6)
			{
				fileName = "6MeV.dat";
			}

			ifstream myfile(fileName.c_str());

			vector<int> distribution;
			vector<double> energy;

			double val = 20 / 16000.0;

			int cnt = 1;
			while (!myfile.eof())
			{
				double x;

				myfile >> x;

				// distribution
				distribution.push_back(x);

				// energy
				energy.push_back((cnt - 1) * val);

				cnt++;
			}
			energy.push_back(cnt * val);

			std::random_device rd;
			bremGenerator = mt19937(rd());
			bremDistribution = piecewise_constant_distribution<>(energy.begin(), energy.end(), distribution.begin());
		}
	}
}

AppsPrimaryGeneratorAction::~AppsPrimaryGeneratorAction()
{

	delete fParticleGun;
}

void AppsPrimaryGeneratorAction::SetDecaySchedule()
{

	/*

	// appends to DecaySchedule the number of events to be processed for one real second

	AppsInput* FInput = AppsInput::Instance();

	vector<G4double> sourceActivity = FInput->GetActivityParameters();
	G4String decayNucleus = FInput->GetDecayNucleus();
	G4double totalEvents = FInput->GetNumberOfEvents()[0];					// number of events for single source type

	G4double sumEvents = 0;

	while(sumEvents < totalEvents){

		G4int nEvents = G4RandGauss::shoot(sourceActivity[0],sourceActivity[1]);

		G4cout << nEvents << endl;

		sumEvents += nEvents; 							// per 10 ms


		if(sumEvents < totalEvents){
			DecaySchedule.push_back(nEvents);
		} else {
			DecaySchedule.push_back(nEvents-(sumEvents-totalEvents));
		}
	}
*/
}

void AppsPrimaryGeneratorAction::GenerateFunctionEvent(G4Event *anEvent)
{

	// set particle energy
	fParticleGun->SetParticleEnergy(GetParticleEnergy());

	// set particle momentum direction
	fParticleGun->SetParticleMomentumDirection(GetParticleDirection());

	// set particle position
	fParticleGun->SetParticlePosition(sourcePos);

	// generate particle
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void AppsPrimaryGeneratorAction::GenerateSimpleEvent(G4Event *anEvent)
{

	// set particle time
	fParticleGun->SetParticleTime(GetParticleTime());

	// set particle energy

	fParticleGun->SetParticleEnergy(GetParticleEnergy());

	// set particle momentum direction

	fParticleGun->SetParticleMomentumDirection(GetParticleDirection());

	// set particle position

	fParticleGun->SetParticlePosition(sourcePos);

	// generate particle
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void AppsPrimaryGeneratorAction::GenerateMultipleEvent(G4Event *anEvent)
{

	// set particle time
	fParticleGun->SetParticleTime(GetParticleTime());

	// set particle energy

	fParticleGun->SetParticleEnergy(GetParticleEnergy());

	// set particle position

	fParticleGun->SetParticlePosition(sourcePos);

	for (int i = 0; i < BunchParticles; i++){
		// set particle momentum direction

		fParticleGun->SetParticleMomentumDirection(GetParticleDirection());

		// generate particle
		fParticleGun->GeneratePrimaryVertex(anEvent);
	}
}

void AppsPrimaryGeneratorAction::GenerateGBSEvent(G4Event *anEvent)
{

	// set particle time
	fParticleGun->SetParticleTime(GetParticleTime());

	// set particle energy

	fParticleGun->SetParticleEnergy(GetParticleEnergy());

	// set particle momentum direction

	fParticleGun->SetParticleMomentumDirection(GetParticleDirection());

	fParticleGun->SetParticlePosition(sourcePos);

	// generate particle

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void AppsPrimaryGeneratorAction::GenerateVEGAEvent(G4Event *anEvent)
{

	MCGamma *gamma = vegaSource->Gamma();

	fParticleGun->SetParticleEnergy(gamma->energy);
	fParticleGun->SetParticleTime(GetParticleTime());
	fParticleGun->SetParticlePosition(gamma->position + sourcePos);
	fParticleGun->SetParticlePolarization(gamma->polarization);
	fParticleGun->SetParticleMomentumDirection(gamma->direction);

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void AppsPrimaryGeneratorAction::GenerateSurfEvent(G4Event *anEvent)
{

	// set particle time

	fParticleGun->SetParticleTime(GetParticleTime()); // checked

	// set particle energy

	fParticleGun->SetParticleEnergy(GetParticleEnergy()); // checked

	// set particle momentum direction

	fParticleGun->SetParticleMomentumDirection(GetParticleDirection()); // checked

	// set particle position

	fParticleGun->SetParticlePosition(GetParticlePosition());

	// generate particle

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void AppsPrimaryGeneratorAction::GenerateSpeckEvent(G4Event *anEvent)
{

	// set particle time

	fParticleGun->SetParticleTime(GetParticleTime()); // checked

	// set particle energy

	fParticleGun->SetParticleEnergy(GetParticleEnergy()); // checked

	// set particle momentum direction

	fParticleGun->SetParticleMomentumDirection(GetParticleDirection()); // checked

	// set particle position

	fParticleGun->SetParticlePosition(GetParticlePosition());

	// generate particle

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void AppsPrimaryGeneratorAction::GeneratePuBeEvent(G4Event *anEvent)
{

	// set particle definition

	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *particle = particleTable->FindParticle("neutron");
	fParticleGun->SetParticleDefinition(particle);

	// set particle time

	fParticleGun->SetParticleTime(GetParticleTime());

	// set particle energy

	fParticleGun->SetParticleEnergy(GetParticleEnergy());

	// set particle momentum direction

	fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

	// set particle position

	fParticleGun->SetParticlePosition(sourcePos);

	// generate particle

	fParticleGun->GeneratePrimaryVertex(anEvent);
		
}

void AppsPrimaryGeneratorAction::GenerateSourceEvent(G4Event *anEvent)
{

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	//G4int histID = analysisManager->GetH1Id("Beam_Edist");

	// set particle position

	if (DecayNucleus == "Co60")
	{

		// set particle time

		fParticleGun->SetParticleTime(GetParticleTime());

		// set particle energy

		G4double energy = 1.332 * MeV;

		fParticleGun->SetParticleEnergy(energy);

		//analysisManager->FillH1(histID, energy);

		// set particle momentum direction

		fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

		// set particle position

		fParticleGun->SetParticlePosition(sourcePos);

		// generate particle

		fParticleGun->GeneratePrimaryVertex(anEvent);

		G4double rand = G4UniformRand();

		if (rand < 0.9985)
		{ // condition to emit the second gamma
			// particle will be emitted at the same time as the first particle
			G4ThreeVector ParticleDirection;
			// set particle energy

			energy = 1.173 * MeV;

			fParticleGun->SetParticleEnergy(energy);

			//analysisManager->FillH1(histID, energy);

			// set particle momentum direction

			G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();

			G4double Theta = GetRandom(AngularDistrCo);

			G4double rndm = G4UniformRand();

			G4double costheta = cos(Theta);
			G4double sintheta = sin(Theta);

			G4double MinPhi = 0;
			G4double MaxPhi = CLHEP::twopi;

			rndm = G4UniformRand();
			G4double Phi = MinPhi + (MaxPhi - MinPhi) * rndm;

			G4double sinphi = std::sin(Phi);
			G4double cosphi = std::cos(Phi);

			G4double px = sintheta * cosphi;
			G4double py = sintheta * sinphi;
			G4double pz = costheta;

			ParticleDirection.setX(px);
			ParticleDirection.setY(py);
			ParticleDirection.setZ(pz);

			dir.rotateUz(ParticleDirection);

			// fParticleGun->SetParticleMomentumDirection(dir);
			fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

			fParticleGun->SetParticlePosition(sourcePos);

			fParticleGun->GeneratePrimaryVertex(anEvent);
		}
	}
}

G4double AppsPrimaryGeneratorAction::GetRandom(vector<double> array)
{

	G4int nbins = array.size();

	if (AngularIntCo.size() == 0)
	{
		AngularIntCo.reserve(nbins + 2);
		G4int ibin = 0;
		AngularIntCo.push_back(0);

		for (G4int i = 1; i <= nbins; i++)
		{
			++ibin;
			AngularIntCo.push_back(AngularIntCo.back() + (array)[i - 1]);
		}

		for (G4int bin = 0; bin <= nbins; bin++)
		{
			AngularIntCo[bin] /= AngularIntCo[nbins];
		}
	}

	G4double r1 = G4UniformRand();

	std::vector<G4double>::iterator pind;

	pind = std::lower_bound((AngularIntCo).begin(), (AngularIntCo).end(), r1);

	int newBin;

	if ((pind != (AngularIntCo).end()) && (*pind == r1))
	{
		newBin = std::distance((AngularIntCo).begin(), pind);
	}
	else
	{
		newBin = std::distance((AngularIntCo).begin(), pind - 1);
	}

	G4double BinWidth = CLHEP::pi / 1000;

	G4double x = newBin * BinWidth;

	if (r1 > (AngularIntCo)[newBin])
	{
		x += 1.0 * BinWidth * (r1 - (AngularIntCo)[newBin]) / ((AngularIntCo)[newBin + 1] - (AngularIntCo)[newBin]);
	}

	return x;
}

void AppsPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{

	if (BeamType == "GBS")
	{
		GenerateGBSEvent(anEvent);
	}

	if (BeamType == "source")
	{
		GenerateSourceEvent(anEvent);
	}

	if (BeamType == "PaNi")
	{
		GeneratePuBeEvent(anEvent);
	}

	if (BeamType == "simple")
	{
		GenerateSimpleEvent(anEvent);
	}

	if (BeamType == "multiple")
	{
		GenerateMultipleEvent(anEvent);
	}

	if (BeamType == "surf")
	{
		GenerateSurfEvent(anEvent);
	}

	if (BeamType == "speck")
	{
		GenerateSpeckEvent(anEvent);
	}

	if (BeamType == "function")
	{
		GenerateFunctionEvent(anEvent);
	}

	if (BeamType == "VEGA")
	{
		GenerateVEGAEvent(anEvent);
	}
}

G4ThreeVector AppsPrimaryGeneratorAction::GetParticleDirection()
{

	G4ThreeVector ParticleDirection;

	if (BeamType == "GBS" || BeamType == "simple" || BeamType == "multiple" || BeamType == "surf" || BeamType == "speck" || BeamType == "function")
	{ // particle direction implemented for beam parameters

		G4double MinTheta = PolarAngle[0];
		G4double MaxTheta = PolarAngle[1];
		G4double MaxPhi = AzimuthalAngle[0];
		G4double MinPhi = AzimuthalAngle[1];

		G4double rndm = G4UniformRand();

		G4double costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta) - std::cos(MaxTheta));

		G4double sintheta = std::sqrt(1. - costheta * costheta);

		rndm = G4UniformRand();
		G4double Phi = MinPhi + (MaxPhi - MinPhi) * rndm;

		G4double sinphi = std::sin(Phi);
		G4double cosphi = std::cos(Phi);

		G4double px = sintheta * cosphi;
		G4double py = sintheta * sinphi;
		G4double pz = costheta;

		ParticleDirection.setX(px);
		ParticleDirection.setY(py);
		ParticleDirection.setZ(pz);

		return ParticleDirection;
	}

	if ((BeamType == "source")||(BeamType == "PaNi"))
	{

		return G4RandomDirection();
	}

	return ParticleDirection;
}

G4double AppsPrimaryGeneratorAction::GetParticleEnergy()
{

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	AppsInput *FInput = AppsInput::Instance();

	if (BeamType == "GBS" || BeamType == "simple" || BeamType == "multiple" || BeamType == "surf" || BeamType == "speck")
	{

		partEnergy = G4RandGauss::shoot(GammaEnergy, GammaBandwidth);
	}

	if (BeamType == "function")
	{

		// partEnergy = 1*MeV;
		if (shapeParameters[0] == 10 || shapeParameters[0] == 0)
		{
			partEnergy = functionDistribution->GetRandom();
		}

		if (shapeParameters[0] == 20)
		{
			partEnergy = bremDistribution(bremGenerator);
		}
	}

	if (BeamType == "PaNi")
	{

		partEnergy = FInput->GetNeutronSpectrum()->GetRandom()/1000.;
		//FInput->GetNeutronSpectrum()->Print("all");

	}

	//G4int histID = analysisManager->GetH1Id("Beam_Edist");

	//analysisManager->FillH1(histID, partEnergy);

	return partEnergy;
}

G4double AppsPrimaryGeneratorAction::GetParticleTime()
{

	G4AutoLock lock(&myLowEPrimGenMutex);

	AppsInput *FInput = AppsInput::Instance();

	if (BeamType == "GBS" || BeamType == "surf" || BeamType == "speck" || BeamType == "VEGA")
	{ // condition that returns the time of the particles in case of GBS time structure
		// with hardcoded photons per bunch
		Macropulse++;
		Micropulse++;

		if (Micropulse == BunchParticles)
		{

			Micropulse = 0;

			Time = Time + TimeStructure[1];
		}

		if (Macropulse == MacroBunchParticles)
		{

			Macropulse = 0;

			Time = Time + TimeStructure[2];

			//	G4cout << " Macropulse time " << FInput->Time/CLHEP::ns << " Macropulse " << FInput->Macropulse << endl;
		}

		return Time;
	}

	if ((BeamType == "source")||(BeamType == "PaNi"))
	{

		// the macropulse is considered a second in real time and counts the number of particles emitted in the current second
		// the time counts the total emitted second
		Micropulse++;

		if (Micropulse > Macropulse)
		{

			Micropulse = 0;
			Time = Time + 10 * CLHEP::ms;
			Macropulse = G4RandGauss::shoot(sourceActivity[0] / 100, sourceActivity[1] / 100);
		}

		G4double start = Time;
		G4double stop = Time + 10 * CLHEP::ms;

		Time = G4RandFlat::shoot(start, stop);

		return Time;
	}

	if (BeamType == "simple" || BeamType == "multiple")
	{

		// one event every 10 ms
		Micropulse++;
		Macropulse++;

		Time = Time + 10 * CLHEP::ms;

		return Time;
	}

	return Time;
}

G4ThreeVector AppsPrimaryGeneratorAction::GetParticlePosition()
{

	G4ThreeVector position;

	G4double x = 0;
	G4double y = 0;

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	if (BeamType == "surf")
	{

		x = G4RandFlat::shoot(-shapeParameters[0], shapeParameters[0]);
		y = G4RandFlat::shoot(-shapeParameters[1], shapeParameters[1]);

		position.setX(x);
		position.setY(y);
		position += sourcePos;

		//G4int histID = analysisManager->GetH2Id("Beam_PosDist");
		//analysisManager->FillH2(histID, x / mm, y / mm);

		return position;
	}

	if (BeamType == "speck")
	{

		const AppsRun *localRun = static_cast<const AppsRun *>(G4RunManager::GetRunManager()->GetCurrentRun());

		if (localRun->GetRunID() > currentRun)
		{

			totalSurf = G4RandFlat::shoot(50, 100);
			totalSurf = 1000;
			refSurf = 5;
			currentRun = localRun->GetRunID();
		}

		if (refSurf > 0)
		{

			if (totalSurf > 0)
			{

				// the square

				// x = G4RandFlat::shoot(SurfR1-1,SurfR1+1);
				// y = G4RandFlat::shoot(SurfR2-1,SurfR2+1);

				// the gaus

				x = G4RandGauss::shoot(SurfR1, 0.7);
				y = G4RandGauss::shoot(SurfR2, 0.7);

				//	cout << x << " " << y << endl;

				position.setX(x);
				position.setY(y);

				sourcePos += position;

				totalSurf--;
			}
			else
			{

				std::random_device rd;	// Will be used to obtain a seed for the random number engine
				std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

				std::uniform_real_distribution<double> distribution(-shapeParameters[0], shapeParameters[0]);

				SurfR1 = distribution(gen);
				SurfR2 = distribution(gen);

				x = G4RandGauss::shoot(SurfR1, 0.7);
				y = G4RandGauss::shoot(SurfR2, 0.7);

				position.setX(x);
				position.setY(y);

				sourcePos += position;

				totalSurf--;

				totalSurf = G4RandFlat::shoot(50, 100);

				totalSurf = 1000;
				refSurf--;
			}
		}
		else
		{

			x = G4RandFlat::shoot(-shapeParameters[0], shapeParameters[0]);
			y = G4RandFlat::shoot(-shapeParameters[1], shapeParameters[1]);

			position.setX(x);
			position.setY(y);

			sourcePos += position;
		};

		//G4int histID = analysisManager->GetH2Id("Beam_PosDist");
		//analysisManager->FillH2(histID, x / mm, y / mm);

		return position;
	}

	return position;
}
