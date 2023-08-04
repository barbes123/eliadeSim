#include <iostream>
#include "G4Elastic.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include <math.h>
#include "TF1.h"
#include "TApplication.h"
#include "TH1.h"
#include <complex>
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <sstream>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <random>
#include "TPad.h"
#include "AppsInput.hh"

typedef std::chrono::high_resolution_clock Clock;

struct Amplitude
{
	double AppReal = 0;
	double AppImag = 0;
	double AparReal = 0;
	double AparImag = 0;
};

G4double c = 2.99792458 * pow(10, 8);					  // speed of light
G4double c2 = c * c;									  // speed of light squared
G4double Amass = 1.660538921 * pow(10, -27);			  // atomic mass
G4double e = 1.602176 * pow(10, -19);					  // elementary charge
G4double e2 = e * e;									  // elementary charge squared
G4double eps0 = 8.854187 * pow(10, -12);				  // vacuum permeability
G4double hbar = 6.5821195 * pow(10, -22);				  // Plank constant
G4double alph = 1 / 137.035999139;						  // solid constant
G4double Elrad = 0.00000000000000281794;				  // classical electron radius
G4double Me = 9.109383 * pow(10, -31);					  // electron mass
G4double Mn = 1.672621 * pow(10, -27);					  // nucleon mass
G4double hMeV = 4.135667662 * pow(10, -15) * pow(10, -6); // Plank constant in MeV
G4double Elrad2 = 0.000000000000281794;					  // classical electron radius

G4Elastic::G4Elastic(const G4ParticleDefinition *, const G4String &nam)
	: G4VEmModel(nam)
{

	theGamma = G4Gamma::Gamma();
	fParticleChange = 0;
}

G4Elastic::~G4Elastic()
{
}

void G4Elastic::Initialise(const G4ParticleDefinition *,
						   const G4DataVector &)
{

	if (!fParticleChange)
	{
		fParticleChange = GetParticleChangeForGamma();
	}

	G4ProductionCutsTable *theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();

	GranP = 361;						// user allocated granularity for phi
	GranT = 181;						// user allocated granularity for theta
	GranE = 1000;						// user allocated granularity for energy
	MaxEnergy = 21;						// maximum energy presented in the simulation can be extracted from the gun max energy
	MinEnergy = 0.1;					// min energy for the validity of the model
	EnInterval = MaxEnergy - MinEnergy; // energy interval between the two

	auto FInput = AppsInput::Instance();

	stokesXi1 = FInput->GetPolarization().getX();

	for (size_t i = 0; i < theCoupleTable->GetTableSize(); i++)
	{

		const G4Material *material = theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
		const G4ElementVector *theElementVector = material->GetElementVector();

		for (size_t j = 0; j < material->GetNumberOfElements(); j++)
		{

			const G4IsotopeVector *theIsotopeVector = theElementVector->at(j)->GetIsotopeVector();

			for (unsigned int l = 0; l < theIsotopeVector->size(); l++)
			{

				int theZ = theIsotopeVector->at(l)->GetZ();
				int theA = theIsotopeVector->at(l)->GetN();

				if (IsotopeVector.size() == 0)
				{
					IsotopeVector.push_back(theIsotopeVector->at(l));
				}

				int check = 0;

				for (unsigned int k = 0; k < IsotopeVector.size(); k++)
				{

					if (IsotopeVector[k]->GetZ() == theZ && IsotopeVector[k]->GetN() == theA)
					{
						check = 1;
					}
				}

				if (check == 0)
				{
					IsotopeVector.push_back(theIsotopeVector->at(l));
				}
			}
		}
	}

	ConstructTable();
}

G4double G4Elastic::CrossSectionPerVolume(const G4Material *material,
										  const G4ParticleDefinition *,
										  G4double ekin,
										  G4double,
										  G4double)
{

	if (ekin > 20 * MeV)
	{
		return 0;
	}

	G4double cross = 0.0;

	const G4ElementVector *theElementVector = material->GetElementVector();
	const G4double *theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
	G4int nelm = material->GetNumberOfElements();

	xsec.clear();

	for (int i = 0; i < nelm; i++)
	{

		const G4IsotopeVector *theIsotopeVector = (*theElementVector)[i]->GetIsotopeVector();
		G4int niso = (*theElementVector)[i]->GetNumberOfIsotopes();
		const G4double *FractionVector = (*theElementVector)[i]->GetRelativeAbundanceVector();

		for (int o = 0; o < niso; o++)
		{
			cross += theAtomNumDensityVector[i] * ComputeCrossSectionPerAtom(ekin, (*theIsotopeVector)[o]) * FractionVector[o];
			xsec.push_back(cross);
		}
	}

	return cross;
}

G4double G4Elastic::ComputeCrossSectionPerAtom(double ekin, G4Isotope *isotope)
{

	G4double gammaEnergy = ekin;

	double crossSec = 0;

	TGraph *gr1 = new TGraph();

	vector<double> *a = GetData(isotope, gammaEnergy);

	for (unsigned int i = 0; i < a->size(); i++)
	{
		gr1->SetPoint(i + 1, (i * pi) / 180, (*a)[i]);
	}

	crossSec = gr1->Integral();

	delete gr1;

	return crossSec * meter2;
}

const G4Isotope *G4Elastic::SelectRandomIsotope(const G4DynamicParticle *aDynamicGamma, const G4Material *aMaterial)
{

	const G4int numElem = aMaterial->GetNumberOfElements();
	const G4ElementVector *theElementVector = aMaterial->GetElementVector();

	if (numElem == 1)
	{

		G4int numIso = (*theElementVector)[0]->GetNumberOfIsotopes();

		if (numIso == 1)
		{
			return (*theElementVector)[0]->GetIsotope(0);
		}
		else
		{

			G4double x = G4UniformRand() * CrossSectionPerVolume(aMaterial, aDynamicGamma->GetParticleDefinition(), aDynamicGamma->GetKineticEnergy(), 0, 0);
			for (G4int i = 0; i < numIso; i++)
			{
				if (x <= xsec[i])
				{
					return (*theElementVector)[0]->GetIsotope(i);
				}
			}
		}
	}
	else
	{

		G4double x = G4UniformRand() * CrossSectionPerVolume(aMaterial, aDynamicGamma->GetParticleDefinition(), aDynamicGamma->GetKineticEnergy(), 0, 0);
		G4int cnt = 0;

		for (G4int i = 0; i < numElem; i++)
		{

			G4int numIso = (*theElementVector)[i]->GetNumberOfIsotopes();

			for (G4int o = 0; o < numIso; o++)
			{
				if (x <= xsec[cnt])
				{
					return (*theElementVector)[i]->GetIsotope(o);
				}
				cnt++;
			}
		}
	}

	return (*theElementVector)[0]->GetIsotope(0);
}

void G4Elastic::SampleSecondaries(std::vector<G4DynamicParticle *> *,
								  const G4MaterialCutsCouple *couple,
								  const G4DynamicParticle *aDynamicGamma,
								  G4double,
								  G4double)
{

	G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();

	if (GammaEnergy < MinEnergy)
	{
		return;
	}

	const G4Isotope *isotope = SelectRandomIsotope(aDynamicGamma, couple->GetMaterial());

	auto crossData = GetData(isotope, GammaEnergy);

	double thet = GetRandom() * pi / 180;

	double crossValue = crossData->at(int(thet * 180 / pi));

	crossValue = crossValue / (0.25 * twopi * sin(thet));

	vector<double> phiIndex(GranP);
	iota(phiIndex.begin(), phiIndex.end(), 0);

	vector<double> phiCross = phiIndex;

	auto xi = stokesXi1;

	transform(phiCross.cbegin(), phiCross.cend(), phiCross.begin(), [crossValue, xi](int i)
			  { return 0.25 * (crossValue + xi * crossValue * cos(2 * i * pi / 180)); });

	// int argc = 0;

	// TApplication test("theApp", &argc, {});

	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine gen(seed1);

	piecewise_linear_distribution<> d{phiIndex.begin(), phiIndex.end(), phiCross.begin()};

	// auto hist = new TH1D("hist", "hist", 360, 0, 360);

	// for (int i = 0; i < 1000000; i++)
	// {
	// 	hist->Fill(d(gen));
	// }

	// auto gr1 = new TGraph();

	// for (int i = 0; i < phiIndex.size(); i++)
	// {
	// 	gr1->SetPoint(i + 1, phiIndex[i], phiCross[i]);
	// }

	// hist->Draw();
	// gPad->WaitPrimitive("ggg");

	// gr1->Draw();
	// gPad->WaitPrimitive("ggg");

	double phi = d(gen) * pi / 180;

	double costhet = cos(thet);
	double sinthet = sin(thet);

	G4ThreeVector gamDirection0 = aDynamicGamma->GetMomentumDirection();
	G4ThreeVector gamDirection1(sinthet * cos(phi), sinthet * sin(phi), costhet);
	gamDirection1.rotateUz(gamDirection0);

	fParticleChange->ProposeMomentumDirection(gamDirection1);
}

vector<double> *G4Elastic::GetData(const G4Isotope *isotope, double E)
{

	for (unsigned int i = 0; i < IsotopeVector.size(); i++)
	{
		if (IsotopeVector[i]->GetN() == isotope->GetN() && IsotopeVector[i]->GetZ() == isotope->GetZ())
		{

			lastM = i;
			lastE = (int)((E - MinEnergy) * GranE * 1.0 / (EnInterval));

			if (ElasticData[lastM][lastE].size() == 0)
			{

				// determine the location of the new energy data
				for (int o = 0; o < GranE; o++)
				{
					if (E >= MinEnergy + o * 1.0 * (EnInterval) / GranE && E < MinEnergy + (o + 1) * 1.0 * (EnInterval) / GranE)
					{
						E = MinEnergy + o * 1.0 * (EnInterval) / GranE;
					}
				}

				// calculate cross section and fill the table
				vector<double> Ttemp;
				for (int u = 0; u < GranT; u++)
				{
					double thet = u * 1.0 * pi / GranT;
					double CrossSec = ComputeCrossS(IsotopeVector[i]->GetZ(), IsotopeVector[i]->GetN(), E, thet);
					Ttemp.push_back(CrossSec);
				}

				ElasticData[lastM][lastE] = Ttemp;

				// calculate the integral
				vector<double> IntTemp;
				IntTemp.reserve(GranT + 2);

				IntTemp.push_back(0);

				for (int u = 1; u <= GranT; u++)
				{
					IntTemp.push_back(IntTemp.back() + ElasticData[lastM][lastE][u - 1]);
				}

				for (int u = 0; u <= GranT; u++)
				{
					IntTemp[u] /= IntTemp[GranT];
				}

				IntegralData[lastM][lastE] = IntTemp;
			}

			return &ElasticData[lastM][lastE];
		}
	}

	return &ElasticData[0][0];
}

void G4Elastic::ConstructTable()
{

	ReadData();

	ElasticData.reserve(IsotopeVector.size());
	IntegralData.reserve(IsotopeVector.size());

	for (unsigned int i = 0; i < ElasticData.capacity(); i++)
	{
		vector<vector<double>> Etemp;
		vector<vector<double>> Itemp;
		Etemp.reserve(GranE);

		for (unsigned int o = 0; o < Etemp.capacity(); o++)
		{
			vector<double> Ttemp;

			Ttemp.reserve(GranT);
			Etemp.push_back(Ttemp);

			Ttemp.reserve(GranT + 2);
			Itemp.push_back(Ttemp);
		}
		ElasticData.push_back(Etemp);
		IntegralData.push_back(Itemp);
	}
}

vector<double> *G4Elastic::GetIntegrals()
{

	return &IntegralData[lastM][lastE];
}

void G4Elastic::CalculateIntegrals()
{

	IntegralData.reserve(ElasticData.capacity());

	for (unsigned int i = 0; i < IntegralData.capacity(); i++)
	{
		vector<vector<double>> InTemp;
		InTemp.reserve(GranE);

		for (unsigned int o = 0; o < InTemp.capacity(); o++)
		{
			vector<double> IntTemp;
			IntTemp.reserve(GranT + 2);

			IntTemp.push_back(0);
			for (int u = 1; u <= GranT; u++)
			{
				IntTemp.push_back(IntTemp.back() + ElasticData[i][o][u - 1]);
			}

			for (int u = 0; u <= GranT; u++)
			{
				IntTemp[u] /= IntTemp[GranT];
			}

			InTemp.push_back(IntTemp);
		}
		IntegralData.push_back(InTemp);
	}
}

double G4Elastic::ComputeCrossS(double Z, double A, double E, double thet)
{

	Amplitude amp1 = Thomson(Z, A, E, thet);
	Amplitude amp2 = Rayleigh(Z, A, E, thet);
	Amplitude amp3 = Delbruck(Z, A, E, thet);
	Amplitude amp4 = GDR(Z, A, E, thet);

	complex<double> AppSum;
	complex<double> AparSum;

	AppSum.real(amp1.AppReal + amp2.AppReal + amp3.AppReal + amp4.AppReal);
	AppSum.imag(amp1.AppImag + amp2.AppImag + amp3.AppImag + amp4.AppImag);

	AparSum.real(amp1.AparReal + amp2.AparReal + amp3.AparReal + amp4.AparReal);
	AparSum.imag(amp1.AparImag + amp2.AparImag + amp3.AparImag + amp4.AparImag);

	double CrossInchGDR = InchGDR(Z, A, E, thet);

	// double CrossSec = ((1.0 / 2.0) * (pow(abs(AppSum), 2) + pow(abs(AparSum), 2)) + CrossInchGDR) * twopi * sin(thet);
	double CrossSec = (1.0 / 4.0) * (pow(abs(AppSum), 2) + pow(abs(AparSum), 2)) * twopi * sin(thet);

	return CrossSec;
}

void G4Elastic::ReadData()
{

	// Read delbruck //

	vector<double> DelData;

	DelAmp.reserve(33);
	DelData.reserve(76);

	double tmp;

	string filename = "DelbruckData.dat";
	ifstream myfile(filename.c_str(), ios::in);
	if (!myfile)
	{
		cout << " Delbruck data file is missing " << endl;
		abort();
	}

	while (!myfile.eof())
	{
		for (int i = 0; i < 76; i++)
		{
			myfile >> tmp;
			DelData.push_back(tmp);
		}
		DelAmp.push_back(DelData);
		DelData.clear();
	}

	myfile.close();

	// Read GDR //

	filename = "GDRData.dat";

	myfile.open(filename.c_str());
	if (!myfile)
	{
		cout << " File " << filename << " not found " << endl;
		abort();
	}

	string useless;
	double Zf;
	double Af;

	for (unsigned int o = 0; o < IsotopeVector.size(); o++)
	{

		vector<double> GDRtemp;
		GDRtemp.reserve(10);

		for (int i = 0; i < 4; i++)
		{
			getline(myfile, useless);
		}

		while (!myfile.eof())
		{

			myfile >> Zf;
			myfile >> Af;

			if (Zf == IsotopeVector[o]->GetZ() && IsotopeVector[o]->GetN() == Af)
			{

				GDRtemp.push_back(Zf);
				GDRtemp.push_back(Af);

				for (int i = 0; i < 8; i++)
				{

					myfile >> Zf;
					GDRtemp.push_back(Zf);
				}

				GDRPar.push_back(GDRtemp);
				break;
			}

			myfile.ignore(numeric_limits<streamsize>::max(), myfile.widen('\n'));
		}

		myfile.clear();
		myfile.seekg(0, ios::beg);
	}

	myfile.close();

	// Read Rayleigh //

	for (unsigned int o = 0; o < IsotopeVector.size(); o++)
	{

		SMEnergy.clear();
		vector<vector<Amplitude>> tps;
		SMAmp.push_back(tps);
		filename = "_cs0sl_sm";
		ostringstream fname;

		if (IsotopeVector[o]->GetZ() < 10)
		{
			fname << "SMData/"
				  << "00" << IsotopeVector[o]->GetZ() << filename;
		}
		else
		{
			fname << "SMData/"
				  << "0" << IsotopeVector[o]->GetZ() << filename;
		}

		myfile.open(fname.str().c_str());
		if (!myfile)
		{
			cout << " File " << fname.str() << " not found " << endl;
			abort();
		}

		int cnt = 0;
		while (!myfile.eof())
		{

			getline(myfile, useless);

			if (useless.find("REMARK BEGIN") != string::npos)
			{

				useless = useless.substr(0, useless.size() - 1);

				SMTheta.clear();
				vector<Amplitude> temps;
				SMAmp[o].push_back(temps);
				SMEnergy.push_back(atof((useless.substr(0, useless.length() - 3).substr(25)).c_str()) * pow(10, -3));

				for (int i = 0; i < 10; i++)
				{
					getline(myfile, useless);
				}

				istringstream ss;
				Amplitude amp;

				while (useless.find("*** END OF DATA ***") == string::npos)
				{

					ss.str(useless);
					ss >> tmp;
					SMTheta.push_back(tmp * pi / 180);

					ss >> tmp;
					ss >> tmp;
					ss >> tmp;

					amp.AparReal = tmp;

					ss >> tmp;
					amp.AparImag = tmp;

					ss >> tmp;
					amp.AppReal = tmp;

					ss >> tmp;
					amp.AppImag = tmp;

					SMAmp[o][cnt].push_back(amp);
					getline(myfile, useless);

					ss.str("");
				}

				cnt++;
				continue;
			}
		}
		myfile.close();
	}

	// read modified

	for (unsigned int o = 0; o < IsotopeVector.size(); o++)
	{

		vector<double> tmp1;
		vector<double> tmp2;

		ostringstream MFfile;

		if (IsotopeVector[o]->GetZ() < 10)
		{
			MFfile << "MFData/00" << IsotopeVector[o]->GetZ() << "_Mf.out";
		}
		else
		{
			MFfile << "MFData/0" << IsotopeVector[o]->GetZ() << "_Mf.out";
		}

		ifstream FromMf(MFfile.str().c_str(), ios::in);
		if (!FromMf)
		{
			cout << " File " << MFfile.str() << " not found " << endl;
			abort();
		}

		string firstVal;
		string secondVal;

		// modified form factor //

		getline(FromMf, useless);

		while (useless.find("X (1/A)  g(q)/ELECTRON") == string::npos)
		{
			getline(FromMf, useless);
		}

		FromMf >> firstVal >> secondVal;

		while (firstVal != "***")
		{

			tmp1.push_back(atof(firstVal.c_str()));
			tmp2.push_back(atof(secondVal.c_str()) * IsotopeVector[o]->GetZ());
			FromMf >> firstVal >> secondVal;
		}

		MXValues.push_back(tmp1);
		MFormFactors.push_back(tmp2);

		tmp1.clear();
		tmp2.clear();
	}
}

Amplitude G4Elastic::Rayleigh(double Z, double, double E, double thet)
{

	Amplitude Amp;

	Amp.AparReal = 0;
	Amp.AparImag = 0;
	Amp.AppReal = 0;
	Amp.AppImag = 0;

	if (E < 6)
	{

		// SMatrix rayleigh

		for (unsigned int i = 0; i < IsotopeVector.size(); i++)
		{
			if (Z == IsotopeVector[i]->GetZ())
			{
				double valp1;
				double valp2;

				for (unsigned int o = 0; o < (SMAmp[i].size() - 1); o++)
				{
					for (unsigned int p = 0; p < (SMAmp[i][o].size() - 1); p++)
					{

						if (SMEnergy[o] <= E && SMEnergy[o + 1] > E && SMTheta[p] <= thet && SMTheta[p + 1] > thet)
						{

							valp1 = Interpolate(SMEnergy[o], SMAmp[i][o][p].AparReal, SMEnergy[o + 1], SMAmp[i][o + 1][p].AparReal, E);
							valp2 = Interpolate(SMEnergy[o], SMAmp[i][o][p + 1].AparReal, SMEnergy[o + 1], SMAmp[i][o + 1][p + 1].AparReal, E);

							Amp.AparReal = Interpolate(SMTheta[p], valp1, SMTheta[p + 1], valp2, thet) * Elrad;

							valp1 = Interpolate(SMEnergy[o], SMAmp[i][o][p].AparImag, SMEnergy[o + 1], SMAmp[i][o + 1][p].AparImag, E);
							valp2 = Interpolate(SMEnergy[o], SMAmp[i][o][p + 1].AparImag, SMEnergy[o + 1], SMAmp[i][o + 1][p + 1].AparImag, E);

							Amp.AparImag = Interpolate(SMTheta[p], valp1, SMTheta[p + 1], valp2, thet) * Elrad;

							valp1 = Interpolate(SMEnergy[o], SMAmp[i][o][p].AppReal, SMEnergy[o + 1], SMAmp[i][o + 1][p].AppReal, E);
							valp2 = Interpolate(SMEnergy[o], SMAmp[i][o][p + 1].AppReal, SMEnergy[o + 1], SMAmp[i][o + 1][p + 1].AppReal, E);

							Amp.AppReal = Interpolate(SMTheta[p], valp1, SMTheta[p + 1], valp2, thet) * Elrad;

							valp1 = Interpolate(SMEnergy[o], SMAmp[i][o][p].AppImag, SMEnergy[o + 1], SMAmp[i][o + 1][p].AppImag, E);
							valp2 = Interpolate(SMEnergy[o], SMAmp[i][o][p + 1].AppImag, SMEnergy[o + 1], SMAmp[i][o + 1][p + 1].AppImag, E);

							Amp.AppImag = Interpolate(SMTheta[p], valp1, SMTheta[p + 1], valp2, thet) * Elrad;
						}
					}
				}
				return Amp;
			}
		}

		// SMatrix rayleigh
	}
	else
	{

		// Modified form factor  rayleigh

		for (unsigned int i = 0; i < IsotopeVector.size(); i++)
		{
			if (Z == IsotopeVector[i]->GetZ())
			{
				TGraph *data1 = new TGraph();
				double lam = hMeV * c / E;
				double x = sin(thet / 2.0) / (lam * pow(10, 10));

				for (unsigned int o = 0; o < MXValues[i].size(); o++)
				{
					data1->SetPoint(o, MXValues[i][o], MFormFactors[i][o]);
				}

				double ff = data1->Eval(x);

				Amp.AppReal = -Elrad * (ff);
				Amp.AparReal = -Elrad * cos(thet) * (ff);

				delete data1;
				return Amp;
			}
		}
	}
	return Amp;

	// Modified form factor  rayleigh
}

Amplitude G4Elastic::GDR(double Z, double A, double E, double thet)
{

	Amplitude Amp;

	Amp.AparReal = 0;
	Amp.AparImag = 0;
	Amp.AppReal = 0;
	Amp.AppImag = 0;

	double fact1 = 0;
	double fact2 = 0;
	double fact3 = 0;
	double fact4 = 0;
	double fact5 = 0;

	string useless;

	if (GDRPar.size() == 0)
	{
		return Amp;
	}

	for (unsigned int i = 0; i < GDRPar.size(); i++)
	{
		if (GDRPar[i][0] == Z && GDRPar[i][1] == A)
		{

			fact1 = (E * E) / (4 * pi * hbar * c);

			fact2 = GDRPar[i][4] * pow(10, -31) * GDRPar[i][5] * (GDRPar[i][3] * GDRPar[i][3] - E * E) /
					((GDRPar[i][3] * GDRPar[i][3] - E * E) * (GDRPar[i][3] * GDRPar[i][3] - E * E) +
					 E * E * GDRPar[i][5] * GDRPar[i][5]);

			fact3 = GDRPar[i][7] * pow(10, -31) * GDRPar[i][8] * (GDRPar[i][6] * GDRPar[i][6] - E * E) /
					((GDRPar[i][6] * GDRPar[i][6] - E * E) * (GDRPar[i][6] * GDRPar[i][6] - E * E) +
					 E * E * GDRPar[i][8] * GDRPar[i][8]);

			if (E < GDRPar[i][2])
			{
				fact4 = 0;
				fact5 = 0;
			}
			else
			{

				fact4 = GDRPar[i][4] * pow(10, -31) * GDRPar[i][5] * (E * GDRPar[i][5]) / ((GDRPar[i][3] * GDRPar[i][3] - E * E) * (GDRPar[i][3] * GDRPar[i][3] - E * E) + E * E * GDRPar[i][5] * GDRPar[i][5]);

				fact5 = GDRPar[i][7] * pow(10, -31) * GDRPar[i][8] * (E * GDRPar[i][8]) / ((GDRPar[i][6] * GDRPar[i][6] - E * E) * (GDRPar[i][6] * GDRPar[i][6] - E * E) + E * E * GDRPar[i][8] * GDRPar[i][8]);
			}

			Amp.AppReal = fact1 * (fact2 + fact3);
			Amp.AppImag = fact1 * (fact4 + fact5);
			Amp.AparReal = Amp.AppReal * cos(thet);
			Amp.AparImag = Amp.AppImag * cos(thet);

			return Amp;
		}
	}

	return Amp;
}

double G4Elastic::InchGDR(double Z, double A, double E, double thet)
{

	double cross = 0;
	double fact1 = 0;
	double fact2 = 0;
	double fact3 = 0;
	double fact4 = 0;
	double fact5 = 0;

	double alph1 = 0;
	double alph2 = 0;
	double beta1 = 0;
	double beta2 = 0;

	if (GDRPar.size() == 0)
	{
		return cross;
	}

	for (unsigned int i = 0; i < GDRPar.size(); i++)
	{
		if (GDRPar[i][0] == Z && GDRPar[i][1] == A)
		{

			fact1 = (E * E) / (4 * pi * hbar * c);

			fact2 = GDRPar[i][4] * pow(10, -31) * GDRPar[i][5] * (GDRPar[i][3] * GDRPar[i][3] - E * E) /
					((GDRPar[i][3] * GDRPar[i][3] - E * E) * (GDRPar[i][3] * GDRPar[i][3] - E * E) +
					 E * E * GDRPar[i][5] * GDRPar[i][5]);

			fact3 = GDRPar[i][7] * pow(10, -31) * GDRPar[i][8] * (GDRPar[i][6] * GDRPar[i][6] - E * E) /
					((GDRPar[i][6] * GDRPar[i][6] - E * E) * (GDRPar[i][6] * GDRPar[i][6] - E * E) +
					 E * E * GDRPar[i][8] * GDRPar[i][8]);

			if (E < GDRPar[i][2])
			{
				fact4 = 0;
				fact5 = 0;
			}
			else
			{

				fact4 = GDRPar[i][4] * pow(10, -31) * GDRPar[i][5] * (E * GDRPar[i][5]) / ((GDRPar[i][3] * GDRPar[i][3] - E * E) * (GDRPar[i][3] * GDRPar[i][3] - E * E) + E * E * GDRPar[i][5] * GDRPar[i][5]);

				fact5 = GDRPar[i][7] * pow(10, -31) * GDRPar[i][8] * (E * GDRPar[i][8]) / ((GDRPar[i][6] * GDRPar[i][6] - E * E) * (GDRPar[i][6] * GDRPar[i][6] - E * E) + E * E * GDRPar[i][8] * GDRPar[i][8]);
			}

			alph1 = fact1 * fact2;
			beta1 = fact1 * fact4;
			alph2 = fact1 * fact3;
			beta2 = fact1 * fact5;

			cross = pow(1, 2) * pow(GDRPar[i][9], 2) * (pow((2 * alph1 - alph2), 2) + pow((2 * beta1 - beta2), 2)) * (1.0 / 40.0) * (13.0 + pow(cos(thet), 2));

			return cross;
		}
	}

	return cross;
}

Amplitude G4Elastic::Thomson(double Z, double A, double, double thet)
{

	Amplitude Amp;

	Amp.AparReal = 0;
	Amp.AparImag = 0;
	Amp.AppReal = 0;
	Amp.AppImag = 0;

	double mass = A * Mn;
	double fact1 = -(Z * Z * Me * Elrad) / mass;

	Amp.AppReal = fact1;
	Amp.AppImag = 0;
	Amp.AparReal = fact1 * cos(thet);
	Amp.AparImag = 0;

	return Amp;
}

Amplitude G4Elastic::Delbruck(double Z, double, double E, double thet)
{

	Amplitude Amp;

	Amp.AparReal = 0;
	Amp.AparImag = 0;
	Amp.AppReal = 0;
	Amp.AppImag = 0;

	if (E < 0.2555)
	{
		return Amp;
	}

	thet = thet * 180.0 / pi;

	double fact1 = Z * Z * alph * alph * Elrad;

	for (unsigned int i = 0; i < DelAmp.size() - 1; i++)
	{

		if (E >= DelAmp[i][0] && E < DelAmp[i + 1][0])
		{

			for (double o = 1; o < 75; o += 5)
			{

				if (thet >= DelAmp[i][o] && thet < DelAmp[i][o + 5])
				{

					double temp1 = Interpolate(DelAmp[i][o], DelAmp[i][o + 1], DelAmp[i][o + 5], DelAmp[i][o + 5 + 1], thet);
					double temp2 = Interpolate(DelAmp[i + 1][o], DelAmp[i + 1][o + 1], DelAmp[i + 1][o + 5], DelAmp[i + 1][o + 5 + 1], thet);
					double temp3 = Interpolate(DelAmp[i][0], temp1, DelAmp[i + 1][0], temp2, E);

					double temp4 = Interpolate(DelAmp[i][o], DelAmp[i][o + 2], DelAmp[i][o + 5], DelAmp[i][o + 5 + 2], thet);
					double temp5 = Interpolate(DelAmp[i + 1][o], DelAmp[i + 1][o + 2], DelAmp[i + 1][o + 5], DelAmp[i + 1][o + 5 + 2], thet);
					double temp6 = Interpolate(DelAmp[i][0], temp4, DelAmp[i + 1][0], temp5, E);

					double temp7 = Interpolate(DelAmp[i][o], DelAmp[i][o + 3], DelAmp[i][o + 5], DelAmp[i][o + 5 + 3], thet);
					double temp8 = Interpolate(DelAmp[i + 1][o], DelAmp[i + 1][o + 3], DelAmp[i + 1][o + 5], DelAmp[i + 1][o + 5 + 3], thet);
					double temp9 = Interpolate(DelAmp[i][0], temp7, DelAmp[i + 1][0], temp8, E);

					double temp10 = Interpolate(DelAmp[i][o], DelAmp[i][o + 4], DelAmp[i][o + 5], DelAmp[i][o + 5 + 4], thet);
					double temp11 = Interpolate(DelAmp[i + 1][o], DelAmp[i + 1][o + 4], DelAmp[i + 1][o + 5], DelAmp[i + 1][o + 5 + 4], thet);
					double temp12 = Interpolate(DelAmp[i][0], temp10, DelAmp[i + 1][0], temp11, E);

					Amp.AppReal = fact1 * temp3;
					Amp.AppImag = fact1 * temp6;
					Amp.AparReal = fact1 * temp9;
					Amp.AparImag = fact1 * temp12;

					return Amp;
				}
			}
		}
	}

	return Amp;
}

double G4Elastic::Interpolate(double x1, double y1, double x2, double y2, double value)
{

	double slope = (y2 - y1) / (x2 - x1);
	double intercept = y1 - slope * x1;
	return (slope * value + intercept);
}

double G4Elastic::GetRandom()
{

	vector<double> *TestIntegrals = GetIntegrals();

	double r1 = gRandom->Rndm();

	std::vector<double>::iterator pind;

	pind = std::lower_bound((*TestIntegrals).begin(), (*TestIntegrals).end(), r1);

	int newBin;

	if ((pind != (*TestIntegrals).end()) && (*pind == r1))
	{
		newBin = std::distance((*TestIntegrals).begin(), pind);
	}
	else
	{
		newBin = std::distance((*TestIntegrals).begin(), pind - 1);
	}

	double BinWidth = 181.0 / GranT;

	double x = newBin * BinWidth;

	if (r1 > (*TestIntegrals)[newBin])
	{
		x += 1.0 * BinWidth * (r1 - (*TestIntegrals)[newBin]) / ((*TestIntegrals)[newBin + 1] - (*TestIntegrals)[newBin]);
	}

	return x;
}

double G4Elastic::Integrate(vector<double> *array)
{

	double integral = 0;

	int first = 0;
	int last = (*array).size();
	int np = last - first + 1;
	double sum = 0;

	for (unsigned int i = 0; i < (*array).size(); i++)
	{
		integral += (*array)[i];
	}

	for (int i = first; i <= last; i++)
	{

		int j = first + (i - first + 1) % np;
		sum += ((*array)[i] + (*array)[j]) * (i - (j));
	}

	return 0.5 * TMath::Abs(sum);
}
