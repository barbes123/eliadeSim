#include "AppsAnalysisObject.hh"
#include "AppsHit.hh"
#include "AppsAnalysisManager.hh"
#include "G4MTRunManager.hh"
#include <chrono>

AppsAnalysisObject::AppsAnalysisObject()
{

	Initialize();
}

AppsAnalysisObject::~AppsAnalysisObject()
{
}

void AppsAnalysisObject::Initialize()
{

	detectorRes.push_back(0);
	detectorRes.push_back(0);
}

void AppsAnalysisObject::RecordSumEvents()
{

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	if (sumEvents.size() == 0)
	{
		return;
	}

	G4double totalEdep = 0;

	G4double time = 0;

	G4double pileFlag = 0;

	for (unsigned int i = 0; i < sumEvents.size(); i++)
	{

		totalEdep += get<1>(sumEvents[i]);
		time += get<0>(sumEvents[i]);
		pileFlag = 1;
	}

	time = (time * 1.0) / sumEvents.size();

	if (digitalFlag)
	{
		rawEvents.push_back(make_tuple(time, totalEdep, pileFlag));
	}

	sumEvents.clear();
	WasHit = false;
}

void AppsAnalysisObject::AddAddBackEvents(vector<tuple<double, double, int>> events)
{

	digitalEvents.insert(digitalEvents.end(), events.begin(), events.end());
}

void AppsAnalysisObject::RecordDigitalEvents()
{

	if (G4Threading::G4GetThreadId() != -1)
	{
		return;
	}

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	G4double AddWindow = 0.02 * CLHEP::ms;

	std::reverse(digitalEvents.begin(), digitalEvents.end());

	while (digitalEvents.size() != 0)
	{

		if (get<1>(digitalEvents.back()) != 0)
		{

			G4double energy = get<1>(digitalEvents.back());

			G4int pileFlag = get<2>(digitalEvents.back());

			analysisManager->FillELIADETreeColumn(0, Mod);
			if (Ch==-1){
				analysisManager->FillELIADETreeColumn(1, Ch);
			}
			else{
				analysisManager->FillELIADETreeColumn(1, Ch+1);
			}
			analysisManager->FillELIADETreeColumn(2, long(get<0>(digitalEvents.back())));
			analysisManager->FillELIADETreeColumn(3, get<0>(digitalEvents.back())*1000);
			analysisManager->FillELIADETreeColumn(4, get<1>(digitalEvents.back())*1000);
			analysisManager->AddELIADETreeRow();
			
			if (Ch==0){
				analysisManager->FillELIADETreeColumn(0, Mod);
				analysisManager->FillELIADETreeColumn(1, 0);
				analysisManager->FillELIADETreeColumn(2, get<0>(digitalEvents.back()));
				analysisManager->FillELIADETreeColumn(3, get<0>(digitalEvents.back())*1000);
				analysisManager->FillELIADETreeColumn(4, get<1>(digitalEvents.back())*1000);
				analysisManager->AddELIADETreeRow();
			}

			//	AddAddBackEvent(digitalEvents.back());
		}

		digitalEvents.pop_back();
	}
}

void AppsAnalysisObject::RecordAddBackEvents()
{

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	sort(digitalEvents.rbegin(), digitalEvents.rend());

	while (digitalEvents.size() != 0)
	{
		G4int pileFlag = get<2>(digitalEvents.back());
		G4double energy = get<1>(digitalEvents.back());

		analysisManager->FillELIADETreeColumn(0, objectId);
		analysisManager->FillELIADETreeColumn(1, objectId);
		analysisManager->FillELIADETreeColumn(2, get<0>(digitalEvents.back()));
		analysisManager->FillELIADETreeColumn(3, get<0>(digitalEvents.back()));
		analysisManager->FillELIADETreeColumn(4, get<1>(digitalEvents.back()));
		analysisManager->FillELIADETreeColumn(5, get<1>(digitalEvents.back()));
		analysisManager->AddELIADETreeRow();
G4cout<<objectId<<" "<<get<1>(digitalEvents.back())<<G4endl;
		digitalEvents.pop_back();
	}
}

void AppsAnalysisObject::AddEvent(tuple<double, double, int> event)
{

	WasHit = true;

	sumEvents.push_back(event);
}

void AppsAnalysisObject::ClearRawEvents()
{

	rawEvents.clear();
}

void AppsAnalysisObject::ClearDigitalEvents()
{

	digitalEvents.clear();
	fullBuffer = false;
}

tuple<double, double, int> AppsAnalysisObject::ProcessEvent(G4VHitsCollection *hitCol)
{

	AppsAnalysisManager *analysisManager = AppsAnalysisManager::Instance();

	G4double totalEdep = 0;
	G4double averageTime = 0;
	vector<G4double> faceEnergy;
	vector<G4ThreeVector> position;

	// extract data from event

	for (unsigned int n = 0; n < hitCol->GetSize(); n++)
	{

		AppsHit *hit = (AppsHit *)hitCol->GetHit(n);
		totalEdep += hit->GetEdep();

		averageTime += hit->GetTime();
		faceEnergy.push_back(hit->GetFaceEnergy());
		position.push_back(hit->GetPosition());
	}

	averageTime /= hitCol->GetSize();

	// add detector resolution

	if (totalEdep != 0)
	{
		totalEdep = DetectorResolution(totalEdep);
		if (totalEdep < 0)
		{
			totalEdep = 0;
		}
	}

	// add monitor specific data to histograms, the monitor data is not passed along with the event

	for (unsigned int i = 0; i < objectHistograms.size(); i++)
	{

		if (objectHistograms[i] == "Eflux")
		{

			//G4int histID = analysisManager->GetH1Id(objectName + "_" + objectHistograms[i]);

			for (unsigned int l = 0; l < faceEnergy.size(); l++)
			{

				//analysisManager->FillH1(histID, faceEnergy[l]);
			}
			continue;
		}

		if (objectHistograms[i] == "Sflux")
		{

			//G4int histID = analysisManager->GetH2Id(objectName + "_" + objectHistograms[i]);

			for (unsigned int l = 0; l < position.size(); l++)
			{

				//analysisManager->FillH2(histID, position[l].getX(), position[l].getY());
			}
			continue;
		}
	}

	tuple<double, double, int> event = make_tuple(averageTime, totalEdep, 1);

	return event;
}

G4double AppsAnalysisObject::DetectorResolution(G4double energy)
{

	G4double menergy;

	if (detectorRes[0] == 0 && detectorRes[1] == 0)
	{

		return energy;
	}
	else
	{

		if (linearRes)
		{

			G4double sigma = detectorRes[0] * energy + detectorRes[1];
			menergy = G4RandGauss::shoot(energy, sigma);
		}
		else
		{

			G4double FWHM = detectorRes[0] * pow(energy, detectorRes[1]);
			G4double sigma = (energy * FWHM) / (2.355 * 100);
			menergy = G4RandGauss::shoot(energy, sigma);
		}
	}

	return menergy;
}

void AppsAnalysisObject::AddDigitalEvents(vector<tuple<double, double, int>> events)
{

	digitalEvents.insert(digitalEvents.end(), events.begin(), events.end());
}

void AppsAnalysisObject::BuildDigitalEvents(unsigned long long currentTime)
{

	if (G4Threading::G4GetThreadId() != -1)
	{
		return;
	}

	if (rawEvents.size() == 0)
	{
		return;
	}

	//	 auto start = std::chrono::high_resolution_clock::now();

	vector<tuple<double, double, int>> workVector = rawEvents;

	vector<tuple<double, double, int>> builtEvents;

	sort(workVector.rbegin(), workVector.rend());

	//	double CorFact = 2.125;
	double CorFact = 1;
	//	double peakHoldOff = 8000;

	// the old digites parameters
	/*
		double peakHoldOff = 1000;
		double trapFlatTop = 1500;
		double trapRiseTime = 3500;
		double inputRiseTime = 200*CorFact;		// value calculated based on the TTF delay from digites
	*/

	// the new compass parameters

	double peakHoldOff = 960;
	double trapFlatTop = 2000;
	double trapRiseTime = 5000;

	// rise time corrections

	double expRiseTime = 400;
	double intercept = -59.6943;
	double slope = 0.51188;

	double inputRiseTime = (1 / slope * expRiseTime - (intercept / slope)) * 1.5;

	// rise time corrections

	//	cout << inputRiseTime << endl;

	//	double inputRiseTime = 850*1.5*CorFact;		// the so called inputrisetime seems to be the TTF delay which is calculated as the rise time of the input times 1.5

	double pkRun = trapRiseTime + trapFlatTop + peakHoldOff;
	double trapTime = trapRiseTime + trapFlatTop;

	double deltaT = 0;

	if (currentTime != 0)
	{
		while (get<0>(workVector.back()) < (currentTime - pkRun) && workVector.size() > 1)
		{

			tuple<double, double, int> currentEvent = workVector.back();
			workVector.pop_back();

			deltaT = get<0>(workVector.back()) - get<0>(currentEvent);

			if (!(deltaT > pkRun))
			{

				if (!(deltaT > trapTime))
				{

					if (!(deltaT > inputRiseTime))
					{

						double time = (get<0>(workVector.back()) + get<0>(currentEvent)) / 2.0;

						double energy = get<1>(workVector.back()) + get<1>(currentEvent);

						/*if (energy > 19)
						{
							energy = 19;
						}*/

						int pileUpFlag = get<2>(workVector.back()) + get<2>(currentEvent);

						workVector.pop_back();
						workVector.push_back(make_tuple(time, energy, pileUpFlag));
					}
					else
					{

						get<1>(currentEvent) = 0;

						builtEvents.push_back(currentEvent);

						get<1>(workVector.back()) = 0;
					}
				}
				else
				{

					builtEvents.push_back(currentEvent);

					get<1>(workVector.back()) = 0;

					get<2>(workVector.back())++;
				}
			}
			else
			{

				builtEvents.push_back(currentEvent);
			}
		}

		rawEvents.clear();

		for (unsigned int i = 0; i < workVector.size(); i++)
		{
			rawEvents.push_back(workVector[i]);
		}

		digitalEvents.clear();
		digitalEvents = builtEvents;
	}
	else
	{

		while (workVector.size() > 1)
		{

			tuple<double, double, int> currentEvent = workVector.back();
			workVector.pop_back();

			deltaT = get<0>(workVector.back()) - get<0>(currentEvent);

			if (!(deltaT > pkRun))
			{

				if (!(deltaT > trapTime))
				{

					if (!(deltaT > inputRiseTime))
					{

						double time = (get<0>(workVector.back()) + get<0>(currentEvent)) / 2.0;

						double energy = get<1>(workVector.back()) + get<1>(currentEvent);

						/*if (energy > 19)
						{
							energy = 19;
						}*/

						int pileUpFlag = get<2>(workVector.back()) + get<2>(currentEvent);

						workVector.pop_back();
						workVector.push_back(make_tuple(time, energy, pileUpFlag));
					}
					else
					{

						get<1>(currentEvent) = 0;

						builtEvents.push_back(currentEvent);

						get<1>(workVector.back()) = 0;
					}
				}
				else
				{

					builtEvents.push_back(currentEvent);

					get<1>(workVector.back()) = 0;

					get<2>(workVector.back())++;
				}
			}
			else
			{

				builtEvents.push_back(currentEvent);
			}
		}

		rawEvents.clear();

		for (unsigned int i = 0; i < workVector.size(); i++)
		{
			builtEvents.push_back(workVector[i]);
		}

		digitalEvents.clear();
		digitalEvents = builtEvents;
	}
}
