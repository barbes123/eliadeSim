#include "G4SDManager.hh"
#include "AppsRun.hh"
#include <iostream>
#include "AppsInput.hh"
#include "AppsHit.hh"
#include "G4MTRunManager.hh"
#include "G4Threading.hh"
#include <limits>

AppsRun::AppsRun()
{

	// construct vector to hold the digital data

	AppsInput *FInput = AppsInput::Instance();

	string BeamType = FInput->GetBeamType();

	if (BeamType == "source")
	{
		mergeTime = 1;
	}

	G4int AddCnt = 0;

	mergeCounter = vector<G4int>(FInput->GetNumberOfCores(), 0);

	BuildAnalysisObjects();
}

AppsRun::~AppsRun()
{
}

void AppsRun::Merge(const G4Run *aRun)
{

	mergeCounter[G4Threading::G4GetThreadId()]++;

	const AppsRun *localRun = static_cast<const AppsRun *>(aRun);

	G4bool doMerge = false;

	for (unsigned int i = 0; i < AnalysisObjects.size(); i++)
	{

		for (unsigned int o = 0; o < AnalysisObjects[i].size(); o++)
		{

			AnalysisObjects[i][o]->rawEvents.insert(AnalysisObjects[i][o]->rawEvents.end(), localRun->AnalysisObjects[i][o]->rawEvents.begin(),
													localRun->AnalysisObjects[i][o]->rawEvents.end());

			localRun->AnalysisObjects[i][o]->ClearRawEvents();

			if (AnalysisObjects[i][o]->rawEvents.size() > rawEventsBufferSize)
			{
				doMerge = true;
			}
		}
	}

	unsigned long long int currentMerge = std::numeric_limits<unsigned long long>::max();

	for (unsigned int j = 0; j < mergeCounter.size(); j++)
	{

		if (mergeCounter[j] < currentMerge)
		{

			currentMerge = mergeCounter[j];
		}
	}

	if (currentMerge > lastMerge && doMerge)
	{

		for (unsigned int i = 0; i < AnalysisObjects.size(); i++)
		{

			for (unsigned int o = 0; o < AnalysisObjects[i].size(); o++)
			{

				//	G4cout << " i "  << i << " o " << o << " before build " << AnalysisObjects[i][o]->rawEvents.size() << endl;
				AnalysisObjects[i][o]->BuildDigitalEvents((currentMerge * mergeIncrement) * CLHEP::s);
				//	G4cout << " i "  << i << " o " << o << " after build " << AnalysisObjects[i][o]->digitalEvents.size() << endl;

				AnalysisObjects[i][o]->RecordDigitalEvents();

			}
		}

		lastMerge = currentMerge;
	}
}

void AppsRun::RecordEvent(const G4Event *evt)
{

	G4int evtNb = evt->GetEventID();

	if (evtNb%10000 == 0) { 
		G4cout << "\n---> end of event: " << evtNb << G4endl;
	}
	
	 
	beamTime = evt->GetPrimaryVertex()->GetT0();

	if (beamTime > (mergeTime * 1) * CLHEP::s)
	{

		G4MTRunManager *mtRM = G4MTRunManager::GetMasterRunManager();

		finalRun = false;

		const AppsRun *localRun = static_cast<const AppsRun *>(G4RunManager::GetRunManager()->GetCurrentRun());
		mtRM->MergeRun(localRun);

		finalRun = true;

		mergeTime += mergeIncrement;
	}

	// Get hits collections of this event
	G4HCofThisEvent *HCE = evt->GetHCofThisEvent();

	for (int i = 0; i < HCE->GetNumberOfCollections(); i++)
	{

		G4VHitsCollection *hc = HCE->GetHC(i);

		if (hc->GetSize() == 0)
		{
			continue;
		}
		else
		{

			tuple<double, double, int> event = AnalysisObjects[i][0]->ProcessEvent(hc);

			if (get<1>(event) != 0)
			{

				for (unsigned int o = 0; o < AnalysisObjects[i].size(); o++)
				{

					// add the event to the list

					if (o == 0 && AnalysisObjects[i][0]->WasHit == false)
					{

						HitList.push_back(i);
					}

					AnalysisObjects[i][o]->AddEvent(event);
				}
			}
		}
	}

	for (unsigned int i = 0; i < HitList.size(); i++)
	{
		for (unsigned int o = 0; o < AnalysisObjects[HitList[i]].size(); o++)
		{
			AnalysisObjects[HitList[i]][o]->RecordSumEvents();
		}
	}
	HitList.clear();
}

void AppsRun::BuildAnalysisObjects()
{

	AppsInput *FInput = AppsInput::Instance();

	vector<AppsGeometryObject> objects = FInput->GetGeometryObjects();

	int nmod{0};

	for (unsigned int i = 0; i < objects.size(); i++)
	{

		if (FInput->GetGeometryObjects()[i].GetObjectActiveStatus())
		{

			if (objects[i].GetObjectType() == "ELIADE")
			{

				G4bool isSegmented = false;
				for (unsigned int o = 0; o < objects[i].GetObjectFeatures().size(); o++)
				{

					if (objects[i].GetObjectFeatures()[o] == "segmented")
					{
						isSegmented = true;
					}
				}

				G4bool isDigital = false;
				for (unsigned int o = 0; o < objects[i].GetObjectHistogramType().size(); o++)
				{

					if (objects[i].GetObjectHistogramType()[o] == "DEdep" || objects[i].GetObjectHistogramType()[o] == "AddBack" || objects[i].GetObjectHistogramType()[o] == "Pileup")
					{
						isDigital = true;
						break;
					}
				}

				AppsAnalysisObject *arranalysis = new AppsAnalysisObject();
				arranalysis->SetObjectName(objects[i].GetObjectName());
				arranalysis->SetDigitalFlag(isDigital);
				arranalysis->SetObjectId(i);

				AppsAnalysisObject *detanalysis = new AppsAnalysisObject();
				detanalysis->SetObjectName(objects[i].GetObjectName() + "_Clover0");
				detanalysis->SetObjectId(i);

				AppsAnalysisObject *cryanalysis = new AppsAnalysisObject();
				cryanalysis->SetObjectName(objects[i].GetObjectName() + "_Clover0" + "_Crystal0");
				cryanalysis->SetObjectId(i);

				if (isSegmented)
				{

					for (unsigned int o = 0; o < objects[i].GetObjectActiveVolumes().size(); o++)
					{

						vector<AppsAnalysisObject *> temp;

						AppsAnalysisObject *seganalysis = new AppsAnalysisObject();
						seganalysis->SetObjectName(objects[i].GetObjectName() + "_Clover" + to_string(o / 32) +
												   "_Crystal" + to_string(o / 8 - 4 * (o / 32)) + "_Segment" + to_string((o % 8) + 1));
						seganalysis->SetObjectId(i);

						if (o % 8 == 0)
						{

							cryanalysis = new AppsAnalysisObject();
							cryanalysis->SetObjectName(objects[i].GetObjectName() + "_Clover" + to_string(o / 32) +
													   "_Crystal" + to_string(o / 8 - 4 * (o / 32)));
							cryanalysis->SetObjectId(i);
							cryanalysis->SetDigitalFlag(isDigital);
						}

						if (o % 32 == 0)
						{

							detanalysis = new AppsAnalysisObject();
							detanalysis->SetObjectName(objects[i].GetObjectName() + "_Clover" + to_string(o / 32));
							detanalysis->SetDigitalFlag(isDigital);
							detanalysis->SetObjectId(i);
						}

						seganalysis->SetDigitalFlag(isDigital);

						temp.push_back(seganalysis);
						temp.push_back(cryanalysis);
						temp.push_back(detanalysis);
						temp.push_back(arranalysis);

						for (unsigned int j = 0; j < temp.size(); j++)
						{

							temp[j]->SetObjectHistograms(objects[i].GetObjectHistogramType());
							temp[j]->SetDetectorResolution(objects[i].GetObjectFWHMParam());
						}

						AnalysisObjects.push_back(temp);
					}
				}
				else
				{

					for (unsigned int o = 0; o < objects[i].GetObjectActiveVolumes().size(); o++)
					{

						vector<AppsAnalysisObject *> temp;

						cryanalysis = new AppsAnalysisObject();
						cryanalysis->SetObjectName(objects[i].GetObjectName() + "_Clover" + to_string(o / 4) +
												   "_Crystal" + to_string(o % 4));
						cryanalysis->SetObjectId(i);

						if (o % 4 == 0)
						{

							detanalysis = new AppsAnalysisObject();
							detanalysis->SetObjectName(objects[i].GetObjectName() + "_Clover" + to_string(o / 4));
							detanalysis->SetDigitalFlag(isDigital);
							detanalysis->SetObjectId(i);
						}

						cryanalysis->SetDigitalFlag(isDigital);

						temp.push_back(cryanalysis);
						temp.push_back(detanalysis);
						temp.push_back(arranalysis);

						for (unsigned int j = 0; j < temp.size(); j++)
						{

							temp[j]->SetObjectHistograms(objects[i].GetObjectHistogramType());
							temp[j]->SetDetectorResolution(objects[i].GetObjectFWHMParam());
						}

						AnalysisObjects.push_back(temp);
					}
				}

				continue;
			}

			if (objects[i].GetObjectType() == "PIXEL")
			{

				G4bool isDigital = false;

				for (unsigned int o = 0; o < objects[i].GetObjectHistogramType().size(); o++)
				{

					if (objects[i].GetObjectHistogramType()[o] == "DEdep" || objects[i].GetObjectHistogramType()[o] == "AddBack" ||
						objects[i].GetObjectHistogramType()[o] == "Pileup")
					{
						isDigital = true;
						break;
					}
				}

				for (unsigned int o = 0; o < objects[i].GetObjectActiveVolumes().size(); o++)
				{

					vector<AppsAnalysisObject *> temp;

					AppsAnalysisObject *pixAnalysis = new AppsAnalysisObject();
					pixAnalysis->SetObjectName(objects[i].GetObjectActiveVolumes()[o]);
					pixAnalysis->SetDigitalFlag(isDigital);

					temp.push_back(pixAnalysis);

					for (unsigned int j = 0; j < temp.size(); j++)
					{

						temp[j]->SetObjectHistograms(objects[i].GetObjectHistogramType());
						temp[j]->SetDetectorResolution(objects[i].GetObjectFWHMParam());
					}

					AnalysisObjects.push_back(temp);
				}

				continue;

			}

			if (objects[i].GetObjectType() == "CloverDetector")
			{

				G4bool isSegmented = false;
				for (unsigned int o = 0; o < objects[i].GetObjectFeatures().size(); o++)
				{

					if (objects[i].GetObjectFeatures()[o] == "segmented")
					{
						isSegmented = true;
					}
				}

				G4bool isShielded = false;
				for (unsigned int o = 0; o < objects[i].GetObjectFeatures().size(); o++)
				{

					if (objects[i].GetObjectFeatures()[o] == "shielded")
					{
						isShielded = true;
					}
				}

				G4bool isDigital = false;
				for (unsigned int o = 0; o < objects[i].GetObjectHistogramType().size(); o++)
				{

					if (objects[i].GetObjectHistogramType()[o] == "DEdep" || objects[i].GetObjectHistogramType()[o] == "AddBack" ||
						objects[i].GetObjectHistogramType()[o] == "Pileup")
					{
						isDigital = true;
						break;
					}
				}

				AppsAnalysisObject *detanalysis = new AppsAnalysisObject();
				detanalysis->SetObjectName(objects[i].GetObjectName());
				detanalysis->SetDigitalFlag(isDigital);

				detanalysis->SetObjectId((i + 1) * pow(10, 9));

				AppsAnalysisObject *cryanalysis = new AppsAnalysisObject();
				cryanalysis->SetObjectName(objects[i].GetObjectName() + "_Crystal0");

				if (isSegmented)
				{
				
					if (isShielded)
					{

						for (unsigned int o = 0; o < objects[i].GetObjectActiveVolumes().size(); o++)
						{
							vector<AppsAnalysisObject *> temp;
							
							if (o % 42 < 32)
							{
cout<<"ici "<<objects[i].AnalysisObjects[o][0]->GetObjectName()<<" "<<objects[i].GetObjectName() + "_Crystal" + to_string(o / 8) + "_Segment" + to_string((o % 8) + 1)<<endl;

								AppsAnalysisObject *seganalysis = new AppsAnalysisObject();
								seganalysis->SetObjectName(objects[i].GetObjectName() + "_Crystal" + to_string(o / 8) + "_Segment" + to_string((o % 8) + 1));

								seganalysis->SetObjectId((i + 1) * pow(10, 9) + (o / 8 + 1) * pow(10, 6) + (o % 8 + 1) * pow(10, 3));
								seganalysis->SetMod(nmod);
								seganalysis->SetCh(o % 8 + 1);

								if (o % 8 == 0)
								{

									cryanalysis = new AppsAnalysisObject();
									cryanalysis->SetObjectName(objects[i].GetObjectName() + "_Crystal" + to_string(o / 8));
									cryanalysis->SetDigitalFlag(isDigital);
									cryanalysis->SetObjectId((i + 1) * pow(10, 9) + (o / 8 + 1) * pow(10, 6));
									cryanalysis->SetMod(nmod);
									cryanalysis->SetCh(0);
								}

								if (o % 8 == 7)
								{
									nmod++;
								}

								seganalysis->SetDigitalFlag(isDigital);

								temp.push_back(seganalysis);
								temp.push_back(cryanalysis);
								temp.push_back(detanalysis);
							}
							else
							{
								
								AppsAnalysisObject *shianalysis = new AppsAnalysisObject();
								
								G4String side;
								
								if (o % 42 < 42){side=string("Back");}
								if (o % 42 < 40){side=string("Side");}
								if (o % 42 < 36){side=string("Front");}
cout<<"ici "<<objects[i].AnalysisObjects[o][0]->GetObjectName()<<" "<<objects[i].GetObjectName() + "_" + side + "Shield" + to_string((o % 4) + 1)<<endl;

								shianalysis->SetObjectName(objects[i].GetObjectName() + "_" + side + "Shield" + to_string((o % 4) + 1));
							
								shianalysis->SetObjectId((i + 1) * pow(10, 9) + (o / 4 + 1) * pow(10, 6) + (o % 4 + 1) * pow(10, 3));
								shianalysis->SetMod(nmod);
								shianalysis->SetCh(o % 8 + 1);
								
								if (o % 42 == 39 || o % 42 == 41 )
								{
									nmod++;
								}

								shianalysis->SetDigitalFlag(isDigital);

								temp.push_back(shianalysis);
								//temp.push_back(cryanalysis);
								temp.push_back(detanalysis);
							}

							AnalysisObjects.push_back(temp);
						}
					}
					else
					{
						for (unsigned int o = 0; o < objects[i].GetObjectActiveVolumes().size(); o++)
						{

							vector<AppsAnalysisObject *> temp;

							AppsAnalysisObject *seganalysis = new AppsAnalysisObject();
							seganalysis->SetObjectName(objects[i].GetObjectName() + "_Crystal" + to_string(o / 8) + "_Segment" + to_string((o % 8) + 1));

							seganalysis->SetObjectId((i + 1) * pow(10, 9) + (o / 8 + 1) * pow(10, 6) + (o % 8 + 1) * pow(10, 3));
							seganalysis->SetMod(nmod);
							seganalysis->SetCh(o % 8 + 1);

							if (o % 8 == 0)
							{

								cryanalysis = new AppsAnalysisObject();
								cryanalysis->SetObjectName(objects[i].GetObjectName() + "_Crystal" + to_string(o / 8));
								cryanalysis->SetDigitalFlag(isDigital);
								cryanalysis->SetObjectId((i + 1) * pow(10, 9) + (o / 8 + 1) * pow(10, 6));
								cryanalysis->SetMod(nmod);
								cryanalysis->SetCh(0);
							}

							if (o % 8 == 7)
							{
								nmod++;
							}

							seganalysis->SetDigitalFlag(isDigital);

							temp.push_back(seganalysis);
							temp.push_back(cryanalysis);
							temp.push_back(detanalysis);

							for (unsigned int j = 0; j < temp.size(); j++)
							{

								temp[j]->SetObjectHistograms(objects[i].GetObjectHistogramType());
								temp[j]->SetDetectorResolution(objects[i].GetObjectFWHMParam());
							}

							AnalysisObjects.push_back(temp);
						}
					}
				}
				else
				{

					for (unsigned int o = 0; o < objects[i].GetObjectActiveVolumes().size(); o++)
					{

						vector<AppsAnalysisObject *> temp;

						cryanalysis = new AppsAnalysisObject();
						cryanalysis->SetObjectName(objects[i].GetObjectName() + "_Crystal" + to_string(o));

						cryanalysis->SetObjectId((i + 1) * pow(10, 9) + (o + 1) * pow(10, 6));

						cryanalysis->SetDigitalFlag(isDigital);

						temp.push_back(cryanalysis);
						temp.push_back(detanalysis);

						for (unsigned int j = 0; j < temp.size(); j++)
						{

							temp[j]->SetObjectHistograms(objects[i].GetObjectHistogramType());
							temp[j]->SetDetectorResolution(objects[i].GetObjectFWHMParam());
						}

						AnalysisObjects.push_back(temp);
					}
				}

				continue;
			}

			G4bool isDigital = false;
			for (unsigned int o = 0; o < objects[i].GetObjectHistogramType().size(); o++)
			{

				if (objects[i].GetObjectHistogramType()[o] == "DEdep" || objects[i].GetObjectHistogramType()[o] == "AddBack" ||
					objects[i].GetObjectHistogramType()[o] == "Pileup")
				{
					isDigital = true;
					break;
				}
			}

			vector<AppsAnalysisObject *> temp;
			AppsAnalysisObject *arranalysis = new AppsAnalysisObject();
			arranalysis->SetObjectName(objects[i].GetObjectName());

			arranalysis->SetObjectId(i);
			arranalysis->SetMod(-2);
			arranalysis->SetObjectHistograms(objects[i].GetObjectHistogramType());
			arranalysis->SetDigitalFlag(isDigital);
			arranalysis->SetDetectorResolution(objects[i].GetObjectFWHMParam());

			if (objects[i].GetObjectType() == "ScintDetector" || objects[i].GetObjectType() == "CeBrDetector")
			{
				arranalysis->SetLinearResFlag(false);
			}

			temp.push_back(arranalysis);
			AnalysisObjects.push_back(temp);
		}
	}
}
