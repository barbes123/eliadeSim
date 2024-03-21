///-------------- Turturica Gabriel  ------- 12.02.2016 --- ////

#include <fstream>
#include "AppsInput.hh"
#include <stdlib.h>
#include <sstream>

using namespace std;

AppsInput *AppsInput::fgInstance = 0;

AppsInput::AppsInput()
{
	Initialize();
}

AppsInput::~AppsInput()
{
}

void AppsInput::Initialize()
{

	InitializeParameters();
	ReadRunParameters();

	//	PrintSimParameters();
	ReadGeometryObjects();
}

void AppsInput::ReadGeometryObjects()
{

	G4String z;
	ifstream inputfile("InputObjects.txt");

	G4bool newobject = false;
	G4bool sameobject = false;

	if (inputfile.is_open())
	{
		while (!inputfile.eof())
		{

			getline(inputfile, z);

			if (z == "")
			{
				newobject = false;
				sameobject = false;
			}

			for (unsigned int i = 0; i < InputGeometryObjects.size(); i++)
			{
				std::size_t Objfound = z.find(InputGeometryObjects[i]);

				if (Objfound != std::string::npos)
				{

					AppsGeometryObject object(InputGeometryObjects[i]);
					GeometryObjects.push_back(object);

					newobject = true;
					sameobject = false;
				}
				else
				{
					if (newobject)
					{
						newobject = false;
						sameobject = true;
					}
				}
			}

			if (newobject || sameobject)
			{

				for (unsigned int o = 0; o < InputObjectParameters.size(); o++)
				{
					std::size_t Parfound = z.find(InputObjectParameters[o]);

					if (Parfound != std::string::npos)
					{

						std::size_t firstpos = z.find("[") + 1;
						std::size_t secondpos = z.find("]");

						stringstream ss;

						ss.str(z.substr(firstpos, (secondpos - firstpos)));

						string x;

						switch (o)
						{

						case 0:
						{ // position

							vector<G4double> position;

							while (std::getline(ss, x, ','))
							{
								position.push_back(stod(x));
							}

							if (GeometryObjects.back().GetObjectType() == "World")
							{

								if (position.size() == 3)
								{

									if ((position[0] || position[1] || position[2]) != 0)
									{

										cout << " " << endl;
										cout << "--------------------------------- Warning ---------------------------------" << endl;
										cout << "\"" << GeometryObjects.back().GetObjectType() << "\""
											 << " position must be defined by 3 zero parameters,";
										cout << " using the default (0,0,0) values." << endl;

										GeometryObjects.back().SetObjectPosition(0 * cm, 0 * cm, 0 * cm);
									}
									else
									{
										GeometryObjects.back().SetObjectPosition(position[0] * cm, position[1] * cm, position[2] * cm);
									}
								}
								else
								{

									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "\"" << GeometryObjects.back().GetObjectType() << "\""
										 << " position must be defined by 3 zero parameters,";
									cout << " using the default (0,0,0) values." << endl;

									GeometryObjects.back().SetObjectPosition(0 * cm, 0 * cm, 0 * cm);
								}
							}

							if (position.size() == 3)
							{
								GeometryObjects.back().SetObjectPosition(position[0] * cm, position[1] * cm, position[2] * cm);
							}
							else
							{
								cout << " " << endl;
								cout << "--------------------------------- Warning ---------------------------------" << endl;
								cout << "\"" << GeometryObjects.back().GetObjectType() << "\""
									 << " position must be defined by 3 parameters,";
								cout << " using the default (0,0,0) values." << endl;

								GeometryObjects.back().SetObjectPosition(0 * cm, 0 * cm, 0 * cm);
							}

							break;
						}

						case 1:
						{ // rotation

							vector<G4double> rotation;

							while (std::getline(ss, x, ','))
							{
								rotation.push_back(stod(x) * deg);
							}

							if (rotation.size() == 3)
							{
								GeometryObjects.back().SetObjectRotation(rotation[0], rotation[1], rotation[2]);
							}
							else
							{
								cout << " " << endl;
								cout << "--------------------------------- Warning ---------------------------------" << endl;
								cout << "\"" << GeometryObjects.back().GetObjectType() << "\""
									 << " rotation must be defined by 3 parameters,";
								cout << " using the default (0,0,0) values." << endl;

								GeometryObjects.back().SetObjectRotation(0 * deg, 0 * deg, 0 * deg);
							}

							break;
						}

						case 2: // material

							GeometryObjects.back().ClearObjectMaterial();

							if (GeometryObjects.back().GetObjectType() == "Monitor")
							{
								cout << " " << endl;
								cout << "--------------------------------- Warning ---------------------------------" << endl;
								cout << "Monitor type object has as mandatory material HighVacuum, the materials provided for the monitor \""
									 << GeometryObjects.back().GetObjectName() << "\" will be ignored." << endl;

								GeometryObjects.back().AddObjectMaterial("HighVacuum");

								break;
							}

							while (std::getline(ss, x, ','))
							{
								GeometryObjects.back().AddObjectMaterial(x);
							}
							break;

						case 3: // mother volume
							ss >> x;

							if (GeometryObjects.back().GetObjectType() == "DiagDetector" || GeometryObjects.back().GetObjectType() == "CloverDetector")
							{
								cout << " " << endl;
								cout << "--------------------------------- Warning ---------------------------------" << endl;
								cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must have World as mother volume"
									 << ", mother volume changed to World" << endl;
								GeometryObjects.back().SetObjectMotherVolume("World");
								break;
							}

							GeometryObjects.back().SetObjectMotherVolume(x);
							break;

						case 4: // name

							ss >> x;

							GeometryObjects.back().SetObjectName(x);
							break;

						case 5:
						{ // shape parameters

							GeometryObjects.back().ClearShapeParameters();

							vector<G4double> shape;

							while (std::getline(ss, x, ','))
							{
								shape.push_back(stod(x));
							}

							if (GeometryObjects.back().GetObjectType() == "Collimator")
							{
								if (shape.size() != 4 || shape[0] == 0 || shape[1] == 0 || shape[2] == 0)
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 4 "
										 << "non-zero shape parameters, using default values (0,1,0,2pi,0,pi)" << endl;

									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(0.3 * cm);
								}
								else
								{
									for (unsigned int j = 0; j < shape.size(); j++)
									{
										GeometryObjects.back().AddShapeParameter(shape[j] * cm);
									}
								}
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "Sphere")
							{
								if (shape.size() != 6 || shape[1] == 0 || shape[0] > shape[1])
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 6 "
										 << "non-zero shape parameters, using default values (0,1,0,2pi,0,pi)" << endl;

									GeometryObjects.back().AddShapeParameter(0 * cm);
									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(0 * deg);
									GeometryObjects.back().AddShapeParameter(360 * deg);
									GeometryObjects.back().AddShapeParameter(0 * deg);
									GeometryObjects.back().AddShapeParameter(180 * deg);
								}
								else
								{

									for (unsigned int j = 0; j < shape.size(); j++)
									{

										if (j < 2)
										{
											GeometryObjects.back().AddShapeParameter(shape[j] * cm);
										}
										else
										{
											GeometryObjects.back().AddShapeParameter(shape[j] * deg);
										}
									}
								}
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "Box")
							{
								if (shape.size() != 3 || shape[0] == 0 || shape[1] == 0 || shape[2] == 0)
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 3 "
										 << "non-zero shape parameters, using default values (1,1,1)" << endl;

									for (unsigned int j = 0; j < 3; j++)
									{
										GeometryObjects.back().AddShapeParameter(1 * cm);
									}
								}
								else
								{
									for (unsigned int j = 0; j < shape.size(); j++)
									{
										GeometryObjects.back().AddShapeParameter(shape[j] * cm);
									}
								}
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "World")
							{
								if (shape.size() != 3 || shape[0] == 0 || shape[1] == 0 || shape[2] == 0)
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 3 "
										 << "non-zero shape parameters, using default values (2000,2000,2000)" << endl;

									for (unsigned int j = 0; j < 3; j++)
									{
										GeometryObjects.back().AddShapeParameter(2000 * cm);
									}
								}
								else
								{
									for (unsigned int j = 0; j < shape.size(); j++)
									{
										GeometryObjects.back().AddShapeParameter(shape[j] * cm);
									}
								}
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "Tube")
							{
								if (shape.size() != 5 || shape[0] >= shape[1] || shape[4] == 0)
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 5 "
										 << "constrained shape parameters, using default values (1,2,3,4,5)" << endl;

									for (unsigned int j = 1; j < 6; j++)
									{
										if (j < 4)
										{
											GeometryObjects.back().AddShapeParameter(j * cm);
										}
										else
										{
											GeometryObjects.back().AddShapeParameter(j * deg);
										}
									}
								}
								else
								{

									for (unsigned int j = 0; j < shape.size(); j++)
									{

										if (j < 4)
										{
											GeometryObjects.back().AddShapeParameter(shape[j] * cm);
										}
										else
										{
											GeometryObjects.back().AddShapeParameter(shape[j] * deg);
										}
									}
								}
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "Monitor")
							{

								if ((shape.size() != 2 && shape.size() != 3) || (shape[0] == 0 && shape.size() == 2) || shape[1] == 0)
								{

									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 2 or 3"
										 << " non-zero shape parameters , using default values (1,1)" << endl;

									for (unsigned int j = 0; j < 2; j++)
									{
										GeometryObjects.back().AddShapeParameter(1 * cm);
									}
									GeometryObjects.back().AddShapeParameter(0.00001 * cm);
									GeometryObjects.back().ClearObjectFeatures();
									GeometryObjects.back().AddObjectFeature("rectangle");
								}
								else
								{

									if (shape.size() == 3)
									{
										GeometryObjects.back().ClearObjectFeatures();
										GeometryObjects.back().AddObjectFeature("circular");

										for (unsigned int j = 0; j < shape.size(); j++)
										{
											if (j < 2)
											{
												GeometryObjects.back().AddShapeParameter(shape[j] * cm);
											}
											else
											{
												GeometryObjects.back().AddShapeParameter(shape[j] * deg);
											}
										}

										GeometryObjects.back().AddShapeParameter(0.00001 * cm);
									}

									if (shape.size() == 2)
									{

										GeometryObjects.back().ClearObjectFeatures();
										GeometryObjects.back().AddObjectFeature("rectangle");

										for (unsigned int j = 0; j < shape.size(); j++)
										{
											GeometryObjects.back().AddShapeParameter(shape[j] * cm);
										}

										GeometryObjects.back().AddShapeParameter(0.00001 * cm);
									}
								}
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "ELIADE")
							{

								if (shape.size() != 3 || shape[0] < 15)
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 3 "
										 << " shape parameters , using default value (20,90,135)" << endl;

									GeometryObjects.back().AddShapeParameter(20 * cm);
									GeometryObjects.back().AddShapeParameter(90 * deg);
									GeometryObjects.back().AddShapeParameter(135 * deg);
								}
								else
								{

									for (unsigned int j = 0; j < shape.size(); j++)
									{
										if (j == 0)
										{
											GeometryObjects.back().AddShapeParameter(shape[j] * cm);
										}
										else
										{
											GeometryObjects.back().AddShapeParameter(shape[j] * deg);
										}
									}
								}
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "PIXEL")
							{

								if (shape.size() != 3 || shape[0] == 0 || shape[1] == 0 || shape[2] == 0)
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 3 "
										 << "non-zero shape parameters, using default values (1,2,2)" << endl;

									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(2);
									GeometryObjects.back().AddShapeParameter(2);
								}
								else
								{

									GeometryObjects.back().AddShapeParameter(shape[0] * cm);
									GeometryObjects.back().AddShapeParameter(shape[1]);
									GeometryObjects.back().AddShapeParameter(shape[2]);
								}

								break;
							}

							if (GeometryObjects.back().GetObjectType() == "Mask")
							{

								if (shape.size() != 5 || shape[0] == 0 || shape[1] == 0 || shape[2] == 0 || shape[3] == 0 || shape[4] == 0)
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 5 "
										 << "non-zero shape parameters, using default values (1,1,1,0.1,1)" << endl;

									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(1 * cm);
									GeometryObjects.back().AddShapeParameter(0.1 * cm);
									GeometryObjects.back().AddShapeParameter(1);
								}
								else
								{

									GeometryObjects.back().AddShapeParameter(shape[0] * cm);
									GeometryObjects.back().AddShapeParameter(shape[1] * cm);
									GeometryObjects.back().AddShapeParameter(shape[2] * cm);
									GeometryObjects.back().AddShapeParameter(shape[3] * cm);
									GeometryObjects.back().AddShapeParameter(shape[4]);
								}

								break;
							}

							if (GeometryObjects.back().GetObjectType() == "VEGA")
							{
								if (shape.size() != 1 || !(1 <= shape[0] && shape[0] <= 7))
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType() << "\" must be defined by 1"
										 << "shape parameters between 1 and 7, using default value (1)" << endl;

									GeometryObjects.back().AddShapeParameter(1);
								}
								else
								{
									GeometryObjects.back().AddShapeParameter(shape[0]);
								}
								break;
							}

							break;
						}
						case 6: // active status

							ss >> x;

							if (GeometryObjects.back().GetObjectType() == "Monitor")
							{

								cout << " " << endl;
								cout << "--------------------------------- Warning ---------------------------------" << endl;
								cout << "Object of type \"" << GeometryObjects.back().GetObjectType()
									 << "\" is mandatory an active object, active value is ignored" << endl;

								GeometryObjects.back().SetObjectActiveStatus(1);
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "World")
							{

								cout << " " << endl;
								cout << "--------------------------------- Warning ---------------------------------" << endl;
								cout << "Object of type \"" << GeometryObjects.back().GetObjectType()
									 << "\" is mandatory an inactive object, active value is ignored" << endl;

								GeometryObjects.back().SetObjectActiveStatus(0);
								break;
							}

							if (GeometryObjects.back().GetObjectType() == "VEGA")
							{
								if (stod(x))
								{
									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType()
										 << "\" is mandatory an inactive object, active value is ignored" << endl;

									GeometryObjects.back().SetObjectActiveStatus(0);
								}
								break;
							}

							GeometryObjects.back().SetObjectActiveStatus(stod(x));
							break;

						case 7: // object overlap
							ss >> x;
							GeometryObjects.back().SetObjectOverlap(stod(x));
							break;

						case 8: // object histogram type

							GeometryObjects.back().ClearObjectHistogramType();

							while (std::getline(ss, x, ','))
							{

								if (x == "Sflux" && GeometryObjects.back().GetObjectType() != "Monitor")
								{

									cout << " " << endl;
									cout << "--------------------------------- Warning ---------------------------------" << endl;
									cout << "Object of type \"" << GeometryObjects.back().GetObjectType()
										 << "\" does not support Sflux histogram, flag is ignored" << endl;

									continue;
								}

								GeometryObjects.back().AddObjectHistogramType(x);
							}

							break;

						case 9: // object features

							GeometryObjects.back().ClearObjectFeatures();

							while (std::getline(ss, x, ','))
							{

								GeometryObjects.back().AddObjectFeature(x);
							}
							break;

						case 10:
						{ // moving parameters

							vector<G4double> shape;

							while (std::getline(ss, x, ','))
							{

								shape.push_back(stod(x));
							}

							if (shape.size() != 4)
							{

								cout << " " << endl;
								cout << "--------------------------------- Warning ---------------------------------" << endl;
								cout << "A moving object must be defined by four parameters, using the default non-rotating object" << endl;

								continue;
							}
							else
							{

								GeometryObjects.back().AddMovingParameters(shape[0]);
								GeometryObjects.back().AddMovingParameters(shape[1] * cm);
								GeometryObjects.back().AddMovingParameters(shape[2]);
								GeometryObjects.back().AddMovingParameters(shape[3] * cm);
							}

							break;
						}

						default:
							break;
						}
					}
				}
			}
		}
	}
	else
	{
		cout << "Unable to open InputObjects.txt" << endl;
		abort();
	}

	CheckGeometryObjects();
}

void AppsInput::CheckGeometryObjects()
{
	// check for world duplicates

	vector<G4int> worldCnt;

	for (unsigned int i = 0; i < GeometryObjects.size(); i++)
	{

		if (GeometryObjects[i].GetObjectType() == "World")
		{
			worldCnt.push_back(i);
		}
	}

	if (worldCnt.size() > 1)
	{

		cout << " " << endl;
		cout << "--------------------------------- Warning ---------------------------------" << endl;
		cout << "Multiple objects of type \"" << GeometryObjects[worldCnt[0]].GetObjectType() << "\" found all but the first will be ignored" << endl;

		for (unsigned int i = 1; i < worldCnt.size(); i++)
		{

			GeometryObjects.erase(GeometryObjects.begin() + worldCnt[i]);
			cout << " deleted object: " << worldCnt[i] << endl;
		}
	}

	// check for objects with undefined name and changes identical name to default "D" + type + count

	for (unsigned int i = 0; i < GeometryObjects.size(); i++)
	{

		if (GeometryObjects[i].GetObjectName() == "undefined")
		{
			DefaultNameCount++;
			GeometryObjects[i].SetObjectName("D" + GeometryObjects[i].GetObjectType() + to_string(DefaultNameCount));
		}
	}

	// check for duplicate names and changes duplicates to  default "D" + type + count

	for (unsigned int i = 0; i < GeometryObjects.size(); i++)
	{

		for (unsigned int o = i + 1; o < GeometryObjects.size(); o++)
		{

			if (GeometryObjects[i].GetObjectName() == GeometryObjects[o].GetObjectName())
			{

				DefaultNameCount++;
				GeometryObjects[o].SetObjectName("D" + GeometryObjects[o].GetObjectType() + to_string(DefaultNameCount));

				cout << "" << endl;
				cout << "--------------------------------- Warning ---------------------------------" << endl;
				cout << "Duplicate name found, changed the name of object " << o << " to the default name \"" << GeometryObjects[o].GetObjectName() << "\"" << endl;
			}
		}
	}

	// check if all the objects have the shape parameters correct, if there is a missmatch between shape parameters and object type default parameters are applied

	for (unsigned int i = 0; i < GeometryObjects.size(); i++)
	{

		if ((GeometryObjects[i].GetObjectType() == "Box" || GeometryObjects[i].GetObjectType() == "World") && GeometryObjects[i].GetShapeParameters().size() != 3)
		{
			cout << " " << endl;
			cout << "--------------------------------- Warning ---------------------------------" << endl;
			cout << "Object of type \"" << GeometryObjects[i].GetObjectType() << "\" must be defined by 3 "
				 << "non-zero shape parameters, using default values (1,1,1)" << endl;

			for (unsigned int j = 0; j < 3; j++)
			{
				GeometryObjects[i].AddShapeParameter(1 * cm);
			}
		}

		if (GeometryObjects[i].GetObjectType() == "Tube" && GeometryObjects[i].GetShapeParameters().size() != 5)
		{
			cout << " " << endl;
			cout << "--------------------------------- Warning ---------------------------------" << endl;
			cout << "Object of type \"" << GeometryObjects[i].GetObjectType() << "\" must be defined by 5 "
				 << "constrained shape parameters, using default values (1,2,3,4,5)" << endl;

			for (unsigned int j = 1; j < 6; j++)
			{
				if (j < 4)
				{
					GeometryObjects[i].AddShapeParameter(j * cm);
				}
				else
				{
					GeometryObjects[i].AddShapeParameter(j * deg);
				}
			}
		}

		if (GeometryObjects[i].GetObjectType() == "Monitor" && GeometryObjects[i].GetShapeParameters().size() != 3)
		{
			if (GeometryObjects[i].GetShapeParameters().size() != 4 || GeometryObjects[i].GetShapeParameters()[0] == GeometryObjects[i].GetShapeParameters()[1])
			{
				cout << " " << endl;
				cout << "--------------------------------- Warning ---------------------------------" << endl;
				cout << "Object of type \"" << GeometryObjects[i].GetObjectType() << "\" must be defined by 2 or 3 "
					 << "non-zero shape parameters, using default values (1,1)" << endl;
				GeometryObjects[i].ClearShapeParameters();

				for (unsigned int j = 0; j < 2; j++)
				{
					GeometryObjects[i].AddShapeParameter(1 * cm);
				}

				GeometryObjects[i].AddShapeParameter(0.00001 * cm);
				GeometryObjects[i].ClearObjectFeatures();
				GeometryObjects[i].AddObjectFeature("rectangle");
			}
		}
	}
}

void AppsInput::ReadRunParameters()
{

	G4String z;
	ifstream inputfile("SimParameters.txt");

	if (inputfile.is_open())
	{

		while (!inputfile.eof())
		{

			getline(inputfile, z);
			for (unsigned int i = 0; i < InputParametersOption.size(); i++)
			{
				std::size_t found = z.find(InputParametersOption[i]);

				if (found != std::string::npos)
				{

					std::size_t firstpos = z.find("[") + 1;
					std::size_t secondpos = z.find("]");

					stringstream ss;

					ss.str(z.substr(firstpos, (secondpos - firstpos)));

					string x;

					switch (i)
					{

					case 0: // "display interface"
						ss >> x;
						SetInterfaceMode(stod(x));
						break;

					case 1: // "number of runs"
						ss >> x;
						SetNumberOfRuns(stod(x));
						break;

					case 2: // "number of events"

						NumberOfEvents.clear();

						while (std::getline(ss, x, ','))
						{
							AddNumberOfEvents(stod(x));
						}

						break;

					case 3: // "gamma energy"

						GammaEnergy.clear();

						while (std::getline(ss, x, ','))
						{
							AddGammaEnergy(stod(x));
						}

						break;

					case 4: // "bandwidth"

						GammaBandwidth.clear();

						while (std::getline(ss, x, ','))
						{
							AddGammaBandwidth(stod(x));
						}

						break;

					case 5: // "energy increment"

						ss >> x;
						SetEnergyIncrement(stod(x));
						break;

					case 6:
					{ // "polar angle"

						PolarAngle.clear();

						vector<G4double> angles;

						while (std::getline(ss, x, ','))
						{
							if (x == "pi")
							{
								x = "3.14159";
							}
							angles.push_back(stod(x));
						}

						if (angles.size() == 2)
						{
							SetPolarAngle(angles[0] * rad, angles[1] * rad);
						}
						else
						{
							cout << " " << endl;
							cout << "--------------------------------- Warning ---------------------------------" << endl;
							cout << "Polar angle must be defined by 2 parameters,";
							cout << " using default (0,pi) values." << endl;

							SetPolarAngle(0 * rad, 0 * rad);
						}

						break;
					}

					case 7:
					{ // "azimuthal angle"

						AzimuthalAngle.clear();

						vector<G4double> angles;

						while (std::getline(ss, x, ','))
						{
							if (x == "2pi")
							{
								x = "6.28318";
							}
							if (x == "pi")
							{
								x = "3.14159";
							}
							angles.push_back(stod(x));
						}

						if (angles.size() == 2)
						{
							SetAzimuthalAngle(angles[0], angles[1]);
						}
						else
						{
							cout << " " << endl;
							cout << "--------------------------------- Warning ---------------------------------" << endl;
							cout << "Azimuthal angle must be defined by 2 parameters,";
							cout << " using default (0,2pi) values." << endl;

							SetAzimuthalAngle(0, 0);
						}

						break;
					}

					case 8:
					{ // "source position"

						vector<G4double> position;

						while (std::getline(ss, x, ','))
						{
							position.push_back(stod(x));
						}

						if (position.size() == 3)
						{
							SetSourcePosition(position[0] * cm, position[1] * cm, position[2] * cm);
						}
						else
						{
							cout << " " << endl;
							cout << "--------------------------------- Warning ---------------------------------" << endl;
							cout << "Source position must be defined by 3 parameters,";
							cout << " using the default (0,0,0) values." << endl;

							SetSourcePosition(0 * cm, 0 * cm, 0 * cm);
						}

						break;
					}

					case 9: // time structure

						ss >> x;
						SetBeamType(x);
						if (x == "PaNi")
						{
							ReadNeutronSpectrum();
						}
						break;

					case 10: // number of particles per bunch

						ss >> x;
						SetBunchParticles(stod(x));
						break;

					case 11: // the activity of the decay source

						ClearActivityParameters();
						while (std::getline(ss, x, ','))
						{
							AddActivityParameter(stod(x));
						}
						break;

					case 12: // the type of nucleus for the decay source

						ss >> x;
						SetDecayNucleus(x);
						break;

					case 13: // the type of nucleus for the decay source

						ss >> x;
						if (x == "all")
						{
							SetNumberOfCores(G4Threading::G4GetNumberOfCores());
						}
						else
						{
							SetNumberOfCores(stod(x));
						}
						break;

					case 14: // the time structure of the beam

						ClearTimeSTructureParameters();

						while (std::getline(ss, x, ','))
						{
							AddTimeStructureParameter(stod(x));
						}
						break;

					case 15: // beam output option

						ClearBeamHistogramType();

						while (std::getline(ss, x, ','))
						{
							AddBeamHistogramType(x);
						}
						break;

					case 16: // beam shape parameters option

						ClearShapeParameters();

						while (std::getline(ss, x, ','))
						{
							AddShapeParameter(stod(x) * cm);
						}
						break;

					case 17: // beam polarization

						vector<G4double> polar;

						while (std::getline(ss, x, ','))
						{
							polar.push_back(stod(x));
						}

						if (polar.size() == 3)
						{
							SetPolarization(polar[0], polar[1], polar[2]);
						}
						else
						{
							cout << " " << endl;
							cout << "--------------------------------- Warning ---------------------------------" << endl;
							cout << "Polarization must be defined by 3 parameters,";
							cout << " using the default (0,0,0) values." << endl;

							SetPolarization(0, 0, 0);
						}
						break;
					}
				}
			}
		}
	}
	else
	{

		cout << "Unable to open RunParameters.txt" << endl;
		abort();
	}
}

void AppsInput::ReadNeutronSpectrum()
{

	G4String z;
	ifstream inputfile("PuBe_neutron_spectrum.dat");

	neutron_spectrum = new TH1D("PuBe_neuttron_spectrum", "PuBe_neuttron_spectrum", 140, 0., 14000.);

	if (inputfile.is_open())
	{

		while (!inputfile.eof())
		{
			getline(inputfile, z);
			
			stringstream ss;

			ss.str(z.substr(1,13));

			string x;
			
			ss>>x;
			
			ss.str(z.substr(14,29));
			
			string y;
			
			ss>>y;
						
			neutron_spectrum->Fill(stod(x),stod(y));

		}
	}
	else
	{

		cout << "Unable to open PuBe_neutron_spectrum.dat" << endl;
		abort();
	}
}


void AppsInput::InitializeParameters()
{

	BunchParticles = 260000;							 // number of particles per bunch
	ActivityParameters.push_back(100000);				 // source activity in Bq
	ActivityParameters.push_back(100);					 // source activity uncertainity
	BeamType = "source";								 // type of source can be GBS or calibration source
	SetInterfaceMode(false);							 // keeps display info
	SetNumberOfRuns(1);									 // keeps the number of runs
	AddNumberOfEvents(10000);							 // keeps the number of events
	AddGammaEnergy(1);									 // energy of the gamma beam
	AddGammaBandwidth(0);								 // FWHM of the gamma beam
	SetEnergyIncrement(0);								 // increment of the energy in multiple energy runs
	SetPolarAngle(0, CLHEP::pi);						 // min angle, max angle
	SetAzimuthalAngle(0, 2 * CLHEP::pi);				 // min angle
	SetSourcePosition(0, 0, 0);							 // source in position [0,0,0]
	SetNumberOfCores(G4Threading::G4GetNumberOfCores()); // max number of cores
	TimeStructure.push_back(32);
	TimeStructure.push_back(16 * ns);
	TimeStructure.push_back(10 * ms);
	ShapeParameters.push_back(1);
	ShapeParameters.push_back(1);
	Polarization = {0, 0, 0};
}

void AppsInput::DisplayMenu()
{
}

void AppsInput::PrintGeometryParameters()
{
	/*
	cout << " " << endl;
	cout << "Geometry objects: " << endl; 
	cout << "--------------------------------" << endl; 	
	for(unsigned int i = 0; i < GeometryObjects.size(); i++){
	
			
	
	}

*/
}

void AppsInput::PrintSimParameters()
{

	cout << " " << endl;
	cout << "Simulation parameters: " << endl;
	cout << "--------------------------------" << endl;
	cout << "User interface: - display interface = [" << (int)GetInterfaceMode() << "]" << endl;
	cout << " " << endl;
	cout << "Run parameters: - number of runs = [" << GetNumberOfRuns() << "]" << endl;
	cout << "                - number of events = [";
	for (unsigned int i = 0; i < GetNumberOfEvents().size(); i++)
	{

		if (i == (GetNumberOfEvents().size() - 1))
		{
			cout << GetNumberOfEvents()[i] << "]" << endl;
		}
		else
		{
			cout << GetNumberOfEvents()[i] << ", ";
		}
	}
	cout << " " << endl;
	cout << "Beam parameters (MeV): - gamma energy = [";
	for (unsigned int i = 0; i < GetGammaEnergy().size(); i++)
	{

		if (i == (GetGammaEnergy().size() - 1))
		{
			cout << GetGammaEnergy()[i] << "]" << endl;
		}
		else
		{
			cout << GetGammaEnergy()[i] << ", ";
		}
	}
	cout << "                       - bandwidth = [";
	for (unsigned int i = 0; i < GetGammaBandwidth().size(); i++)
	{

		if (i == (GetGammaBandwidth().size() - 1))
		{
			cout << GetGammaBandwidth()[i] << "]" << endl;
		}
		else
		{
			cout << GetGammaBandwidth()[i] << ", ";
		}
	}
	cout << "		       - energy increment = [" << GetEnergyIncrement() << "]" << endl;
	cout << " 		       - polar angle (theta min-max) = [";
	for (unsigned int i = 0; i < GetPolarAngle().size(); i++)
	{

		if (i == (GetPolarAngle().size() - 1))
		{
			cout << GetPolarAngle()[i] << "]" << endl;
		}
		else
		{
			cout << GetPolarAngle()[i] << ", ";
		}
	}
	cout << " 		       - azimuthal angle (phi min-max) = [";
	for (unsigned int i = 0; i < GetAzimuthalAngle().size(); i++)
	{

		if (i == (GetAzimuthalAngle().size() - 1))
		{
			cout << GetAzimuthalAngle()[i] << "]" << endl;
		}
		else
		{
			cout << GetAzimuthalAngle()[i] << ", ";
		}
	}
	cout << "		       - source position = [";
	cout << GetSourcePosition().getX() << ", ";
	cout << GetSourcePosition().getY() << ", ";
	cout << GetSourcePosition().getZ() << "]" << endl;
	cout << " " << endl;
}

AppsInput *AppsInput::Instance()
{

	if (fgInstance == 0)
	{
		fgInstance = new AppsInput();
	}

	return fgInstance;
}


