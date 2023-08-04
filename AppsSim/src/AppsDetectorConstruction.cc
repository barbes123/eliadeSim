#include "AppsDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "DiagDetector.hh"
#include "ScintDetector.hh"
#include "CeBrDetector.hh"
#include "AppsInput.hh"
#include "AppsSD.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"
#include "CloverDetector.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "AppsAnalysisObject.hh"
#include "AppsInput.hh"
#include "AppsRun.hh"
#include "G4MultiUnion.hh"
#include "vegaCollimator.hh"

AppsDetectorConstruction::AppsDetectorConstruction()
	: G4VUserDetectorConstruction()
{
}

AppsDetectorConstruction::~AppsDetectorConstruction()
{
}

void AppsDetectorConstruction::ConstructWorld()
{

	AppsInput *FInput = AppsInput::Instance();

	G4int inputWorldId = -1;

	for (unsigned int i = 0; i < FInput->GetGeometryObjects().size(); i++)
	{

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "World")
		{

			inputWorldId = i;
			break;
		}
	}

	if (inputWorldId >= 0)
	{

		ConstructWorldObject(FInput->GetGeometryObjects()[inputWorldId]); // build the object defined in the input file
	}
	else
	{

		AppsGeometryObject DefaultWorldObject;
		DefaultWorldObject.ClearObjectMaterial();
		DefaultWorldObject.AddObjectMaterial("AIR"); // default intialization material is AIR
		DefaultWorldObject.SetObjectName("World");

		DefaultWorldObject.ClearShapeParameters();

		for (unsigned int j = 0; j < 3; j++)
		{

			DefaultWorldObject.AddShapeParameter(80000);
		}

		ConstructWorldObject(DefaultWorldObject); // build default world objects
	}
}

void AppsDetectorConstruction::ConstructWorldObject(AppsGeometryObject object)
{

	world_mat = BuildMaterial(object.GetObjectMaterial());

	G4Box *solidWorld = new G4Box(object.GetObjectName(), object.GetShapeParameters()[0], object.GetShapeParameters()[1], object.GetShapeParameters()[2]);

	logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");

	LogicalVolumes.push_back(logicWorld);

	physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, object.GetObjectOverlap());
}

void AppsDetectorConstruction::ConstructPaNiSphere()
{

	AppsInput *FInput = AppsInput::Instance();

	G4Material *PC = BuildMaterial({"POLYCARBONATE"});
	G4Material *PaNi1 = BuildMaterial({"PaNi1"});
	G4Material *PaNi2 = BuildMaterial({"PaNi2"});
	G4Material *PaNi3 = BuildMaterial({"PaNi3"});
	
	G4double Rout = 8.*cm;
	G4double TShell = 0.175*cm;
	G4double Rhole = 1.5*cm;
	G4double Rcyltop = 4.*cm;
	G4double Hcyltop = 10.*cm;
	G4double Rcylfill = 1.*cm;
	G4double Hcylfill = 11.*cm;
	G4double Hcylfilled = Hcylfill - 4.*cm;
	G4double Xcylfill = Rcyltop + Rcylfill + 0.5*cm;

	//Shell
	G4VSolid *ShellSphere = new G4Sphere("ShellSphere", Rout-TShell, Rout, 0., twopi, 0., pi);
	
	G4VSolid *CylTop = new G4Tubs("CylTop", Rhole, Rcyltop, Hcyltop/2., 0., twopi);
	G4VSolid *CylTopVoid = new G4Tubs("CylTopVoid", Rhole+TShell, Rcyltop-TShell, Hcyltop/2.-TShell, 0., twopi);
	
	G4VSolid *CylFill = new G4Tubs("CylTop", Rcylfill-TShell, Rcylfill, Hcylfill/2., 0., twopi);
	
	G4RotationMatrix *yRot1 = new G4RotationMatrix;
	yRot1->rotateY(0 * deg);	
	
	G4VSolid *Shell0 = new G4UnionSolid("Shell0", ShellSphere, CylTop, yRot1, G4ThreeVector(0., 0., Hcyltop/2.));
	G4VSolid *Shell00 = new G4SubtractionSolid("Shell00", Shell0, CylTopVoid, yRot1, G4ThreeVector(0., 0., Hcyltop/2.));
	G4VSolid *Shell000 = new G4UnionSolid("Shell000", Shell00, CylFill, yRot1, G4ThreeVector(Xcylfill, 0., Hcylfill/2.));

	//Inner Sphere
	G4VSolid *Sphere = new G4Sphere("Sphere", 0., Rout-TShell, 0., twopi, 0., pi);
	
	G4VSolid *Hole = new G4Tubs("Hole", 0., Rhole, Hcyltop/2., 0., twopi);
	
	G4VSolid *Filled = new G4Tubs("Filled", 0., Rcylfill-TShell, Hcylfilled/2., 0., twopi);
	
	G4VSolid *PaNiSphere0 = new G4SubtractionSolid("PaNiSphere0", Sphere, Hole, yRot1, G4ThreeVector(0., 0., Hcyltop/2.));
	G4VSolid *PaNiSphere = new G4UnionSolid("PaNiSphere", PaNiSphere0, Filled, yRot1, G4ThreeVector(Xcylfill, 0., Hcylfilled/2.));
	
	G4VSolid *Shell = new G4SubtractionSolid("Shell", Shell000, PaNiSphere, yRot1, G4ThreeVector(0., 0., 0.));

	G4VSolid *CutBox = new G4Box("CutBox", Rout, Rout, Rout);

	G4double ZCut1 = 5.5*cm;
	G4double ZCut2 = -5.5*cm;
	
	G4VSolid *PaNiSphere1 = new G4SubtractionSolid("PaNiSphere1", PaNiSphere, CutBox, yRot1, G4ThreeVector(0., 0., ZCut1 - Rout));
	G4VSolid *PaNiSphere2_0 = new G4SubtractionSolid("PaNiSphere2_0", PaNiSphere, CutBox, yRot1, G4ThreeVector(0., 0., ZCut1 + Rout));
	G4VSolid *PaNiSphere2 = new G4SubtractionSolid("PaNiSphere2", PaNiSphere2_0, CutBox, yRot1, G4ThreeVector(0., 0., ZCut2 - Rout));
	G4VSolid *PaNiSphere3 = new G4SubtractionSolid("PaNiSphere3", PaNiSphere, CutBox, yRot1, G4ThreeVector(0., 0., ZCut2 + Rout));

	G4LogicalVolume *logicPaNiShell = new G4LogicalVolume(Shell, PC, "logicPaNiShell");
	G4VisAttributes* vis = new G4VisAttributes(false);
	logicPaNiShell->SetVisAttributes(vis);
	G4LogicalVolume *logicPaNiSphere1 = new G4LogicalVolume(PaNiSphere1, PaNi1, "logicPaNiSphere1");
	G4LogicalVolume *logicPaNiSphere2 = new G4LogicalVolume(PaNiSphere2, PaNi2, "logicPaNiSphere2");
	G4LogicalVolume *logicPaNiSphere3 = new G4LogicalVolume(PaNiSphere3, PaNi3, "logicPaNiSphere3");

	LogicalVolumes.push_back(logicPaNiShell);
	LogicalVolumes.push_back(logicPaNiSphere1);
	LogicalVolumes.push_back(logicPaNiSphere2);
	LogicalVolumes.push_back(logicPaNiSphere3);

	G4PVPlacement *shellPlace = new G4PVPlacement(yRot1, FInput->GetSourcePosition(), logicPaNiShell, "logicPaNiShell", logicWorld, false, 0, true);
	G4PVPlacement *PaNiSphere1Place = new G4PVPlacement(yRot1, FInput->GetSourcePosition(), logicPaNiSphere1, "logicPaNiSphere1", logicWorld, false, 0, true);
	G4PVPlacement *PaNiSphere2Place = new G4PVPlacement(yRot1, FInput->GetSourcePosition(), logicPaNiSphere2, "logicPaNiSphere2", logicWorld, false, 0, true);
	G4PVPlacement *PaNiSphere3Place = new G4PVPlacement(yRot1, FInput->GetSourcePosition(), logicPaNiSphere3, "logicPaNiSphere3", logicWorld, false, 0, true);

	if (shellPlace->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"Shell\"" << endl;
		cout << " " << endl;
	}
	if (PaNiSphere1Place->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"PaNiSphere1\"" << endl;
		cout << " " << endl;
	}
	if (PaNiSphere2Place->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"PaNiSphere2\"" << endl;
		cout << " " << endl;
	}
	if (PaNiSphere3Place->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"PaNiSphere3\"" << endl;
		cout << " " << endl;
	}

	/*FInput->GeometryObjects[objectId].AddObjectActiveVolume(FInput->GetGeometryObjects()[objectId].GetObjectName());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);*/
}

void AppsDetectorConstruction::ConstructSphereObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4Material *material = BuildMaterial(object.GetObjectMaterial());

	G4VSolid *Sphere = new G4Sphere(object.GetObjectName(), object.GetShapeParameters()[0], object.GetShapeParameters()[1], object.GetShapeParameters()[2],
									object.GetShapeParameters()[3], object.GetShapeParameters()[4], object.GetShapeParameters()[5]);

	G4LogicalVolume *logicSphere = new G4LogicalVolume(Sphere, material, object.GetObjectName());

	LogicalVolumes.push_back(logicSphere);

	G4PVPlacement *boxPlace = new G4PVPlacement(object.GetObjectRotation(), object.GetObjectPosition(), logicSphere, object.GetObjectName(),
												GetLogicalVolume(object.GetObjectMotherVolume()), false, 0, object.GetObjectOverlap());

	if (boxPlace->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"" << object.GetObjectName() << "\"" << endl;
		cout << " " << endl;
	}

	FInput->GeometryObjects[objectId].AddObjectActiveVolume(FInput->GetGeometryObjects()[objectId].GetObjectName());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::ConstructCollimatorObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4Material *material = BuildMaterial(object.GetObjectMaterial());

	G4VSolid *Box = new G4Box("Box", object.GetShapeParameters()[0], object.GetShapeParameters()[1], object.GetShapeParameters()[2]);

	G4VSolid *Tubs = new G4Tubs("Tubs", 0, object.GetShapeParameters()[3], object.GetShapeParameters()[2] * 1.1, 0 * deg, 360 * deg);

	G4VSolid *Collimator = new G4SubtractionSolid(object.GetObjectName(), Box, Tubs, 0, G4ThreeVector(0, 0, 0));

	G4LogicalVolume *logicBox = new G4LogicalVolume(Collimator, material, object.GetObjectName());

	LogicalVolumes.push_back(logicBox);

	bool isVisible = true;

	double Red = 102;
	double Green = 102;
	double Blue = 102;

	G4VisAttributes *tempVisAttr = new G4VisAttributes(isVisible,
													   G4Colour(Red / 255., Green / 255., Blue / 255., 0.7));

	logicBox->SetVisAttributes(tempVisAttr);

	G4PVPlacement *boxPlace = new G4PVPlacement(object.GetObjectRotation(), object.GetObjectPosition(), logicBox, object.GetObjectName(),
												GetLogicalVolume(object.GetObjectMotherVolume()), false, 0, object.GetObjectOverlap());

	if (boxPlace->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"" << object.GetObjectName() << "\"" << endl;
		cout << " " << endl;
	}

	FInput->GeometryObjects[objectId].AddObjectActiveVolume(FInput->GetGeometryObjects()[objectId].GetObjectName());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::ConstructBoxObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4Material *material = BuildMaterial(object.GetObjectMaterial());

	G4VSolid *Box = new G4Box(object.GetObjectName(), object.GetShapeParameters()[0], object.GetShapeParameters()[1], object.GetShapeParameters()[2]);
	G4LogicalVolume *logicBox = new G4LogicalVolume(Box, material, object.GetObjectName());

	LogicalVolumes.push_back(logicBox);

	bool isVisible = true;

	double Red = 0;
	double Green = 0;
	double Blue = 0;

	G4PVPlacement *boxPlace = new G4PVPlacement(object.GetObjectRotation(), object.GetObjectPosition(), logicBox, object.GetObjectName(),
												GetLogicalVolume(object.GetObjectMotherVolume()), false, 0, object.GetObjectOverlap());

	PlacementVolumes.push_back(boxPlace);

	if (boxPlace->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"" << object.GetObjectName() << "\"" << endl;
		cout << " " << endl;
	}

	FInput->GeometryObjects[objectId].AddObjectActiveVolume(FInput->GetGeometryObjects()[objectId].GetObjectName());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());

	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::ConstructMaskObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4Material *material = BuildMaterial(object.GetObjectMaterial());

	G4VSolid *Mask = new G4Box(object.GetObjectName(), object.GetShapeParameters()[0], object.GetShapeParameters()[1], object.GetShapeParameters()[2]);
	G4VSolid *Tubs = new G4Tubs("Tubs", 0, object.GetShapeParameters()[3], object.GetShapeParameters()[2] * 1.1, 0 * deg, 360 * deg);

	//G4VSolid * Mask;

	std::random_device rd;	//Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(-object.GetShapeParameters()[0] + object.GetShapeParameters()[3], object.GetShapeParameters()[0] - object.GetShapeParameters()[3]);

	std::random_device rd1;	  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen1(rd1()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis1(-object.GetShapeParameters()[1] + object.GetShapeParameters()[3], object.GetShapeParameters()[1] - object.GetShapeParameters()[3]);

	for (unsigned int i = 0; i < object.GetShapeParameters()[4]; i++)
	{

		double xrand = dis(gen);
		double yrand = dis1(gen1);

		Mask = new G4SubtractionSolid(object.GetObjectName(), Mask, Tubs, 0, G4ThreeVector(xrand, yrand, 0));
	}

	G4LogicalVolume *logicBox = new G4LogicalVolume(Mask, material, object.GetObjectName());

	LogicalVolumes.push_back(logicBox);

	G4PVPlacement *boxPlace = new G4PVPlacement(object.GetObjectRotation(), object.GetObjectPosition(), logicBox, object.GetObjectName(),
												GetLogicalVolume(object.GetObjectMotherVolume()), false, 0, object.GetObjectOverlap());

	PlacementVolumes.push_back(boxPlace);

	if (boxPlace->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"" << object.GetObjectName() << "\"" << endl;
		cout << " " << endl;
	}

	FInput->GeometryObjects[objectId].AddObjectActiveVolume(FInput->GetGeometryObjects()[objectId].GetObjectName());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());

	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::ConstructScintDetectorObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4ThreeVector position = object.GetObjectPosition();

	ScintDetector *detector = new ScintDetector(objectId);

	G4AssemblyVolume *assembly = detector->GetDetectorAssembly();

	assembly->MakeImprint(logicWorld, position, object.GetObjectRotation(),
						  0, object.GetObjectOverlap());

	auto iterator = assembly->GetVolumesIterator();

	for (unsigned int it = 0; it < assembly->TotalImprintedVolumes(); it++)
	{

		if ((*iterator)->CheckOverlaps())
		{

			cout << " " << endl;
			cout << "------------------- Warning -------------------" << endl;
			cout << "Overlap detected for volume \"" << (*iterator)->GetName() << "\" part of assembly \""
				 << object.GetObjectName() << "\"" << endl;
			cout << " " << endl;
		}

		iterator++;
	}

	FInput->GeometryObjects[objectId].SetObjectActiveVolumes(detector->GetActiveVolumes());
	FInput->GeometryObjects[objectId].SetObjectFWHMParam(detector->GetFWHMParameters());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetDetectorResolution(detector->GetFWHMParameters());

	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());

	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::ConstructCeBrDetectorObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4ThreeVector position = object.GetObjectPosition();

	CeBrDetector *detector = new CeBrDetector(objectId);

	G4AssemblyVolume *assembly = detector->GetDetectorAssembly();

	assembly->MakeImprint(logicWorld, position, object.GetObjectRotation(),
						  0, object.GetObjectOverlap());

	auto iterator = assembly->GetVolumesIterator();

	for (unsigned int it = 0; it < assembly->TotalImprintedVolumes(); it++)
	{

		if ((*iterator)->CheckOverlaps())
		{

			cout << " " << endl;
			cout << "------------------- Warning -------------------" << endl;
			cout << "Overlap detected for volume \"" << (*iterator)->GetName() << "\" part of assembly \""
				 << object.GetObjectName() << "\"" << endl;
			cout << " " << endl;
		}

		iterator++;
	}

	FInput->GeometryObjects[objectId].SetObjectActiveVolumes(detector->GetActiveVolumes());
	FInput->GeometryObjects[objectId].SetObjectFWHMParam(detector->GetFWHMParameters());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetDetectorResolution(detector->GetFWHMParameters());

	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());

	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::ConstructDiagDetectorObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4ThreeVector position = object.GetObjectPosition();

	DiagDetector *detector = new DiagDetector(objectId);

	G4AssemblyVolume *assembly = detector->GetDetectorAssembly();

	assembly->MakeImprint(logicWorld, position, object.GetObjectRotation(),
						  0, object.GetObjectOverlap());

	auto iterator = assembly->GetVolumesIterator();

	for (unsigned int it = 0; it < assembly->TotalImprintedVolumes(); it++)
	{

		if ((*iterator)->CheckOverlaps())
		{

			cout << " " << endl;
			cout << "------------------- Warning -------------------" << endl;
			cout << "Overlap detected for volume \"" << (*iterator)->GetName() << "\" part of assembly \""
				 << object.GetObjectName() << "\"" << endl;
			cout << " " << endl;
		}

		iterator++;
	}

	FInput->GeometryObjects[objectId].SetObjectActiveVolumes(detector->GetActiveVolumes());
	FInput->GeometryObjects[objectId].SetObjectFWHMParam(detector->GetFWHMParameters());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetDetectorResolution(detector->GetFWHMParameters());
	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());

	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::ConstructPIXELObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	vector<G4double> shapePar = object.GetShapeParameters();

	double xLength = shapePar[0] * shapePar[1];
	double yLength = shapePar[0] * shapePar[2];

	G4AssemblyVolume *pixelAssembly = new G4AssemblyVolume();

	for (unsigned int i = 0; i < shapePar[1]; i++)
	{

		double posX = -xLength / 2.0 + shapePar[0] / 2.0 + i * shapePar[0];

		for (unsigned int o = 0; o < shapePar[2]; o++)
		{

			double posY = -yLength / 2.0 + shapePar[0] / 2.0 + o * shapePar[0];

			string pixelName = object.GetObjectName() + "_pixel_" + to_string(i) + "_" + to_string(o);

			G4ThreeVector pos = G4ThreeVector(posX, posY, 0);

			G4RotationMatrix *yRot2 = new G4RotationMatrix;
			yRot2->rotateY(0 * deg);

			G4Material *material = BuildMaterial(object.GetObjectMaterial());

			G4VSolid *Box = new G4Box(object.GetObjectName(), shapePar[0] / 2.0 - shapePar[0] * 0.01 * mm, shapePar[0] / 2.0 - shapePar[0] * 0.01 * mm, shapePar[0] * 2);

			G4LogicalVolume *logicBox = new G4LogicalVolume(Box, material, pixelName);

			pixelAssembly->AddPlacedVolume(logicBox, pos, yRot2);

			vector<AppsAnalysisObject *> temp;

			AppsAnalysisObject *analysis = new AppsAnalysisObject();
			analysis->SetObjectName(pixelName);
			analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
			temp.push_back(analysis);

			FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);

			FInput->GeometryObjects[objectId].AddObjectActiveVolume(pixelName);
		}
	}

	G4ThreeVector position = object.GetObjectPosition();
	G4RotationMatrix rotation = *object.GetObjectRotation();

	pixelAssembly->MakeImprint(logicWorld, position, &rotation,
							   0, object.GetObjectOverlap());
}

void AppsDetectorConstruction::ConstructVEGACollimator(G4int objectId)
{
	AppsInput *FInput = AppsInput::Instance();
	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	auto collimator = new vegaCollimator(objectId);

	G4AssemblyVolume *assembly = collimator->GetAssembly();

	G4ThreeVector position = FInput->GetSourcePosition();

	assembly->MakeImprint(logicWorld, position, 0, 0, object.GetObjectOverlap());

	auto iterator = assembly->GetVolumesIterator();

	for (unsigned int it = 0; it < assembly->TotalImprintedVolumes(); it++)
	{

		if ((*iterator)->CheckOverlaps())
		{

			cout << " " << endl;
			cout << "------------------- Warning -------------------" << endl;
			cout << "Overlap detected for volume \"" << (*iterator)->GetName() << "\" part of assembly \""
				 << object.GetObjectName() << "\"" << endl;
			cout << " " << endl;
		}

		iterator++;
	}
}

void AppsDetectorConstruction::ConstructELIADEObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	vector<G4double> shapePar = object.GetShapeParameters();

	G4AssemblyVolume *ELIADEassembly = new G4AssemblyVolume();

	G4ThreeVector position = object.GetObjectPosition();
	G4RotationMatrix rotation = *object.GetObjectRotation();

	G4bool isSegmented = false;

	for (unsigned int o = 0; o < object.GetObjectFeatures().size(); o++)
	{
		if (object.GetObjectFeatures()[o] == "segmented")
		{
			isSegmented = true;
		}
	}

	// obtain digital signature and add is as digital object

	G4bool isDigital = false;

	for (unsigned int i = 0; i < object.GetObjectHistogramType().size(); i++)
	{
		if (object.GetObjectHistogramType()[i] == "digital")
		{
			isDigital = true;
		}
	}

	CloverDetector *detector = new CloverDetector(objectId);

	for (unsigned int i = 0; i < 8; i++)
	{

		G4ThreeVector cPosition;
		G4RotationMatrix cRotation;

		if (i == 0)
		{
			cRotation.rotateY(shapePar[1]);
			cPosition = G4ThreeVector(shapePar[0] * sin(shapePar[1]), 0, shapePar[0] * cos(shapePar[1]));
		}

		if (i == 1)
		{
			cRotation.rotateY(-shapePar[1]);
			cPosition = G4ThreeVector(-shapePar[0] * sin(shapePar[1]), 0, shapePar[0] * cos(shapePar[1]));
		}

		if (i == 2)
		{
			cPosition = G4ThreeVector(0, shapePar[0] * sin(shapePar[1]), shapePar[0] * cos(shapePar[1]));
			cRotation.rotateY(-shapePar[1]);
			cRotation.rotateZ(-CLHEP::pi / 2);
		}

		if (i == 3)
		{
			cPosition = G4ThreeVector(0, -shapePar[0] * sin(shapePar[1]), shapePar[0] * cos(shapePar[1]));
			cRotation.rotateY(-shapePar[1]);
			cRotation.rotateZ(CLHEP::pi / 2);
		}

		if (i == 4)
		{
			cRotation.rotateY(shapePar[2]);
			cPosition = G4ThreeVector(shapePar[0] * sin(shapePar[2]), 0, shapePar[0] * cos(shapePar[2]));
		}

		if (i == 5)
		{
			cRotation.rotateY(-shapePar[2]);
			cPosition = G4ThreeVector(-shapePar[0] * sin(shapePar[2]), 0, shapePar[0] * cos(shapePar[2]));
		}

		if (i == 6)
		{
			cPosition = G4ThreeVector(0, shapePar[0] * sin(shapePar[2]), shapePar[0] * cos(shapePar[2]));
			cRotation.rotateY(-shapePar[2]);
			cRotation.rotateZ(-CLHEP::pi / 2);
		}

		if (i == 7)
		{
			cPosition = G4ThreeVector(0, -shapePar[0] * sin(shapePar[2]), shapePar[0] * cos(shapePar[2]));
			cRotation.rotateY(-shapePar[2]);
			cRotation.rotateZ(CLHEP::pi / 2);
		}

		G4AssemblyVolume *assembly = new G4AssemblyVolume();

		if (isSegmented)
		{

			assembly = detector->ConstructSegmCloverAssembly(object.GetObjectName() + "_" + "Clover" + to_string(i));
		}
		else
		{

			assembly = detector->ConstructCloverAssembly(object.GetObjectName() + "_" + "Clover" + to_string(i));
		}

		ELIADEassembly->AddPlacedAssembly(assembly, cPosition, &cRotation);

		auto iterator = assembly->GetVolumesIterator();

		for (unsigned int it = 0; it < assembly->TotalImprintedVolumes(); it++)
		{

			if ((*iterator)->CheckOverlaps())
			{

				cout << " " << endl;
				cout << "------------------- Warning -------------------" << endl;
				cout << "Overlap detected for volume \"" << (*iterator)->GetName() << "\" part of assembly \""
					 << object.GetObjectName() << "\"" << endl;
				cout << " " << endl;
			}

			iterator++;
		}

		vector<G4String> activeVolume = detector->GetActiveVolumes();

		detector->ClearActiveVolumes();

		for (unsigned int l = 0; l < activeVolume.size(); l++)
		{
			FInput->GeometryObjects[objectId].AddObjectActiveVolume(activeVolume[l]);
		}

		FInput->GeometryObjects[objectId].SetObjectFWHMParam(detector->GetFWHMParameters());
	}

	ELIADEassembly->MakeImprint(logicWorld, position, &rotation,
								0, object.GetObjectOverlap());

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());

	AppsAnalysisObject *clovanalysis = new AppsAnalysisObject();
	clovanalysis->SetObjectName(object.GetObjectName() + "_" + "Clover0");

	for (unsigned int i = 0; i < FInput->GetGeometryObjects()[objectId].AnalysisObjects.size(); i++)
	{
		if (isSegmented)
		{

			if (i % 32 == 0)
			{

				clovanalysis = new AppsAnalysisObject();
				clovanalysis->SetObjectName(object.GetObjectName() + "_" + "Clover" + to_string(i / 32));
			}

			FInput->GeometryObjects[objectId].AnalysisObjects[i].push_back(clovanalysis);
			FInput->GeometryObjects[objectId].AnalysisObjects[i].push_back(analysis);
		}
		else
		{

			if (i % 4 == 0)
			{

				clovanalysis = new AppsAnalysisObject();
				clovanalysis->SetObjectName(object.GetObjectName() + "_" + "Clover" + to_string(i / 4));
			}

			FInput->GeometryObjects[objectId].AnalysisObjects[i].push_back(clovanalysis);
			FInput->GeometryObjects[objectId].AnalysisObjects[i].push_back(analysis);
		}
	}

	/*
	for(int i = 0; i < FInput->GetGeometryObjects()[objectId].AnalysisObjects.size(); i++){
		for(int o = 0; o < FInput->GetGeometryObjects()[objectId].AnalysisObjects[i].size(); o++){
			cout << FInput->GetGeometryObjects()[objectId].AnalysisObjects[i][o]->GetObjectName() << endl;
		}
	}
*/
}

void AppsDetectorConstruction::ConstructCloverDetectorObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	CloverDetector *detector = new CloverDetector(objectId);

	G4ThreeVector position = object.GetObjectPosition();

	G4AssemblyVolume *assembly = new G4AssemblyVolume();

	G4bool isSegmented = false;

	for (unsigned int o = 0; o < object.GetObjectFeatures().size(); o++)
	{
		if (object.GetObjectFeatures()[o] == "segmented")
		{
			isSegmented = true;
		}
	}

	if (isSegmented)
	{
		assembly = detector->ConstructSegmCloverAssembly(object.GetObjectName());
	}
	else
	{
		assembly = detector->ConstructCloverAssembly(object.GetObjectName());
	}

	assembly->MakeImprint(logicWorld, position, object.GetObjectRotation(),
						  0, object.GetObjectOverlap());

	auto iterator = assembly->GetVolumesIterator();

	for (unsigned int it = 0; it < assembly->TotalImprintedVolumes(); it++)
	{

		if ((*iterator)->CheckOverlaps())
		{

			cout << " " << endl;
			cout << "------------------- Warning -------------------" << endl;
			cout << "Overlap detected for volume \"" << (*iterator)->GetName() << "\" part of assembly \""
				 << object.GetObjectName() << "\"" << endl;
			cout << " " << endl;
		}

		iterator++;
	}

	FInput->GeometryObjects[objectId].SetObjectActiveVolumes(detector->GetActiveVolumes());
	FInput->GeometryObjects[objectId].SetObjectFWHMParam(detector->GetFWHMParameters());

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());

	for (unsigned int i = 0; i < FInput->GetGeometryObjects()[objectId].AnalysisObjects.size(); i++)
	{

		FInput->GeometryObjects[objectId].AnalysisObjects[i].push_back(analysis);
	}

	/*
	for(int i = 0; i < FInput->GetGeometryObjects()[objectId].AnalysisObjects.size(); i++){
		for(int o = 0; o < FInput->GetGeometryObjects()[objectId].AnalysisObjects[i].size(); o++){
			cout << FInput->GetGeometryObjects()[objectId].AnalysisObjects[i][o]->GetObjectName() << endl;
		}
	}
*/
}

void AppsDetectorConstruction::ConstructTubeObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4Material *material = BuildMaterial(object.GetObjectMaterial());

	G4VSolid *Tubs = new G4Tubs(object.GetObjectName(), object.GetShapeParameters()[0], object.GetShapeParameters()[1], object.GetShapeParameters()[2],
								object.GetShapeParameters()[3], object.GetShapeParameters()[4]);
	G4LogicalVolume *logicTube = new G4LogicalVolume(Tubs, material, object.GetObjectName());

	LogicalVolumes.push_back(logicTube);

	bool isVisible = true;

	double Red = 0;
	double Green = 0;
	double Blue = 0;

	G4VisAttributes *tempVisAttr = new G4VisAttributes(isVisible,
													   G4Colour(Red / 255., Green / 255., Blue / 255.));

	//logicTube->SetVisAttributes(tempVisAttr);

	G4PVPlacement *tubePlace = new G4PVPlacement(object.GetObjectRotation(), object.GetObjectPosition(), logicTube, object.GetObjectName(),
												 GetLogicalVolume(object.GetObjectMotherVolume()), false, 0, object.GetObjectOverlap());

	PlacementVolumes.push_back(tubePlace);

	if (tubePlace->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"" << object.GetObjectName() << "\"" << endl;
		cout << " " << endl;
	}

	FInput->GeometryObjects[objectId].AddObjectActiveVolume(object.GetObjectName());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

G4LogicalVolume *AppsDetectorConstruction::GetLogicalVolume(G4String volumeName)
{

	G4LogicalVolume *logVolume = NULL;

	for (unsigned int i = 0; i < LogicalVolumes.size(); i++)
	{

		if (LogicalVolumes[i]->GetName() == volumeName)
		{

			return LogicalVolumes[i];
		}
	}

	return logVolume;
}

void AppsDetectorConstruction::ConstructMonitorObject(G4int objectId)
{

	AppsInput *FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[objectId];

	G4bool isCircular = false;

	for (unsigned int i = 0; i < object.GetObjectFeatures().size(); i++)
	{
		if (object.GetObjectFeatures()[i] == "circular")
		{
			isCircular = true;
		}
	}

	G4Material *material = BuildMaterial(object.GetObjectMaterial());

	G4VSolid *monitor;

	if (isCircular)
	{

		monitor = new G4Tubs(object.GetObjectName(), object.GetShapeParameters()[0], object.GetShapeParameters()[1],
							 object.GetShapeParameters()[3], 0, object.GetShapeParameters()[2]);
	}
	else
	{

		monitor = new G4Box(object.GetObjectName(), object.GetShapeParameters()[0], object.GetShapeParameters()[1], object.GetShapeParameters()[2]);
	}

	G4LogicalVolume *logicMonitor = new G4LogicalVolume(monitor, material, object.GetObjectName());

	LogicalVolumes.push_back(logicMonitor);

	G4PVPlacement *monitorPlace = new G4PVPlacement(object.GetObjectRotation(), object.GetObjectPosition(), logicMonitor, object.GetObjectName(),
													GetLogicalVolume(object.GetObjectMotherVolume()), false, 0, object.GetObjectOverlap());

	if (monitorPlace->CheckOverlaps())
	{

		cout << " " << endl;
		cout << "------------------- Warning -------------------" << endl;
		cout << "Overlap detected for object \"" << object.GetObjectName() << "\"" << endl;
		cout << " " << endl;
	}

	FInput->GeometryObjects[objectId].AddObjectActiveVolume(FInput->GetGeometryObjects()[objectId].GetObjectName());

	vector<AppsAnalysisObject *> temp;

	AppsAnalysisObject *analysis = new AppsAnalysisObject();
	analysis->SetObjectName(object.GetObjectName());
	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());
	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);
}

void AppsDetectorConstruction::MoveObjects()
{

	AppsInput *FInput = AppsInput::Instance();

	auto run = static_cast<const AppsRun *>(G4RunManager::GetRunManager()->GetCurrentRun());

	G4int currentRunId;

	if (run == NULL)
	{
		currentRunId = 0;
	}
	else
	{
		currentRunId = run->GetRunID() + 1;
	}

	for (unsigned int i = 0; i < FInput->GetGeometryObjects().size(); i++)
	{

		if (FInput->GetGeometryObjects()[i].GetMovingParameters().size() != 0)
		{

			vector<G4double> movPar = FInput->GetGeometryObjects()[i].GetMovingParameters();
			G4ThreeVector currentPos = FInput->GetGeometryObjects()[i].GetObjectPosition();
			auto objName = FInput->GetGeometryObjects()[i].GetObjectName();

			int xStep = (currentRunId % int(movPar[0]));
			int yStep = (currentRunId / int(movPar[0]));

			double xPos = -(movPar[0] * movPar[1] - movPar[1]) / 2.0 + xStep * movPar[1];
			double yPos = -(movPar[2] * movPar[3] - movPar[3]) / 2.0 + yStep * movPar[3];

			G4ThreeVector someVect(xPos, yPos, 0);

			currentPos = currentPos + someVect;

			for (unsigned int o = 0; o < PlacementVolumes.size(); o++)
			{

				if (PlacementVolumes[o]->GetName() == objName)
				{

					PlacementVolumes[o]->SetTranslation(currentPos);
				}
			}
		}
	}
}

void AppsDetectorConstruction::RotateObject(G4String objName, G4double angle)
{

	for (unsigned int i = 0; i < PlacementVolumes.size(); i++)
	{

		if (PlacementVolumes[i]->GetName() == objName)
		{

			G4RotationMatrix *theRot = new G4RotationMatrix;
			theRot->rotateZ(angle * deg);

			G4RotationMatrix *currentRot = PlacementVolumes[i]->GetRotation();
			currentRot->rotateZ(angle * deg);

			PlacementVolumes[i]->SetRotation(currentRot);
		}
	}

	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void AppsDetectorConstruction::ConstructSDandField()
{

	AppsInput *FInput = AppsInput::Instance();

	for (unsigned int i = 0; i < FInput->GetGeometryObjects().size(); i++)
	{
		if (FInput->GetGeometryObjects()[i].GetObjectActiveStatus())
		{

			for (unsigned int o = 0; o < FInput->GetGeometryObjects()[i].GetAnalysisObjects().size(); o++)
			{

				G4String SD = FInput->GetGeometryObjects()[i].GetAnalysisObjects()[o][0]->GetObjectName() + "SD";
				AppsSD *appsSD = new AppsSD(SD, FInput->GetGeometryObjects()[i].GetAnalysisObjects()[o][0]->GetObjectName() + "HitsCollection");
				G4SDManager::GetSDMpointer()->AddNewDetector(appsSD);
				SetSensitiveDetector(FInput->GetGeometryObjects()[i].GetAnalysisObjects()[o][0]->GetObjectName(), appsSD, true);
			}
		}
	}
}

G4VPhysicalVolume *AppsDetectorConstruction::Construct()
{

	AppsInput *FInput = AppsInput::Instance();

	ConstructWorld();
	
	if (FInput->GetBeamType() == "PaNi")
	{
		
		ConstructPaNiSphere();
		
	}

	for (unsigned int i = 0; i < FInput->GetGeometryObjects().size(); i++)
	{

		// part of the code in charge of constructing the box type object

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "Box")
		{

			ConstructBoxObject(i);
		}

		// part of the code in charge of constructing the tube type object

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "Tube")
		{

			ConstructTubeObject(i);
		}

		// part of the code in charge of constructing the sphere type object

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "Sphere")
		{

			ConstructSphereObject(i);
		}

		// part of the code in charge of constructing the 150 % relative efficiency germanium part of the diagnostics setup

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "DiagDetector")
		{

			ConstructDiagDetectorObject(i);
		}

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "ScintDetector")
		{

			ConstructScintDetectorObject(i);
		}

		// part of the code in charge of constructing the clover germanium part of eliade,
		// can be run in segmented ( - features = [segmented]  or unsegmented mode (default) )

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "CloverDetector")
		{

			ConstructCloverDetectorObject(i);
		}

		// part of code in charge of construting a monitor, thin box like object made out of vacuum in order of probe the photons that pass in that specific possition

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "Monitor")
		{

			ConstructMonitorObject(i);
		}

		// part of code in charge of constructing eliade

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "ELIADE")
		{
			ConstructELIADEObject(i);
		}

		// part of code in charge of constructing rectangular collimator with round hole

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "Collimator")
		{

			ConstructCollimatorObject(i);
		}

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "PIXEL")
		{

			ConstructPIXELObject(i);
		}

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "Mask")
		{

			ConstructMaskObject(i);
		}

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "CeBrDetector")
		{

			ConstructCeBrDetectorObject(i);
		}

		if (FInput->GetGeometryObjects()[i].GetObjectType() == "VEGA")
		{
			ConstructVEGACollimator(i);
		}
	}
	return physWorld;
}

G4Material *AppsDetectorConstruction::BuildMaterial(vector<G4String> objectMaterial)
{

	G4NistManager *nist = G4NistManager::Instance();

	G4Material *material;

	if (objectMaterial.size() == 1)
	{ // multiple materials implementation is not written

		material = nist->FindOrBuildMaterial("G4_" + objectMaterial[0]);
		if (material)
		{
			return material;
		}
		else
		{

			if (objectMaterial[0] == "NaI")
			{

				material = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
				return material;
			}

			if (objectMaterial[0] == "HighVacuum")
			{ // defined high vacuum

				G4double density = universe_mean_density;
				G4double pressure = 1.e-9 * pascal;
				G4double temp = 0.1 * kelvin;

				material = new G4Material("HighVacuum", 1., 1.01 * g / mole, density, kStateGas, temp, pressure);

				return material;
			}

			if (objectMaterial[0] == "DepletedUranium")
			{ // define depleted uranium

				G4Isotope *iU235 = new G4Isotope("U235", 92, 235, 235.01 * g / mole);
				G4Isotope *iU238 = new G4Isotope("U238", 92, 238, 238.03 * g / mole);

				G4Element *elU238 = new G4Element("Uranium 238", "U238", 1);
				elU238->AddIsotope(iU238, 100. * perCent);

				G4Element *elU235 = new G4Element("Uranium 235", "U235", 1);
				elU235->AddIsotope(iU235, 100. * perCent);

				G4double density = 19.05 * g / cm3;
				material = new G4Material("depletedUranium", density, 2);
				material->AddElement(elU238, 1.0);
				material->AddElement(elU235, 0.0);

				return material;
			}

			if (objectMaterial[0] == "SomeUranium")
			{ // define some uranium

				G4Isotope *iU235 = new G4Isotope("U235", 92, 235, 235.01 * g / mole);
				G4Isotope *iU238 = new G4Isotope("U238", 92, 238, 238.03 * g / mole);

				G4Element *uranium = new G4Element("Some Uranium", "U238", 2);
				uranium->AddIsotope(iU238, 10. * perCent);
				uranium->AddIsotope(iU235, 90. * perCent);

				G4double density = 19.05 * g / cm3;
				material = new G4Material("someUranium", density, 1);
				material->AddElement(uranium, 1.0);

				return material;
			}

			if (objectMaterial[0] == "EnrichedUranium")
			{ // define enriched uranium

				G4Isotope *iU235 = new G4Isotope("U235", 92, 235, 235.01 * g / mole);
				G4Isotope *iU238 = new G4Isotope("U238", 92, 238, 238.03 * g / mole);

				G4Element *elU238 = new G4Element("Uranium 238", "U238", 1);
				elU238->AddIsotope(iU238, 100. * perCent);

				G4Element *elU235 = new G4Element("Uranium 235", "U235", 1);
				elU235->AddIsotope(iU235, 100. * perCent);

				G4double density = 19.05 * g / cm3;
				material = new G4Material("depletedUranium", density, 2);
				material->AddElement(elU238, 0.01);
				material->AddElement(elU235, 0.99);

				return material;
			}

			if (objectMaterial[0] == "Thorium")
			{ // define thorium

				G4Isotope *iTh232 = new G4Isotope("Th232", 90, 232, 232.038 * g / mole);

				G4Element *elTh232 = new G4Element("Thorium", "Th232", 1);
				elTh232->AddIsotope(iTh232, 100. * perCent);

				G4double density = 11.72 * g / cm3;
				material = new G4Material("Thorium", density, 1);
				material->AddElement(elTh232, 1.0);

				return material;
			}

			if (objectMaterial[0] == "Boron")
			{ // define thorium

				G4Isotope *iB11 = new G4Isotope("B11", 5, 11, 11.009 * g / mole);
				G4Isotope *iB10 = new G4Isotope("B10", 5, 10, 10.013 * g / mole);

				G4Element *elB = new G4Element("Boron", "B11", 2);
				elB->AddIsotope(iB11, 80. * perCent);
				elB->AddIsotope(iB10, 20. * perCent);

				G4double density = 2.08 * g / cm3;
				material = new G4Material("Boron", density, 1);
				material->AddElement(elB, 1.0);

				return material;
			}

			if (objectMaterial[0] == "Sulfur")
			{ // define S32

				G4Isotope *iS32 = new G4Isotope("S32", 16, 32, 32.065 * g / mole);

				G4Element *elS32 = new G4Element("Sulfur", "S32", 1);
				elS32->AddIsotope(iS32, 100. * perCent);

				G4double density = 2.07 * g / cm3;
				material = new G4Material("Sulfur", density, 1);
				material->AddElement(elS32, 1.0);

				return material;
			}

			if (objectMaterial[0] == "Iron")
			{ // define Fe56

				G4Isotope *iFe56 = new G4Isotope("Fe56", 26, 56, 55.845 * g / mole);

				G4Element *elFe56 = new G4Element("Iron", "Fe56", 1);
				elFe56->AddIsotope(iFe56, 100. * perCent);

				G4double density = 7.874 * g / cm3;
				material = new G4Material("Iron", density, 1);
				material->AddElement(elFe56, 1.0);

				return material;
			}

			if (objectMaterial[0] == "EnrichedZirconia96")
			{ // define zirconium
				double fractionmass;

				G4Isotope *iZr96 = new G4Isotope("Zr96", 40, 96, 95.908 * g / mole);
				G4Isotope *iZr90 = new G4Isotope("Zr90", 40, 90, 89.905 * g / mole);
				G4Isotope *iO = new G4Isotope("O", 8, 16, 15.999 * g / mole);

				G4Element *elZr96 = new G4Element("Zirconium 96", "Zr96", 1);
				elZr96->AddIsotope(iZr96, 100. * perCent);

				G4Element *elZr90 = new G4Element("Zirconium 90", "Zr90", 1);
				elZr90->AddIsotope(iZr90, 100. * perCent);

				G4Element *elO = new G4Element("oxygen", "O", 1);
				elO->AddIsotope(iO, 100. * perCent);

				G4double density = 5.68 * g / cm3;
				material = new G4Material("EnrichedZirconia96", density, 3);
				material->AddElement(elZr96, fractionmass = 0.794342);
				material->AddElement(elZr90, fractionmass = 0.125035);
				material->AddElement(elO, fractionmass = 0.080623);

				return material;
			}

			if (objectMaterial[0] == "EnrichedPb208")
			{ // define enriched Pb208

				G4Isotope *iPb208 = new G4Isotope("Pb208", 82, 208, 207.9766521 * g / mole);

				G4Element *elPb208 = new G4Element("Lead 208", "Pb208", 1);
				elPb208->AddIsotope(iPb208, 100. * perCent);

				G4double density = 11.34 * g / cm3;
				material = new G4Material("enrichedPb", density, 1);
				material->AddElement(elPb208, 1.0);

				return material;
			}

			if (objectMaterial[0] == "EnrichedPb206")
			{ // define enriched Pb206

				G4Isotope *iPb206 = new G4Isotope("Pb206", 82, 206, 205.974465 * g / mole);

				G4Element *elPb206 = new G4Element("Lead 206", "Pb206", 1);
				elPb206->AddIsotope(iPb206, 100. * perCent);

				G4double density = 11.34 * g / cm3;
				material = new G4Material("enrichedPb206", density, 1);
				material->AddElement(elPb206, 1.0);

				return material;
			}

			if (objectMaterial[0] == "StainlessSteel")
			{ // define enriched uranium

				G4double density;
				G4int ncomponents;
				G4double fractionmass;

				G4Element *C = nist->FindOrBuildElement("C");
				G4Element *Si = nist->FindOrBuildElement("Si");
				G4Element *Cr = nist->FindOrBuildElement("Cr");
				G4Element *Mn = nist->FindOrBuildElement("Mn");
				G4Element *Fe = nist->FindOrBuildElement("Fe");
				G4Element *Ni = nist->FindOrBuildElement("Ni");

				material = new G4Material("StainlessSteel", density = 8.06 * g / cm3, ncomponents = 6);
				material->AddElement(C, fractionmass = 0.001);
				material->AddElement(Si, fractionmass = 0.007);
				material->AddElement(Cr, fractionmass = 0.18);
				material->AddElement(Mn, fractionmass = 0.01);
				material->AddElement(Fe, fractionmass = 0.712);
				material->AddElement(Ni, fractionmass = 0.09);

				return material;
			}

			if (objectMaterial[0] == "Densimet")
			{
				G4int mcomponents;

				G4Element *W = nist->FindOrBuildElement("W");
				G4Element *Fe = nist->FindOrBuildElement("Fe");
				G4Element *Ni = nist->FindOrBuildElement("Ni");

				material = new G4Material("Densimet", 18.55 * g / cm3, mcomponents = 3);
				material->AddElement(W, 0.97);
				material->AddElement(Fe, 0.015);
				material->AddElement(Ni, 0.015);
			}

			if (objectMaterial[0] == "Explosive")
			{
				G4int mcomponents;

				G4Element *C = nist->FindOrBuildElement("C");
				G4Element *N = nist->FindOrBuildElement("N");
				G4Element *O = nist->FindOrBuildElement("O");

				material = new G4Material("Explosive", 1.46 * g / cm3, mcomponents = 3);
				material->AddElement(C, 0.33);
				material->AddElement(N, 0.33);
				material->AddElement(O, 0.33);
			}

			if (objectMaterial[0] == "Melamine")
			{
				G4int mcomponents;

				G4Element *C = nist->FindOrBuildElement("C");
				G4Element *N = nist->FindOrBuildElement("N");
				G4Element *H = nist->FindOrBuildElement("H");

				material = new G4Material("Melamine", 1.57 * g / cm3, mcomponents = 3);
				material->AddElement(C, 0.29);
				material->AddElement(N, 0.67);
				material->AddElement(H, 0.04);
			}
			if (objectMaterial[0] == "Water")
			{
				G4int mcomponents;
				G4int natoms;

				G4Element *O = nist->FindOrBuildElement("O");

				G4Element *H = nist->FindOrBuildElement("H");

				material = new G4Material("Water", 0.997 * g / cm3, mcomponents = 2);
				material->AddElement(O, natoms = 1);
				material->AddElement(H, natoms = 2);
			}

			if (objectMaterial[0] == "CZT")
			{
				G4int mcomponents;
				G4int natoms;

				G4Element *Cd = nist->FindOrBuildElement("Cd");

				G4Element *Zn = nist->FindOrBuildElement("Zn");
				G4Element *Te = nist->FindOrBuildElement("Te");

				material = new G4Material("CZT", 5.78 * g / cm3, mcomponents = 3);
				material->AddElement(Cd, natoms = 1);
				material->AddElement(Zn, natoms = 1);
				material->AddElement(Te, natoms = 1);
			}

			if (objectMaterial[0] == "LYSO")
			{
				G4int ncomponents;
				G4int natoms;
				G4Element *Lu = nist->FindOrBuildElement("Lu");
				G4Element *Y = nist->FindOrBuildElement("Y");
				G4Element *Si = nist->FindOrBuildElement("Si");
				G4Element *O = nist->FindOrBuildElement("O");
				G4Element *Ce = nist->FindOrBuildElement("Ce");

				G4Material *LYSO = new G4Material("LYSO", 7.1 * g / cm3, ncomponents = 4);
				LYSO->AddElement(Lu, natoms = 1);
				LYSO->AddElement(Y, natoms = 1);
				LYSO->AddElement(Si, natoms = 1);
				LYSO->AddElement(O, natoms = 5);
				material = new G4Material("LYSO_Ce", 7.1 * g / cm3,
										  ncomponents = 2);
				material->AddMaterial(LYSO, 99.5 * perCent);
				material->AddElement(Ce, 0.5 * perCent);
			}

			if (objectMaterial[0] == "PaNi1")
			{
				G4int ncomponents;
				G4Material *Pa = nist->FindOrBuildMaterial("G4_PARAFFIN");
				G4Element *Ni = nist->FindOrBuildElement("Ni");
				
				G4double density = 1.5 * g / cm3;
				material = new G4Material("PaNi1", density, ncomponents = 2);
				material->AddMaterial(Pa, (density/g*cm3-0.93)/(8.9-0.93));
				material->AddElement(Ni, (8.9-density/g*cm3)/(8.9-0.93));
			}

			if (objectMaterial[0] == "PaNi2")
			{
				G4int ncomponents;
				G4Material *Pa = nist->FindOrBuildMaterial("G4_PARAFFIN");
				G4Element *Ni = nist->FindOrBuildElement("Ni");
				
				G4double density = 1.087 * g / cm3;
				material = new G4Material("PaNi2", density, ncomponents = 2);
				material->AddMaterial(Pa, (density/g*cm3-0.93)/(8.9-0.93));
				material->AddElement(Ni, (8.9-density/g*cm3)/(8.9-0.93));
			}

			if (objectMaterial[0] == "PaNi3")
			{
				G4int ncomponents;
				G4Material *Pa = nist->FindOrBuildMaterial("G4_PARAFFIN");
				G4Element *Ni = nist->FindOrBuildElement("Ni");
				
				G4double density = 2. * g / cm3;
				material = new G4Material("PaNi3", density, ncomponents = 2);
				material->AddMaterial(Pa, (density/g*cm3-0.93)/(8.9-0.93));
				material->AddElement(Ni, (8.9-density/g*cm3)/(8.9-0.93));
			}
		}
	}
	else
	{

		if (objectMaterial[0] == "Uranium238")
		{ // define some uranium

			double Proc238 = stod(objectMaterial[1]);

			G4Isotope *iU235 = new G4Isotope("U235", 92, 235, 235.01 * g / mole);
			G4Isotope *iU238 = new G4Isotope("U238", 92, 238, 238.03 * g / mole);

			G4Element *uranium = new G4Element("Some Uranium", "U238", 2);
			uranium->AddIsotope(iU238, Proc238 * perCent);
			uranium->AddIsotope(iU235, (100 - Proc238) * perCent);

			G4double density = 19.05 * g / cm3;
			material = new G4Material("someUranium", density, 1);
			material->AddElement(uranium, 1.0);

			return material;
		}

		/*
		material = new G4Material("material"+to_string(NewMaterialCount), stod(objectMaterial.back()),objectMaterial.size()-1);

		G4double fSum =0;
		vector<G4double> frac;
		for(unsigned int i = 0 ; i < objectMaterial.size()-1; i++){
			std::size_t firstpos = objectMaterial[i].find("-")+1;
			frac.push_back(stod(objectMaterial[i].substr(firstpos)));
			fSum+= stod(objectMaterial[i].substr(firstpos));
		}

		for(unsigned int i = 0 ; i < objectMaterial.size()-1; i++){
		
			std::size_t firstpos = objectMaterial[i].find("-");
			G4Element * elem = nist->FindOrBuildElement(objectMaterial[i].substr(0,firstpos));
			material->AddElement(elem, frac[i]/fSum);

		}

		NewMaterialCount++;
		*/
	}

	if (!material)
	{
		material = nist->FindOrBuildMaterial("G4_AIR");
	}

	return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
