#include "ScintDetector.hh"
#include "G4NistManager.hh"
#include <fstream>
#include "G4SystemOfUnits.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "AppsInput.hh"

using namespace std;

ScintDetector::ScintDetector(G4int id)
{

	AppsInput *FInput = AppsInput::Instance();

	object = FInput->GetGeometryObjects()[id];

	objectId = id;

	ReadInputParameters();
	ConstructMaterials();

	// set linear FWHM parameters

	SetFWHMParameters(7.23, -0.52);
}

ScintDetector::~ScintDetector()
{
}

void ScintDetector::ConstructCrystal(G4String name)
{

	G4VSolid *NaICrystal = new G4Tubs("NaICrystal", 0 * cm, CrystalRadius, CrystalLength / 2.0, 0. * deg, 360. * deg);

	LogicalVolume1 = new G4LogicalVolume(NaICrystal, NaI_mat, name);

	G4VisAttributes *vis = new G4VisAttributes(false);
	LogicalVolume1->SetVisAttributes(vis);
}

void ScintDetector::ConstructQuartzWindow()
{

	G4VSolid *Quartz1 = new G4Tubs("Quartz1", 0 * cm, QuartzWindowRadius, QuartzWindowLength / 2.0, 0. * deg, 360. * deg);

	LogicalVolume3 = new G4LogicalVolume(Quartz1, Quartz_mat, "QuartzWindow");

	G4VisAttributes *vis = new G4VisAttributes(false);
	LogicalVolume3->SetVisAttributes(vis);
}

void ScintDetector::ConstructHolder()
{

	G4VSolid *Holder1 = new G4Tubs("Holder1", HoldeInnerRadius, HoldeOuterRadius, HolderLength / 2.0, 0. * deg, 360. * deg);
	G4VSolid *Holder2 = new G4Tubs("Holder2", 0 * cm, HoldeOuterRadius, HolderFrontThickness / 2.0, 0. * deg, 360. * deg);
	G4VSolid *Holder3 = new G4UnionSolid("Holder3", Holder1, Holder2, 0, G4ThreeVector(0., 0., -(HolderLength / 2.0 + HolderFrontThickness / 2.0)));
	G4VSolid *Holder4 = new G4Tubs("Holder4", FirstRingInnerRadius, FirstRingOuterRadius, FirstRingThickness / 2.0, 0. * deg, 360. * deg);
	G4VSolid *Holder5 = new G4UnionSolid("Holder5", Holder3, Holder4, 0, G4ThreeVector(0., 0., HolderLength / 2.0 + FirstRingThickness / 2.0));
	G4VSolid *Holder6 = new G4Tubs("Holder6", SecondRingInnerRadius, SecondRingOuterRadius, SecondRingThickness / 2.0, 0. * deg, 360. * deg);
	G4VSolid *Holder7 = new G4UnionSolid("Holder7", Holder5, Holder6, 0, G4ThreeVector(0., 0., HolderLength / 2.0 + FirstRingThickness + SecondRingThickness / 2.0));
	G4VSolid *Holder8 = new G4Tubs("Holder8", FirstPhotoRingInnerRadius, FirstPhotoRingOuterRadius, FirstPhotoRingLength / 2.0, 0. * deg, 360. * deg);
	G4VSolid *Holder9 = new G4UnionSolid("Holder9", Holder7, Holder8, 0, G4ThreeVector(0., 0., HolderLength / 2.0 + FirstRingThickness + FirstPhotoRingLength / 2.0));
	G4VSolid *Holder10 = new G4Cons("Holder10", FirstPhotoRingInnerRadius, FirstPhotoRingOuterRadius, SecondPhotoRingInnerRadius,
									SecondPhotoRingOuterRadius, ConeLength / 2.0, 0. * deg, 360. * deg);

	G4VSolid *Holder11 = new G4UnionSolid("Holder11", Holder9, Holder10, 0, G4ThreeVector(0., 0., HolderLength / 2.0 + FirstRingThickness + FirstPhotoRingLength + ConeLength / 2.0));
	G4VSolid *Holder12 = new G4Tubs("Holder12", SecondPhotoRingInnerRadius, SecondPhotoRingOuterRadius, SecondPhotoRingLength / 2.0, 0. * deg, 360. * deg);
	G4VSolid *Holder13 = new G4UnionSolid("Holder13", Holder11, Holder12, 0, G4ThreeVector(0., 0., HolderLength / 2.0 + FirstRingThickness + FirstPhotoRingLength + ConeLength + SecondPhotoRingLength / 2.0));
	G4VSolid *Holder14 = new G4Tubs("Holder14", 0 * cm, SecondPhotoRingOuterRadius, SecondPhotoRingBackThickness / 2.0, 0. * deg, 360. * deg);
	G4VSolid *Holder15 = new G4UnionSolid("Holder15", Holder13, Holder14, 0, G4ThreeVector(0., 0., HolderLength / 2.0 + FirstRingThickness + FirstPhotoRingLength + ConeLength + SecondPhotoRingLength + SecondPhotoRingBackThickness / 2.0));
	LogicalVolume2 = new G4LogicalVolume(Holder15, Al_mat, "Holder");
}

void ScintDetector::ReadInputParameters()
{

	G4String z;
	G4double a[23];

	ifstream myfile("ScintDetParameters.txt");

	if (myfile.is_open())
	{
		getline(myfile, z);
		getline(myfile, z);
		getline(myfile, z);

		for (G4int i = 0; i < 23; i++)
		{

			myfile >> a[i];
			getline(myfile, z);
			getline(myfile, z);
		}
	}
	else
	{
		cout << "Unable to open DiagDetParameters.txt" << G4endl;
	}

	CrystalLength = a[0] * cm;
	CrystalRadius = a[1] * cm;
	HolderLength = a[2] * cm;
	HoldeInnerRadius = a[3] * cm;
	HoldeOuterRadius = a[4] * cm;
	HolderFrontThickness = a[5] * cm;
	FirstRingOuterRadius = a[6] * cm;
	FirstRingInnerRadius = a[7] * cm;
	FirstRingThickness = a[8] * cm;
	SecondRingOuterRadius = a[9] * cm;
	SecondRingInnerRadius = a[10] * cm;
	SecondRingThickness = a[11] * cm;
	QuartzWindowRadius = a[12] * cm;
	QuartzWindowLength = a[13] * cm;
	FirstPhotoRingOuterRadius = a[14] * cm;
	FirstPhotoRingInnerRadius = a[15] * cm;
	FirstPhotoRingLength = a[16] * cm;
	SecondPhotoRingInnerRadius = a[17] * cm;
	SecondPhotoRingOuterRadius = a[18] * cm;
	SecondPhotoRingLength = a[19] * cm;
	SecondPhotoRingBackThickness = a[20] * cm;
	CrystalToHolderGap = a[21] * cm;
	ConeLength = a[22] * cm;
}

void ScintDetector::ConstructMaterials()
{

	G4NistManager *nist = G4NistManager::Instance();

	G4Element *O = new G4Element("Oxygen", "O", 8, 16.00 * g / mole);
	G4Element *Si = new G4Element("Silicon", "Si", 14., 28.09 * g / mole);

	Quartz_mat = new G4Material("Quartz", 2.32 * g / cm3, 2);
	Quartz_mat->AddElement(Si, 1);
	Quartz_mat->AddElement(O, 2);

	Ge_mat = nist->FindOrBuildMaterial("G4_Ge");
	Al_mat = nist->FindOrBuildMaterial("G4_Al");
	NaI_mat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
}

G4AssemblyVolume *ScintDetector::GetDetectorAssembly()
{

	ConstructCrystal(object.GetObjectName());
	ConstructHolder();
	ConstructQuartzWindow();

	G4AssemblyVolume *assemblyDetector1 = new G4AssemblyVolume();

	G4ThreeVector pos0 = G4ThreeVector(0, 0, CrystalLength / 2.0 + CrystalToHolderGap + HolderFrontThickness);
	G4ThreeVector pos1 = G4ThreeVector(0, 0, 0) + pos0;
	G4ThreeVector pos2 = G4ThreeVector(0, 0, HolderLength + HolderFrontThickness - CrystalLength - CrystalToHolderGap) + pos0;
	G4ThreeVector pos3 = G4ThreeVector(0, 0, CrystalLength / 2.0 + QuartzWindowLength / 2.0) + pos0;

	G4RotationMatrix *yRot2 = new G4RotationMatrix;
	yRot2->rotateY(0 * deg);

	assemblyDetector1->AddPlacedVolume(LogicalVolume1, pos1, yRot2);

	bool isVisible = true;

	double Red = 0;
	double Green = 126;
	double Blue = 128;

	G4VisAttributes *tempVisAttr = new G4VisAttributes(isVisible,
													   G4Colour(Red / 255., Green / 255., Blue / 255.));

	LogicalVolume2->SetVisAttributes(tempVisAttr);

	assemblyDetector1->AddPlacedVolume(LogicalVolume2, pos2, yRot2);
	assemblyDetector1->AddPlacedVolume(LogicalVolume3, pos3, yRot2);

	AddActiveVolumeName(object.GetObjectName());

	return assemblyDetector1;
}

void ScintDetector::SetFWHMParameters(G4double slope, G4double intercept)
{

	FWHMParameters.push_back(slope);
	FWHMParameters.push_back(intercept);
}
