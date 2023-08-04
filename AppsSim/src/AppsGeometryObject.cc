#include "AppsGeometryObject.hh"

AppsGeometryObject::AppsGeometryObject(G4String type)
{

	SetObjectType(type);
	Initialize();
}

AppsGeometryObject::AppsGeometryObject()
{
	Initialize();
}

AppsGeometryObject::~AppsGeometryObject()
{
}

void AppsGeometryObject::SetObjectRotation(G4double rotX, G4double rotY, G4double rotZ)
{

	ObjectRotation = new G4RotationMatrix;
	ObjectRotation->rotateX(rotX);
	ObjectRotation->rotateY(rotY);
	ObjectRotation->rotateZ(rotZ);

	ObjectRotationAngle.clear();
	ObjectRotationAngle.push_back(rotX);
	ObjectRotationAngle.push_back(rotY);
	ObjectRotationAngle.push_back(rotZ);
}

void AppsGeometryObject::Initialize()
{

	SetObjectName("undefined");
	SetObjectMotherVolume("World");
	SetObjectPosition(0 * cm, 0 * cm, 0 * cm);
	SetObjectRotation(0 * deg, 0 * deg, 0 * deg);
	SetObjectOverlap(1);
	AddObjectMaterial("H"); // default intialization material is H
	SetObjectActiveStatus(true);
	AddObjectHistogramType("Edep");

	AddObjectFeature("undefined");

	FWHMParameters.push_back(0); // slope
	FWHMParameters.push_back(0); // intercept

	if (GetObjectType() == "ELIADE")
	{

		AddShapeParameter(15 * cm);
		AddShapeParameter(90 * deg);
		AddShapeParameter(135 * deg);
	}

	if (GetObjectType() == "Collimator")
	{

		AddShapeParameter(1 * cm);
		AddShapeParameter(1 * cm);
		AddShapeParameter(1 * cm);
		AddShapeParameter(0.3 * cm);
	}

	if (GetObjectType() == "Monitor")
	{
		ClearObjectHistogramType();
		AddObjectHistogramType("Eflux");
		AddObjectHistogramType("Sflux");
	}

	if (GetObjectType() == "PIXEL")
	{

		AddShapeParameter(1 * cm);
		AddShapeParameter(1 * cm);
		AddShapeParameter(2);
		AddShapeParameter(2);
	}

	if (GetObjectType() == "VEGA")
	{
		SetObjectActiveStatus(false);
	}
}
