#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "AppsGeometryObject.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class AppsDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	AppsDetectorConstruction();
	virtual ~AppsDetectorConstruction();

	virtual G4VPhysicalVolume *Construct();

	void ConstructSDandField();

	void ConstructWorld(); // function in charge of constructing the world
	void ConstructBoxObject(G4int objectId);
	void ConstructTubeObject(G4int objectId);
	void ConstructMonitorObject(G4int objectId);
	void ConstructWorldObject(AppsGeometryObject object);
	void ConstructDiagDetectorObject(G4int objectId);
	void ConstructScintDetectorObject(G4int objectId);
	void ConstructCloverDetectorObject(G4int objectId);
	void ConstructELIADEObject(G4int objectId);
	void ConstructSphereObject(G4int objectId);
	void ConstructPIXELObject(G4int objectId);
	void ConstructMaskObject(G4int objectId);
	void ConstructCeBrDetectorObject(G4int objectId);
	void ConstructCollimatorObject(G4int objectId);
	void ConstructVEGACollimator(G4int objectId);
	void ConstructPaNiSphere();

	void RotateObject(G4String objName, G4double angle);
	void MoveObjects();

	G4LogicalVolume *GetLogicalVolume(G4String volumeName);

	G4Material *world_mat;

	G4Material *BuildMaterial(vector<G4String> objectMaterial);

	G4VPhysicalVolume *physWorld;
	G4LogicalVolume *logicWorld;

	vector<G4LogicalVolume *> LogicalVolumes;
	vector<G4PVPlacement *> PlacementVolumes;

	G4int NewMaterialCount = 0;

protected:
};
