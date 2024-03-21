#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "AppsGeometryObject.hh"
#include "AppsAnalysisObject.hh"
#include "AppsInput.hh"

class CloverDetector
{

public:
	CloverDetector(G4int id);
	~CloverDetector();

	void ConstructMaterials();
	void ReadInputParameters();
	void ConstructCrystal();

	void ConstructBackDeadLayer();
	void ConstructDeadLayer();
	void ConstructDetectorCrystal();
	void ConstructInsideDeadLayer();
	void ConstructCrystalSegments();
	void ConstructHolder(G4String name);
	void ConstructBackHolder();
	void ConstructShield();

	G4AssemblyVolume *ConstructCrystalAssembly(G4String name);
	G4AssemblyVolume *ConstructSegmCloverAssembly(G4String name);
	G4AssemblyVolume *ConstructCloverAssembly(G4String name);
	G4AssemblyVolume *ConstructSegCrystalAssembly(G4String name);
	G4AssemblyVolume *ConstructShieldAssembly(G4String name);

	G4LogicalVolume *logicShape4;
	G4LogicalVolume *logicShape1;
	G4LogicalVolume *logicShape3;
	G4LogicalVolume *logicShape6;
	G4LogicalVolume *logicShape2;
	G4LogicalVolume *logicShape5;
	G4LogicalVolume *logicShape7 = NULL;

	G4LogicalVolume *Segment1;
	G4LogicalVolume *Segment2;
	G4LogicalVolume *Segment3;
	G4LogicalVolume *Segment4;
	G4LogicalVolume *Segment5;
	G4LogicalVolume *Segment6;
	G4LogicalVolume *Segment7;
	G4LogicalVolume *Segment8;

	G4LogicalVolume *FrontShield1;
	G4LogicalVolume *FrontShield2;
	G4LogicalVolume *FrontShield3;
	G4LogicalVolume *FrontShield4;
	G4LogicalVolume *SideShield1;
	G4LogicalVolume *SideShield2;
	G4LogicalVolume *SideShield3;
	G4LogicalVolume *SideShield4;
	G4LogicalVolume *BackShield1;
	G4LogicalVolume *BackShield2;
	G4LogicalVolume *FrontShieldb1;
	G4LogicalVolume *FrontShieldb2;
	G4LogicalVolume *FrontShieldb3;
	G4LogicalVolume *FrontShieldb4;
	G4LogicalVolume *SideShieldb1;
	G4LogicalVolume *SideShieldb2;
	G4LogicalVolume *SideShieldb3;
	G4LogicalVolume *SideShieldb4;
	G4LogicalVolume *BackShieldb1;
	G4LogicalVolume *BackShieldb2;
	G4LogicalVolume *FrontShieldc1;
	G4LogicalVolume *FrontShieldc2;
	G4LogicalVolume *FrontShieldc3;
	G4LogicalVolume *FrontShieldc4;

	// the active volumes of the detector

	void AddActiveVolumeName(G4String active) { ActiveVolumeNames.push_back(active); }
	void ClearActiveVolumes() { ActiveVolumeNames.clear(); }
	vector<G4String> GetActiveVolumes() { return ActiveVolumeNames; }

	// the FWHM parameters for the resolution

	void SetFWHMParameters(G4double slope, G4double intercept);
	vector<G4double> GetFWHMParameters() { return FWHMParameters; }

private:
	G4int objectId;

	vector<G4String> ActiveVolumeNames;
	vector<G4double> FWHMParameters;

	G4Material *Air_mat;
	G4Material *Ge_mat;
	G4Material *Al_mat;
	G4Material *BGO_mat;
	G4Material *CsI_mat;
	G4Material *Pb_mat;
	G4Material *StainlessS_mat;
	G4Material *Densimet_mat;

	G4double CrystalRadius;
	G4double CrystalLength;
	G4double CrystalInsideLength;
	G4double CrystalOutsideLength;
	G4double FrontDeadLayer;
	G4double BackDeadLayer;
	G4double LateralDeadLayer;
	G4double WellRoundRadius;
	G4double WellRadius;
	G4double WellBottomDeadLayer;
	G4double WellLateralDeadLayer;
	G4double WellLength;
	G4double HolderLength;
	G4double HolderWidth;
	G4double HolderRoundRadius;
	G4double HolderFrontRoundRadius;
	G4double HolderFrontLength;
	G4double HolderFrontWidth;
	G4double HolderWallThickness;
	G4double HolderFrontWallThickness;
	G4double SpaceBetweenCrystals;
	G4double HolderBackWidth;
	G4double HolderBackThickness;
	G4double HolderBackPieceALength;
	G4double HolderBackPieceAWidth;
	G4double HolderBackPieceAThickness;
	G4double HolderBackPieceBWidth;
	G4double HolderBackPieceCRadius;
	G4double HolderBackPieceCThickness;
	G4double HolderBackHoleRadius;
	G4double HolderBackHolePosition;
	G4double HolderBackPieceLength2;
	G4double InsideAluminiumCylinderLength;
	G4double CrystalCutAngle;
	G4double HolderAngle;
	G4double HolderToCrystalDistance;
	G4double CrystalFrontLength;
	G4double FrontToBackHolderDistance;
	G4double CrystalCutDistance;
	G4double CrystalCutDeadLayer;
	G4double sideShieldaL = 12.5*cm;
	G4double sideShieldaH = 5.6*cm;
	G4double sideShieldbL = 17.5*cm;
	G4double sideShieldbH = 14.2*cm;
	G4double sideShieldaaL = (12.5+2*0.05)*cm;
	G4double sideShieldaaH = (5.6-0.4/tan(22.5*deg))*cm;
	G4double sideShieldabL = (12.5+2.*0.05+2.*2.)*cm;
	G4double sideShieldabH = 6.7*cm;
	G4double sideShieldZ;
	G4double backShieldalT = 0.075*cm;
	G4double backShieldaH = 3.75*cm;
	G4double backShieldL = 12.35*cm;
	G4double backShieldH = 12.15*cm;
	G4double backShieldRH = 5.3*cm;
	G4double backShieldR1 = 2.55*cm;
	G4double backShieldR2 = 3.1*cm;
	G4double frontShieldaH = (13.0-1.3)*cm;
	G4double frontShieldaT = 1.0*cm;
	G4double frontShieldcH = 1.3*cm;
	G4double frontShieldcT = 1.3*cm;
	G4double frontShieldcl = 10.63*cm;
	G4double frontShieldaL;
	G4double frontShieldcL;
	G4double frontShieldal;
	G4double frontShieldbH = 12.6*cm;
	G4double frontShieldbT = 2.35*cm;
	G4double frontShieldbL = 22.3*cm;
	G4double frontShieldH = 24.1*cm;
	G4double frontShieldZ;
	G4double length;
	G4double err;

	G4VSolid *GeCrystal6;
	G4VSolid *GeCrystal7;
	G4VSolid *InsideCylinder1;
	G4VSolid *DeadGeCrystal8;
	G4VSolid *BackDeadLayer5;
	G4VSolid *Seg1;
	G4VSolid *Seg2;
	G4VSolid *Seg3;
	G4VSolid *Seg4;
	G4VSolid *Seg5;
	G4VSolid *Seg6;
	G4VSolid *Seg7;
	G4VSolid *Seg8;
	G4VSolid *FrontShield;
	G4VSolid *FrontShieldb;
	G4VSolid *FrontShieldc;
	G4VSolid *SideShield;
	G4VSolid *SideShieldb;
	G4VSolid *BackShield;
	G4VSolid *BackShieldb;
};
