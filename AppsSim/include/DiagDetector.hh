#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "AppsGeometryObject.hh"


class DiagDetector 
{


	public:
		DiagDetector(G4int id);
		~DiagDetector();

		AppsGeometryObject object;
		

		void ConstructMaterials();
		G4AssemblyVolume*  GetDetectorAssembly();
		void ReadInputParameters();		
		void ConstructCrystal();
		void ConstructDeadLayer();
		void ConstructInsideDeadLayer();
		void ConstructHolder();
		void ConstructAlWinEndcap();
		void ConstructDetectorCrystal(G4String name);
		
		void AddActiveVolumeName(G4String active){ ActiveVolumeNames.push_back(active); }
		vector<G4String> GetActiveVolumes() { return ActiveVolumeNames; } 

		// the FWHM parameters for the resolution

		void SetFWHMParameters(G4double slope, G4double intercept);	
		vector<G4double> GetFWHMParameters() { return FWHMParameters; }
		


	private:

		G4int objectId; 

		vector<G4double> FWHMParameters;
		vector<G4String> ActiveVolumeNames;

		G4Material* Ge_mat;
  		G4Material* Al_mat;

		G4double CrystalLength;
  		G4double CrystalRadius;
   		G4double FullCrystalLength;
 		G4double CrystalRoundRadius;
   		G4double WellRoundRadius;
   		G4double WellRadius;
   		G4double FrontDeadLayer;
   		G4double LateralDeadLayer;  
   		G4double BackDeadLayerThickness;
   		G4double WellLateralDeadLayer;
   		G4double WellBottomDeadLayer;
   		G4double EndcapInsideRadius;
   		G4double EndcapOutsideRadius;
   		G4double EndcapLength;
   		G4double BeWindowThickness;
   		G4double FrontToGeGap;
   		G4double EndcapRoundRadius; 
   		G4double HolderInsideRadius; 
   		G4double HolderBackInsideRadius; 
   		G4double HolderOutsideRadius; 
   		G4double HolderDiskRadius; 
   		G4double HolderBackOutsideRadius; 
   		G4double HolderLength; 
   		G4double FrontDiskLength;  
   		G4double DiskLength; 
   		G4double FrontToDiskOne; 
   		G4double FrontToDiskTwo; 
   		G4double HolderBackWallThickness; 
   		G4double HolderBackLength; 
   		G4double HolderBackCaseThickness;
 		G4double BackRiftDepth;
 		G4double BackRiftInsideRadius;
 		G4double BackRiftOutsideRadius;     		
		G4double Err;



		G4LogicalVolume* LogicalVolume1;
  		G4LogicalVolume* LogicalVolume2;
  		G4LogicalVolume* LogicalVolume3; 
   		G4LogicalVolume* LogicalVolume4;
		G4LogicalVolume* LogicalVolume5;		
		G4LogicalVolume* LogicalVolume6;

		// Objects from crystal 
		G4VSolid* GeCrystalHead1; 
		G4VSolid* GuardRing1; 
 		G4VSolid* InsideHead2;
   		G4VSolid * GeCrystal; 

};



