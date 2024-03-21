#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "AppsGeometryObject.hh"


class ScintDetector 
{


	public:
		ScintDetector(G4int id);
		~ScintDetector();

		AppsGeometryObject object;
		

		void ConstructMaterials();
		G4AssemblyVolume*  GetDetectorAssembly();
		void ReadInputParameters();		
		void ConstructCrystal(G4String name);
		void ConstructHolder();
		void ConstructQuartzWindow();

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
  		G4Material* NaI_mat;

		G4Material* Quartz_mat;


		G4double CrystalLength;
  		G4double CrystalRadius;

		G4double HolderLength;
		G4double HoldeRadius;

		G4double HoldeInnerRadius;
		G4double HoldeOuterRadius;
		G4double HolderFrontThickness;

		G4double FirstRingOuterRadius;
		G4double FirstRingInnerRadius;
		G4double FirstRingThickness;


		G4double SecondRingOuterRadius;
		G4double SecondRingInnerRadius;
		G4double SecondRingThickness;

		G4double QuartzWindowRadius;
 		G4double QuartzWindowLength;
		G4double FirstPhotoRingOuterRadius;
		G4double FirstPhotoRingInnerRadius;
		G4double FirstPhotoRingLength;
		G4double SecondPhotoRingInnerRadius;
		G4double SecondPhotoRingOuterRadius;
		G4double SecondPhotoRingLength;
  
		G4double SecondPhotoRingBackThickness;
		G4double CrystalToHolderGap;
		G4double ConeLength;

		G4double Err;



		G4LogicalVolume* LogicalVolume1;
		G4LogicalVolume* LogicalVolume2;
		G4LogicalVolume* LogicalVolume3;


};



