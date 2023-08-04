#include "DiagDetector.hh"
#include "G4NistManager.hh"
#include <fstream>
#include "G4SystemOfUnits.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "AppsInput.hh"


using namespace std;	

DiagDetector::DiagDetector(G4int id){


	AppsInput* FInput = AppsInput::Instance();

	object = FInput->GetGeometryObjects()[id];

	objectId = id;

	ReadInputParameters();
	ConstructMaterials();	

	// set linear FWHM parameters

	SetFWHMParameters(0.000228,0.000731);

}	



DiagDetector::~DiagDetector(){
}

void DiagDetector::ConstructCrystal(){

	G4VSolid* GeTorus1 = new G4Torus ("GeTorus1", 0.*mm, CrystalRoundRadius-FrontDeadLayer,CrystalRadius-CrystalRoundRadius,0.*deg,360.*deg);
	G4VSolid* GeCylinder1 = new G4Tubs ("GeCylinder1", 0*mm, CrystalRadius-CrystalRoundRadius, CrystalRoundRadius-FrontDeadLayer, 0.*deg, 360.*deg);
	GeCrystalHead1 = new G4UnionSolid ("GeCrystalHead1", GeTorus1, GeCylinder1); 
	G4VSolid* GeCylinder2 = new G4Tubs ("GeCylinder2",0*cm,CrystalRadius-LateralDeadLayer,(CrystalLength-CrystalRoundRadius-BackDeadLayerThickness)/2, 0.*deg, 360.*deg);
	GeCrystal= new G4UnionSolid ("GeCrystal",GeCrystalHead1, GeCylinder2,0, G4ThreeVector(0.,0.,(CrystalLength-CrystalRoundRadius-BackDeadLayerThickness)/2+Err)); 
  
}

void DiagDetector::ConstructDeadLayer() {

	G4VSolid* DeadTorus1 = new G4Torus ("DeadTorus1", 0.*mm, CrystalRoundRadius,CrystalRadius-CrystalRoundRadius,0.*deg,360.*deg);
     	G4VSolid* DeadCylinder1 = new G4Tubs ("DeadCylinder1",0*mm, CrystalRadius-CrystalRoundRadius, CrystalRoundRadius, 0.*deg, 360.*deg);
  	G4VSolid* DeadHead1 = new G4UnionSolid ("DeadHead1", DeadTorus1, DeadCylinder1); 
     	G4VSolid* DeadCylinder2 = new G4Tubs ("DeadCylinder2",0*cm,CrystalRadius+Err,CrystalLength/2, 0.*deg, 360.*deg);
     	G4VSolid* DeadHead2 = new G4SubtractionSolid ("DeadHead2", DeadHead1, DeadCylinder2,0,G4ThreeVector(0,0,CrystalLength/2));
     	G4VSolid* DeadHead3 = new G4SubtractionSolid ("DeadHead3", DeadHead2, GeCrystalHead1,0,G4ThreeVector(0,0,0));
     	G4VSolid* DeadCylinder3 = new G4Tubs ("DeadCylinder3",CrystalRadius-LateralDeadLayer,CrystalRadius,(CrystalLength-CrystalRoundRadius-BackDeadLayerThickness)/2,0.*deg,360.*deg);
     	G4VSolid* DeadLayer1 = new G4UnionSolid ("DeadLayer1",DeadHead3,DeadCylinder3,0,G4ThreeVector(0,0,(CrystalLength-CrystalRoundRadius-BackDeadLayerThickness)/2+Err)); 
     	G4VSolid* BackDeadLayer1 = new G4Tubs ("BackDeadLayer1",WellRadius+Err,CrystalRadius+Err,BackDeadLayerThickness/2, 0.*deg, 360.*deg);
	GuardRing1 = new G4Tubs ("GuardRing1",BackRiftInsideRadius, BackRiftOutsideRadius , BackRiftDepth/2, 0.*deg, 360.*deg); 
	G4VSolid* BackDeadLayer2 = new G4SubtractionSolid ("BackDeadLayer2", BackDeadLayer1, GuardRing1,0,G4ThreeVector(0,0,0));
   
  	LogicalVolume2 =   new G4LogicalVolume(DeadLayer1, Ge_mat, "DeadLayer");	
	LogicalVolume4 =   new G4LogicalVolume(BackDeadLayer2, Ge_mat, "BackDeadLayer");	


	G4VisAttributes* vis = new G4VisAttributes(false);
	LogicalVolume2->SetVisAttributes(vis);
	LogicalVolume4->SetVisAttributes(vis);

}


void DiagDetector::ConstructInsideDeadLayer() {

	G4VSolid* InsideTorus1 = new G4Torus ("InsideTorus1", 0.*mm, WellRoundRadius,WellRadius-WellRoundRadius+WellBottomDeadLayer,0.*deg,360.*deg);
    	G4VSolid* InsideCylinder1 = new G4Tubs ("InsideCylinder1",0*mm, WellRadius-WellRoundRadius+WellBottomDeadLayer,WellRoundRadius, 0.*deg, 360.*deg);
     	G4VSolid* InsideHead1  = new G4UnionSolid ("InsideHead1", InsideTorus1, InsideCylinder1); 
     	G4VSolid* InsideCylinder2 = new G4Tubs ("InsideCylinder2",0*cm,CrystalRadius+Err,CrystalLength/2,0.*deg,360.*deg);
 	InsideHead2 = new G4SubtractionSolid ("InsideHead2",InsideHead1,InsideCylinder2,0,G4ThreeVector(0,0,CrystalLength/2));
  
     	G4VSolid* InsideTorus2 = new G4Torus ("InsideTorus2", 0.*mm, WellRoundRadius-WellBottomDeadLayer,WellRadius-WellRoundRadius+WellBottomDeadLayer,0.*deg,360.*deg);
    	G4VSolid* InsideCylinder3 = new G4Tubs ("InsideCylinder3",0*mm, WellRadius-WellRoundRadius+WellBottomDeadLayer, WellRoundRadius-WellBottomDeadLayer, 0.*deg, 360.*deg);
     	G4VSolid* InsideHead3 = new G4UnionSolid ("InsideHead3", InsideTorus2, InsideCylinder3);  	
     	G4VSolid* InsideHead4 = new G4SubtractionSolid ("InsideHead4", InsideHead2, InsideHead3,0,G4ThreeVector(0,0,0));
     	G4VSolid* InsideCylinder4 = new G4Tubs ("InsideCylinder4",WellRadius,WellRadius+WellLateralDeadLayer,(CrystalLength-FullCrystalLength-BackDeadLayerThickness
			-WellRoundRadius+WellBottomDeadLayer)/2,0.*deg,360.*deg);
     	G4VSolid* InsideDeadlayer1 = new G4UnionSolid ("InsideDeadlayer1",InsideHead4,InsideCylinder4,0,G4ThreeVector(0,0,(CrystalLength-FullCrystalLength
			-BackDeadLayerThickness-WellRoundRadius+WellBottomDeadLayer)/2));
 
  	LogicalVolume3 =   new G4LogicalVolume(InsideDeadlayer1, Ge_mat, "InsideDeadLayer");	

	G4VisAttributes* vis = new G4VisAttributes(false);
	LogicalVolume3->SetVisAttributes(vis);

}

void DiagDetector::ConstructHolder(){
 
     	G4VSolid* CrystalEndcapCylinder1 = new G4Tubs ("CrystalEndcapCylinder1",HolderInsideRadius, HolderOutsideRadius , HolderLength/2, 0.*deg, 360.*deg);
     	G4VSolid* CrystalEndcapCylinder2 = new G4Tubs ("CrystalEndcapCylinder2",HolderOutsideRadius,HolderDiskRadius , FrontDiskLength/2, 0.*deg, 360.*deg);
     	G4VSolid* CrystalHolder1 = new G4UnionSolid ("CrystalHolder1",CrystalEndcapCylinder2,CrystalEndcapCylinder1,0,G4ThreeVector(0,0,(HolderLength/2-FrontDiskLength/2)));
     	G4VSolid* CrystalEndcapCylinder3 = new G4Tubs ("CrystalEndcapCylinder3",HolderOutsideRadius,HolderDiskRadius ,DiskLength/2 , 0.*deg, 360.*deg); 
     	G4VSolid* CrystalHolder2 = new G4UnionSolid ("CrystalHolder2",CrystalHolder1,CrystalEndcapCylinder3,0,G4ThreeVector(0,0,(DiskLength/2-FrontDiskLength/2+FrontToDiskOne)));
     	G4VSolid* CrystalHolder3 = new G4UnionSolid ("CrystalHolder3",CrystalHolder2,CrystalEndcapCylinder3,0,G4ThreeVector(0,0,(DiskLength/2-FrontDiskLength/2+FrontToDiskTwo)));      
     	G4VSolid* CrystalEndcapCylinder4 = new G4Tubs ("CrystalEndcapCylinder4",HolderBackInsideRadius,HolderOutsideRadius ,HolderBackWallThickness/2 , 0.*deg, 360.*deg); 
     	G4VSolid* CrystalHolder4 = new G4UnionSolid ("CrystalHolder4",CrystalHolder3,CrystalEndcapCylinder4,0,G4ThreeVector(0,0,(HolderBackWallThickness/2-
					FrontDiskLength/2+HolderLength)));   
     	G4VSolid* CrystalEndcapCylinder5 = new G4Tubs ("CrystalEndcapCylinder5",HolderBackInsideRadius,HolderBackOutsideRadius,
					(HolderBackLength-HolderBackCaseThickness)/2 , 0.*deg, 360.*deg); 
     	G4VSolid* CrystalHolderCylinder6 = new G4Tubs ("CrystalHolderCylinder6",0,HolderBackOutsideRadius ,HolderBackCaseThickness/2 , 0.*deg, 360.*deg); 
     	G4VSolid* CrystalHolderBack1 = new G4UnionSolid ("CrystalHolderBack1",CrystalEndcapCylinder5,CrystalHolderCylinder6,0,G4ThreeVector(0,0,((HolderBackLength-			
					HolderBackCaseThickness)/2+HolderBackCaseThickness/2)));   
     	G4VSolid* CrystalHolder5 = new G4UnionSolid ("CrystalHolder5",CrystalHolder4,CrystalHolderBack1,0,G4ThreeVector(0,0,((HolderBackLength-HolderBackCaseThickness)/2-
					FrontDiskLength/2+HolderLength+HolderBackWallThickness)));    


	

  	LogicalVolume5 =   new G4LogicalVolume(CrystalHolder5, Al_mat, "Holder");	

}

void DiagDetector::ConstructAlWinEndcap(){

      	G4VSolid* EndcapTorus1 = new G4Torus ("EndcapTorus1", 0.*mm, EndcapRoundRadius,EndcapOutsideRadius-EndcapRoundRadius,0.*deg,360.*deg);
      	G4VSolid* EndcapCylinder1 = new G4Tubs ("EndcapCylinder1",0, EndcapOutsideRadius-EndcapRoundRadius,EndcapRoundRadius, 0.*deg, 360.*deg);
   	G4VSolid* EndcapHead1 = new G4UnionSolid ("EndcapHead1", EndcapTorus1, EndcapCylinder1); 
   	G4VSolid* EndcapCylinder2 = new G4Tubs ("EndcapCylinder2",0*cm,EndcapOutsideRadius+Err,CrystalLength/2,0.*deg,360.*deg);
   	G4VSolid* EndcapHead2 = new G4SubtractionSolid ("EndcapHead1",EndcapHead1,EndcapCylinder2,0,G4ThreeVector(0,0,CrystalLength/2));
    	G4VSolid* EndcapTorus2 = new G4Torus ("EndcapTorus2", 0.*mm, EndcapRoundRadius-(EndcapOutsideRadius-EndcapInsideRadius),EndcapOutsideRadius-EndcapRoundRadius,0.*deg,360.*deg);
      	G4VSolid* EndcapCylinder3 = new G4Tubs ("EndcapCylinder3",0, EndcapOutsideRadius-EndcapRoundRadius, EndcapRoundRadius-
					(EndcapOutsideRadius-EndcapInsideRadius), 0.*deg, 360.*deg);
      	G4VSolid* EndcapHead3 = new G4UnionSolid ("EndcapHead3",EndcapTorus2, EndcapCylinder3);
      	G4VSolid* EndcapHead4 = new G4SubtractionSolid ("EndcapHead4",EndcapHead2,EndcapHead3,0,G4ThreeVector(0,0,0));
      	G4VSolid* EndcapCylinder5 = new G4Tubs ("EndcapCylinder5",EndcapInsideRadius, EndcapOutsideRadius , (EndcapLength-EndcapRoundRadius)/2, 0.*deg, 360.*deg);
      	G4VSolid* Endcap1 = new G4UnionSolid ("Endcap1", EndcapHead4, EndcapCylinder5, 0, G4ThreeVector(0,0,(EndcapLength-EndcapRoundRadius)/2));
      	G4VSolid* Endcap2 = new G4Tubs ("Endcap2",0*cm,EndcapOutsideRadius,(EndcapOutsideRadius-EndcapInsideRadius)/2,0.*deg,360.*deg);
      	G4VSolid* Endcap3 = new G4UnionSolid ("Endcap3", Endcap1, Endcap2, 0, G4ThreeVector(0,0,EndcapLength-EndcapRoundRadius));

	LogicalVolume6 =   new G4LogicalVolume(Endcap3, Al_mat, "Endcap");	

}

void DiagDetector::ConstructDetectorCrystal(G4String name){
	
   	G4VSolid* InsideCylinder5 = new G4Tubs ("InsideCylinder5",0,WellRadius+WellLateralDeadLayer,(CrystalLength-FullCrystalLength-BackDeadLayerThickness-
				WellRoundRadius+WellBottomDeadLayer)/2,0.*deg,360.*deg);
     	G4VSolid* InsideDeadlayer2 = new G4UnionSolid ("InsideDeadlayer2",InsideHead2,InsideCylinder5,0,G4ThreeVector(0,0,(CrystalLength-FullCrystalLength-
				BackDeadLayerThickness-WellRoundRadius+WellBottomDeadLayer)/2)); 
     	G4VSolid* GeCrystal2 = new G4SubtractionSolid ("GeCrystal2",GeCrystal,InsideDeadlayer2,0,G4ThreeVector(0,0,FullCrystalLength-CrystalRoundRadius+WellRoundRadius-
				WellBottomDeadLayer+2*Err));
	G4VSolid* GeCrystal3 = new G4SubtractionSolid ("GeCrystal3",GeCrystal2,GuardRing1,0,G4ThreeVector(0,0,CrystalLength-CrystalRoundRadius
				-BackDeadLayerThickness+6*Err-BackRiftDepth/2));
 
  	LogicalVolume1 =   new G4LogicalVolume(GeCrystal3, Ge_mat, name);	
 
} 


void DiagDetector::ReadInputParameters(){
	

	G4String z; 
	G4double a[33];

	ifstream myfile ("DiagDetParameters.txt");
	if(myfile.is_open()){
		getline(myfile, z);
		getline(myfile, z);
   		getline(myfile, z);
   
		for(G4int i=0;i<33;i++){
   
			myfile >> a[i];
  	 		getline(myfile, z);
  			getline(myfile, z); 

		}
	} else {
   		cout << "Unable to open DiagDetParameters.txt" << G4endl;
	}
   
 
  	CrystalLength=a[0]*cm;
  	CrystalRadius=a[1]*cm;
  	FullCrystalLength=a[2]*cm;
  	CrystalRoundRadius=a[3]*cm;
  	WellRoundRadius=a[4]*cm;
 	WellRadius=a[5]*cm;
 	FrontDeadLayer=a[6]*cm;
  	LateralDeadLayer=a[7]*cm;  
  	BackDeadLayerThickness=a[8]*cm;
  	WellLateralDeadLayer=a[9]*cm;
  	WellBottomDeadLayer=a[10]*cm;
  	EndcapInsideRadius= a[11]*cm;
  	EndcapOutsideRadius=a[12]*cm;
  	EndcapLength= a[13]*cm;
  	BeWindowThickness=a[14]*cm;
  	FrontToGeGap=a[15]*cm;
  	EndcapRoundRadius=a[16]*cm; 
  	HolderInsideRadius= a[17]*cm; 
  	HolderBackInsideRadius= a[18]*cm; 
 	HolderOutsideRadius= a[19]*cm; 
  	HolderDiskRadius= a[20]*cm;   
  	HolderBackOutsideRadius=a[21]*cm; 
  	HolderLength=a[22]*cm; 
  	FrontDiskLength=a[23]*cm;  
  	DiskLength=a[24]*cm; 
  	FrontToDiskOne=a[25]*cm; 
  	FrontToDiskTwo=a[26]*cm; 
  	HolderBackWallThickness=a[27]*cm; 
  	HolderBackLength=a[28]*cm; 
  	HolderBackCaseThickness=a[29]*cm; 
    	BackRiftDepth = a[30]*cm;
  	BackRiftInsideRadius = a[31]*cm;
 	BackRiftOutsideRadius= a[32]*cm;
	Err=0.010*cm; 

}


void DiagDetector::ConstructMaterials(){

	G4NistManager* nist = G4NistManager::Instance();

	Ge_mat = nist->FindOrBuildMaterial("G4_Ge");
	Al_mat = nist->FindOrBuildMaterial("G4_Al");

}


G4AssemblyVolume*  DiagDetector::GetDetectorAssembly(){
 

	ConstructCrystal();
	ConstructDeadLayer();
	ConstructInsideDeadLayer();
	ConstructHolder();
	ConstructAlWinEndcap();
	ConstructDetectorCrystal(object.GetObjectName());



	G4AssemblyVolume* assemblyDetector1 = new G4AssemblyVolume();

	G4ThreeVector pos1 = G4ThreeVector(0,0,CrystalRoundRadius+FrontToGeGap+EndcapRoundRadius);
	G4ThreeVector pos2 = G4ThreeVector(0,0,FrontToGeGap+EndcapRoundRadius+CrystalLength-BackDeadLayerThickness/2+2*Err);
	G4ThreeVector pos3 = G4ThreeVector(0,0,FrontToGeGap+EndcapRoundRadius+FullCrystalLength+WellRoundRadius-WellBottomDeadLayer+2*Err);
	G4ThreeVector pos4 = G4ThreeVector(0,0,FrontToGeGap+EndcapRoundRadius);
	G4ThreeVector pos5 = G4ThreeVector(0,0,EndcapRoundRadius+5*Err);

	G4RotationMatrix *  yRot2 = new G4RotationMatrix;
	yRot2->rotateY(0*deg);

     
  	assemblyDetector1->AddPlacedVolume(LogicalVolume2,pos1,yRot2);
  	assemblyDetector1->AddPlacedVolume(LogicalVolume4,pos2,yRot2);
  	assemblyDetector1->AddPlacedVolume(LogicalVolume3,pos3,yRot2);
  	assemblyDetector1->AddPlacedVolume(LogicalVolume5,pos4,yRot2);
  	assemblyDetector1->AddPlacedVolume(LogicalVolume6,pos5,yRot2);
  	assemblyDetector1->AddPlacedVolume(LogicalVolume1,pos1,yRot2);



	AddActiveVolumeName(object.GetObjectName());


	return assemblyDetector1;

}


void DiagDetector::SetFWHMParameters(G4double slope, G4double intercept){
	
	FWHMParameters.push_back(slope);
	FWHMParameters.push_back(intercept);

}	


