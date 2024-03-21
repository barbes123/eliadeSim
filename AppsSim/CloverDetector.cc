#include "CloverDetector.hh"
#include "G4NistManager.hh"
#include <fstream>
#include "G4SystemOfUnits.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VisAttributes.hh"




using namespace std;


CloverDetector::CloverDetector(G4int id){

	AppsInput* FInput = AppsInput::Instance();

	AppsGeometryObject object = FInput->GetGeometryObjects()[id];
	
	objectId = id;

	ReadInputParameters();
	ConstructMaterials();	
	
	ConstructCrystal();
	ConstructDetectorCrystal();
	ConstructInsideDeadLayer();
	ConstructDeadLayer();
	ConstructBackDeadLayer();	
	
	ConstructHolder(object.GetObjectName());
	ConstructBackHolder();
	ConstructCrystalSegments();
	
	

	for(unsigned int i = 0; i < object.GetObjectFeatures().size(); i++){
		if(object.GetObjectFeatures()[i] == "segmented"){
			ConstructCrystalSegments();
			continue;
		}

		if(object.GetObjectFeatures()[i] == "shielded"){
			ConstructShield();
			continue;
		}
	
	} 

	
	// set linear FWHM parameters

	SetFWHMParameters(0.000228,0.000731);



}

CloverDetector::~CloverDetector(){
}


G4AssemblyVolume* CloverDetector::ConstructSegmCloverAssembly(G4String name){

	G4AssemblyVolume* assemblySegClover = new G4AssemblyVolume();

	
	G4ThreeVector pos0 = G4ThreeVector(0,0,0);
	G4ThreeVector pos1 = G4ThreeVector(0,0,HolderFrontWallThickness+FrontDeadLayer+HolderToCrystalDistance);
  	G4ThreeVector pos2 = G4ThreeVector(0,0,HolderFrontLength/2);
 	G4ThreeVector pos3 = G4ThreeVector(0,0,HolderFrontLength/2-HolderFrontWallThickness-HolderToCrystalDistance-CrystalLength/2-3*err+CrystalLength/2);
	G4ThreeVector pos4 = G4ThreeVector(0,0,FrontToBackHolderDistance);
 

	G4RotationMatrix *  yRot1 = new G4RotationMatrix;
	yRot1->rotateY(0*deg);

	G4RotationMatrix *  yRot2 = new G4RotationMatrix;
	yRot2->rotateZ(270*deg);

	G4RotationMatrix *  yRot3 = new G4RotationMatrix;
	yRot3->rotateZ(180*deg);

	G4RotationMatrix *  yRot4 = new G4RotationMatrix;
	yRot4->rotateZ(90*deg);

	G4RotationMatrix *  yRot5 = new G4RotationMatrix;
	yRot5->rotateY(90*deg);

	if(!logicShape7){			// shielded condition



		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal0"),pos1,yRot1);
		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal1"),pos1,yRot2);
		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal2"),pos1,yRot3);
		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal3"),pos1,yRot4);
		
		assemblySegClover->AddPlacedVolume(logicShape6,pos2,yRot2);
		assemblySegClover->AddPlacedVolume(logicShape5,pos4,yRot2);

	} else {

		G4double shieldTodet = 10*cm;
		
		G4ThreeVector pos5 = G4ThreeVector(0,0,-shieldTodet);


		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal0"),pos1,yRot1);
		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal1"),pos1,yRot2);
		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal2"),pos1,yRot3);
		assemblySegClover->AddPlacedAssembly(ConstructSegCrystalAssembly(name+"_Crystal3"),pos1,yRot4);
		
		assemblySegClover->AddPlacedVolume(logicShape6,pos2,yRot2);
		assemblySegClover->AddPlacedVolume(logicShape5,pos4,yRot2);

		assemblySegClover->AddPlacedAssembly(ConstructShieldAssembly(name+"_Shield"),pos0,yRot1);

		bool isVisible = true;
 
		double Red = 0;
		double Green = 97;	
		double Blue = 128;

		G4VisAttributes * tempVisAttr = new G4VisAttributes(isVisible,
			G4Colour(Red / 255., Green / 255., Blue / 255.));

		logicShape7->SetVisAttributes(tempVisAttr);

	}


	return assemblySegClover;
}


G4AssemblyVolume* CloverDetector::ConstructCloverAssembly(G4String name){


	

	G4AssemblyVolume* assemblyClover = new G4AssemblyVolume();


	G4ThreeVector pos1 = G4ThreeVector(0,0,HolderFrontWallThickness+FrontDeadLayer+HolderToCrystalDistance);
	G4ThreeVector pos2 = G4ThreeVector(0,0,HolderFrontLength/2);
 	G4ThreeVector pos3 = G4ThreeVector(0,0,HolderFrontLength/2-HolderFrontWallThickness-HolderToCrystalDistance-CrystalLength/2-3*err+CrystalLength/2);
	G4ThreeVector pos4 = G4ThreeVector(0,0,FrontToBackHolderDistance);


	G4RotationMatrix *  yRot1 = new G4RotationMatrix;
	yRot1->rotateY(0*deg);

	G4RotationMatrix *  yRot2 = new G4RotationMatrix;
	yRot2->rotateZ(90*deg);

	G4RotationMatrix *  yRot3 = new G4RotationMatrix;
	yRot3->rotateZ(180*deg);

	G4RotationMatrix *  yRot4 = new G4RotationMatrix;
	yRot4->rotateZ(270*deg);

	if(!logicShape7){			// shielded condition



		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal0"),pos1,yRot1);
		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal1"),pos1,yRot2);
		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal2"),pos1,yRot3);
		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal3"),pos1,yRot4);
		
		assemblyClover->AddPlacedVolume(logicShape6,pos2,yRot2);
		assemblyClover->AddPlacedVolume(logicShape5,pos4,yRot2);

	} else {
		
		G4double shieldTodet = 1*cm;
		
		G4ThreeVector pos5 = G4ThreeVector(0,0,-shieldTodet);


		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal0"),pos1,yRot1);
		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal1"),pos1,yRot2);
		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal2"),pos1,yRot3);
		assemblyClover->AddPlacedAssembly(ConstructCrystalAssembly(name+"_Crystal3"),pos1,yRot4);
		
		assemblyClover->AddPlacedVolume(logicShape6,pos2,yRot2);
		assemblyClover->AddPlacedVolume(logicShape5,pos4,yRot2);




		bool isVisible = true;
 
		double Red = 0;
		double Green = 97;	
		double Blue = 128;

		G4VisAttributes * tempVisAttr = new G4VisAttributes(isVisible,
			G4Colour(Red / 255., Green / 255., Blue / 255.));

		logicShape7->SetVisAttributes(tempVisAttr);

		assemblyClover->AddPlacedVolume(logicShape7,pos5,yRot1);	
		
	
	}


	return assemblyClover;
}


G4AssemblyVolume* CloverDetector::ConstructShieldAssembly(G4String name){

	G4AssemblyVolume* shieldAssembly = new G4AssemblyVolume();
	
	FrontShield1 = new G4LogicalVolume(FrontShield, BGO_mat, name + "_FrontShield1");
	FrontShield2 = new G4LogicalVolume(FrontShield, BGO_mat, name + "_FrontShield2");
	FrontShield3 = new G4LogicalVolume(FrontShield, BGO_mat, name + "_FrontShield3");
	FrontShield4 = new G4LogicalVolume(FrontShield, BGO_mat, name + "_FrontShield4");
	SideShield1 = new G4LogicalVolume(SideShield, BGO_mat, name + "_SideShield1");
	SideShield2 = new G4LogicalVolume(SideShield, BGO_mat, name + "_SideShield2");
	SideShield3 = new G4LogicalVolume(SideShield, BGO_mat, name + "_SideShield3");
	SideShield4 = new G4LogicalVolume(SideShield, BGO_mat, name + "_SideShield4");
	BackShield1 = new G4LogicalVolume(BackShield, CsI_mat, name + "_BackShield1");
	BackShield2 = new G4LogicalVolume(BackShield, CsI_mat, name + "_BackShield2");
	FrontShieldb1 = new G4LogicalVolume(FrontShieldb, Air_mat, name + "_FrontShieldb1");
	FrontShieldb2 = new G4LogicalVolume(FrontShieldb, Air_mat, name + "_FrontShieldb2");
	FrontShieldb3 = new G4LogicalVolume(FrontShieldb, Air_mat, name + "_FrontShieldb3");
	FrontShieldb4 = new G4LogicalVolume(FrontShieldb, Air_mat, name + "_FrontShieldb4");
	SideShieldb1 = new G4LogicalVolume(SideShieldb, Air_mat, name + "_SideShield1");
	SideShieldb2 = new G4LogicalVolume(SideShieldb, Air_mat, name + "_SideShield2");
	SideShieldb3 = new G4LogicalVolume(SideShieldb, Air_mat, name + "_SideShield3");
	SideShieldb4 = new G4LogicalVolume(SideShieldb, Air_mat, name + "_SideShield4");
	BackShieldb1 = new G4LogicalVolume(BackShieldb, Air_mat, name + "_BackShieldb1");
	BackShieldb2 = new G4LogicalVolume(BackShieldb, Air_mat, name + "_BackShieldb2");
	FrontShieldc1 = new G4LogicalVolume(FrontShieldc, Pb_mat, name + "_FrontShieldc1");
	FrontShieldc2 = new G4LogicalVolume(FrontShieldc, Pb_mat, name + "_FrontShieldc2");
	FrontShieldc3 = new G4LogicalVolume(FrontShieldc, Pb_mat, name + "_FrontShieldc3");
	FrontShieldc4 = new G4LogicalVolume(FrontShieldc, Pb_mat, name + "_FrontShieldc4");

	G4RotationMatrix *  yRot1 = new G4RotationMatrix;
	yRot1->rotateY(90.*deg);
	G4RotationMatrix *  yRot2 = new G4RotationMatrix;
	yRot2->rotateZ(90.*deg);
	yRot2->rotateY(22.5*deg);
	G4RotationMatrix *  yRot3 = new G4RotationMatrix;
	yRot3->rotateY(-90.*deg);
	yRot3->rotateX(180.*deg);
	
	sideShieldZ = HolderFrontLength+sideShieldaH+2*err+0.5*cm;
	frontShieldZ = frontShieldaH/2.*cos(22.5*deg)-frontShieldaT/2.*sin(22.5*deg)-0.5*cm;
	
	G4ThreeVector pos1 = G4ThreeVector(0., sideShieldaL/2.+err, sideShieldZ);
	G4ThreeVector pos1a = G4ThreeVector(0., sideShieldaL/2.+err+0.5*mm, sideShieldZ-sideShieldaH+sideShieldaaH+2.*mm);
	G4ThreeVector pos2 = G4ThreeVector((frontShieldaL+frontShieldal)/4.-frontShieldaT/2.+50*err, 0., frontShieldZ);
	G4ThreeVector pos2c = G4ThreeVector((frontShieldcl+frontShieldcL)/4.-frontShieldcT/2.+50*err, 0., frontShieldZ-(frontShieldaH+frontShieldcH)/2.*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/4.*sin(22.5*deg)+err);
	G4ThreeVector pos3 = G4ThreeVector(-backShieldL/4.-err, 0., HolderLength+backShieldH/2.+err+1.*mm);
	G4ThreeVector pos3a = G4ThreeVector(-backShieldL/4.-err, 0., HolderLength+backShieldH/2.+err+1.*mm-(backShieldH-backShieldaH)/2.+backShieldalT);
	G4ThreeVector pos11 = G4ThreeVector(frontShieldbL/2.-frontShieldbT/2., 0., frontShieldH-(frontShieldaH+2*frontShieldcH)/2.*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/2*sin(22.5*deg)-frontShieldbH/2.+frontShieldZ);

	shieldAssembly->AddPlacedVolume(SideShieldb1,pos1,yRot1);
	shieldAssembly->AddPlacedVolume(SideShield1,pos1a,yRot1);
	shieldAssembly->AddPlacedVolume(FrontShield1,pos2,yRot2);
	shieldAssembly->AddPlacedVolume(FrontShieldc1,pos2c,yRot2);
	shieldAssembly->AddPlacedVolume(BackShield1,pos3a,yRot1);
	shieldAssembly->AddPlacedVolume(BackShieldb1,pos3,yRot1);
	shieldAssembly->AddPlacedVolume(FrontShieldb1,pos11,yRot3);
	
	yRot1->rotateZ(90.*deg);
	yRot2->rotateZ(90.*deg);
	yRot3->rotateZ(90.*deg);

	G4ThreeVector pos4 = G4ThreeVector(-sideShieldaL/2.-err, 0., sideShieldZ);
	G4ThreeVector pos4a = G4ThreeVector(-sideShieldaL/2.-err-0.5*mm, 0., sideShieldZ-sideShieldaH+sideShieldaaH+2.*mm);
	G4ThreeVector pos5 = G4ThreeVector(0., (frontShieldaL+frontShieldal)/4.-frontShieldaT/2.+50*err, frontShieldZ);
	G4ThreeVector pos5c = G4ThreeVector(0., (frontShieldcl+frontShieldcL)/4.-frontShieldcT/2.+50*err, frontShieldZ-(frontShieldaH+frontShieldcH)/2.*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/4.*sin(22.5*deg)+err);
	G4ThreeVector pos12 = G4ThreeVector(0., frontShieldbL/2.-frontShieldbT/2., frontShieldH-(frontShieldaH+2*frontShieldcH)/2*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/2*sin(22.5*deg)-frontShieldbH/2.+frontShieldZ);

	shieldAssembly->AddPlacedVolume(SideShieldb2,pos4,yRot1);
	shieldAssembly->AddPlacedVolume(SideShield2,pos4a,yRot1);
	shieldAssembly->AddPlacedVolume(FrontShield2,pos5,yRot2);
	shieldAssembly->AddPlacedVolume(FrontShieldc2,pos5c,yRot2);
	shieldAssembly->AddPlacedVolume(FrontShieldb2,pos12,yRot3);
	
	yRot1->rotateZ(90.*deg);
	yRot2->rotateZ(90.*deg);
	yRot3->rotateZ(90.*deg);

	G4ThreeVector pos6 = G4ThreeVector(0., -sideShieldaL/2.-err, sideShieldZ);
	G4ThreeVector pos6a = G4ThreeVector(0., -sideShieldaL/2.-err-0.5*mm, sideShieldZ-sideShieldaH+sideShieldaaH+2.*mm);
	G4ThreeVector pos7 = G4ThreeVector(-(frontShieldaL+frontShieldal)/4.+frontShieldaT/2.-50*err, 0., frontShieldZ);
	G4ThreeVector pos7c = G4ThreeVector(-(frontShieldcl+frontShieldcL)/4.+frontShieldcT/2.-50*err, 0., frontShieldZ-(frontShieldaH+frontShieldcH)/2.*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/4.*sin(22.5*deg)+err);
	G4ThreeVector pos8 = G4ThreeVector(backShieldL/4.+err, 0., HolderLength+backShieldH/2.+err+1.*mm);
	G4ThreeVector pos8a = G4ThreeVector(backShieldL/4.+err, 0., HolderLength+backShieldH/2.+err+1.*mm-(backShieldH-backShieldaH)/2.+backShieldalT);
	G4ThreeVector pos13 = G4ThreeVector(-frontShieldbL/2.+frontShieldbT/2., 0., frontShieldH-(frontShieldaH+2*frontShieldcH)/2*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/2*sin(22.5*deg)-frontShieldbH/2.+frontShieldZ);

	shieldAssembly->AddPlacedVolume(SideShieldb3,pos6,yRot1);
	shieldAssembly->AddPlacedVolume(SideShield3,pos6a,yRot1);
	shieldAssembly->AddPlacedVolume(FrontShield3,pos7,yRot2);
	shieldAssembly->AddPlacedVolume(FrontShieldc3,pos7c,yRot2);
	shieldAssembly->AddPlacedVolume(BackShield2,pos8a,yRot1);
	shieldAssembly->AddPlacedVolume(BackShieldb2,pos8,yRot1);
	shieldAssembly->AddPlacedVolume(FrontShieldb3,pos13,yRot3);
	
	yRot1->rotateZ(90.*deg);
	yRot2->rotateZ(90.*deg);
	yRot3->rotateZ(90.*deg);

	G4ThreeVector pos9 = G4ThreeVector(sideShieldaL/2.+err,0, sideShieldZ);
	G4ThreeVector pos9a = G4ThreeVector(sideShieldaL/2.+err+0.5*mm,0, sideShieldZ-sideShieldaH+sideShieldaaH+2.*mm);
	G4ThreeVector pos10 = G4ThreeVector(0,-(frontShieldaL+frontShieldal)/4.+frontShieldaT/2.-50*err, frontShieldZ);
	G4ThreeVector pos10c = G4ThreeVector(0,-(frontShieldcl+frontShieldcL)/4.+frontShieldcT/2.-50*err, frontShieldZ-(frontShieldaH+frontShieldcH)/2.*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/4.*sin(22.5*deg)+err);
	G4ThreeVector pos14 = G4ThreeVector(0., -frontShieldbL/2.+frontShieldbT/2., frontShieldH-(frontShieldaH+2*frontShieldcH)/2*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/2*sin(22.5*deg)-frontShieldbH/2.+frontShieldZ);

	shieldAssembly->AddPlacedVolume(SideShieldb4,pos9,yRot1);
	shieldAssembly->AddPlacedVolume(SideShield4,pos9a,yRot1);
	shieldAssembly->AddPlacedVolume(FrontShield4,pos10,yRot2);
	shieldAssembly->AddPlacedVolume(FrontShieldc4,pos10c,yRot2);
	shieldAssembly->AddPlacedVolume(FrontShieldb4,pos14,yRot3);

	AddActiveVolumeName(name + "_SSh1");	
	AddActiveVolumeName(name + "_SSh2");	
	AddActiveVolumeName(name + "_SSh3");	
	AddActiveVolumeName(name + "_SSh4");	
	AddActiveVolumeName(name + "_FSh1");	
	AddActiveVolumeName(name + "_FSh2");	
	AddActiveVolumeName(name + "_FSh3");	
	AddActiveVolumeName(name + "_FSh4");	
	AddActiveVolumeName(name + "_BSh1");	
	AddActiveVolumeName(name + "_BSh2");
	
	AppsInput* FInput = AppsInput::Instance();
	

	AppsAnalysisObject * shianalysis = new AppsAnalysisObject();
	shianalysis->SetObjectName(name+"_Shield");


	for(int i = 0; i < 4; i++){

		vector<AppsAnalysisObject *> temp;		

		AppsAnalysisObject * analysis = new AppsAnalysisObject();

		analysis->SetObjectName(name + "_FrontShield" + to_string(i+1));
		analysis->SetDetectorResolution(0.000228,0.000731);		
		analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());

		temp.push_back(analysis);
		temp.push_back(shianalysis);

		FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);

	}

	for(int i = 0; i < 4; i++){

		vector<AppsAnalysisObject *> temp;		

		AppsAnalysisObject * analysis = new AppsAnalysisObject();

		analysis->SetObjectName(name + "_SideShield" + to_string(i+1));
		analysis->SetDetectorResolution(0.000228,0.000731);		
		analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());

		temp.push_back(analysis);
		temp.push_back(shianalysis);

		FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);

	}

	for(int i = 0; i < 2; i++){

		vector<AppsAnalysisObject *> temp;		

		AppsAnalysisObject * analysis = new AppsAnalysisObject();

		analysis->SetObjectName(name + "_BackShield" + to_string(i+1));
		analysis->SetDetectorResolution(0.000228,0.000731);		
		analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());

		temp.push_back(analysis);
		temp.push_back(shianalysis);

		FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);

	}

	return shieldAssembly;
	
}


G4AssemblyVolume* CloverDetector::ConstructSegCrystalAssembly(G4String name){

	G4AssemblyVolume* segmentsAssembly = new G4AssemblyVolume();

	logicShape2 = new G4LogicalVolume(GeCrystal7, Ge_mat, name);			// build the logical volume of the detector

	logicShape3 =   new G4LogicalVolume(InsideCylinder1, Ge_mat, name + "_DeadLayer");	// build the inside dead layer
	logicShape1 =   new G4LogicalVolume(DeadGeCrystal8, Ge_mat, name + "_DeadLayer");	// build the inside dead layer
	logicShape4 =   new G4LogicalVolume(BackDeadLayer5, Ge_mat, name + "_DeadLayer");	// build the back dead layer

	Segment1 =   new G4LogicalVolume(Seg3, Ge_mat, name + "_Segment1");	
  	Segment2 =   new G4LogicalVolume(Seg4, Ge_mat, name + "_Segment2");	
  	Segment3 =   new G4LogicalVolume(Seg1, Ge_mat, name + "_Segment3");	
  	Segment4 =   new G4LogicalVolume(Seg2, Ge_mat, name + "_Segment4");	
  	Segment5 =   new G4LogicalVolume(Seg7, Ge_mat, name + "_Segment5");	
  	Segment6 =   new G4LogicalVolume(Seg8, Ge_mat, name + "_Segment6");	
  	Segment7 =   new G4LogicalVolume(Seg5, Ge_mat, name + "_Segment7");	
  	Segment8 =   new G4LogicalVolume(Seg6, Ge_mat, name + "_Segment8");





  	G4ThreeVector pos1 = G4ThreeVector((CrystalInsideLength+err)+SpaceBetweenCrystals/2,(CrystalInsideLength+err)+SpaceBetweenCrystals/2,CrystalLength/2);
	G4ThreeVector pos2 = G4ThreeVector((CrystalInsideLength+err+SpaceBetweenCrystals/2),(CrystalInsideLength+err+SpaceBetweenCrystals/2),
				CrystalLength/2-BackDeadLayer/2+2*err+CrystalLength/2);	
	G4ThreeVector pos3 = G4ThreeVector((CrystalInsideLength+err+SpaceBetweenCrystals/2),(CrystalInsideLength+err+SpaceBetweenCrystals/2),	
				-((WellLength-BackDeadLayer))/2+CrystalLength/2-BackDeadLayer+CrystalLength/2);


	G4RotationMatrix *  yRot2 = new G4RotationMatrix;
	yRot2->rotateY(0*deg);

	segmentsAssembly->AddPlacedVolume(Segment1,pos1,yRot2);
	segmentsAssembly->AddPlacedVolume(Segment2,pos1,yRot2);	
	segmentsAssembly->AddPlacedVolume(Segment3,pos1,yRot2);	
	segmentsAssembly->AddPlacedVolume(Segment4,pos1,yRot2);
	segmentsAssembly->AddPlacedVolume(Segment5,pos1,yRot2);
	segmentsAssembly->AddPlacedVolume(Segment6,pos1,yRot2);	
	segmentsAssembly->AddPlacedVolume(Segment7,pos1,yRot2);	
	segmentsAssembly->AddPlacedVolume(Segment8,pos1,yRot2);	
	segmentsAssembly->AddPlacedVolume(logicShape1,pos1,yRot2);	
	segmentsAssembly->AddPlacedVolume(logicShape4,pos2,yRot2);	
	segmentsAssembly->AddPlacedVolume(logicShape3,pos3,yRot2);

		
	AddActiveVolumeName(name + "_Seg1");	
	AddActiveVolumeName(name + "_Seg2");	
	AddActiveVolumeName(name + "_Seg3");	
	AddActiveVolumeName(name + "_Seg4");	
	AddActiveVolumeName(name + "_Seg5");	
	AddActiveVolumeName(name + "_Seg6");	
	AddActiveVolumeName(name + "_Seg7");	
	AddActiveVolumeName(name + "_Seg8");	


	AppsInput* FInput = AppsInput::Instance();
	

	AppsAnalysisObject * cryanalysis = new AppsAnalysisObject();
	cryanalysis->SetObjectName(name);


	for(int i = 0; i < 8; i++){

		vector<AppsAnalysisObject *> temp;		

		AppsAnalysisObject * analysis = new AppsAnalysisObject();



		analysis->SetObjectName(name + "_Segment" + to_string(i+1));
		analysis->SetDetectorResolution(0.000228,0.000731);		
		analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());



		temp.push_back(analysis);
		temp.push_back(cryanalysis);

		FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);

	}


	return segmentsAssembly;

}



G4AssemblyVolume* CloverDetector::ConstructCrystalAssembly(G4String name){

	logicShape2 =   new G4LogicalVolume(GeCrystal7, Ge_mat, name);				// build the logical volume of the detector
	logicShape3 =   new G4LogicalVolume(InsideCylinder1, Ge_mat, name + "_DeadLayer");	// build the inside dead layer
	logicShape1 =   new G4LogicalVolume(DeadGeCrystal8, Ge_mat, name + "_DeadLayer");	// build the inside dead layer
	logicShape4 =   new G4LogicalVolume(BackDeadLayer5, Ge_mat, name + "_DeadLayer");	// build the back dead layer


	G4VisAttributes* vis = new G4VisAttributes(false);
//	logicShape2->SetVisAttributes(vis);
	logicShape3->SetVisAttributes(vis);
	logicShape1->SetVisAttributes(vis);
	logicShape4->SetVisAttributes(vis);


	G4AssemblyVolume* assemblyDetector1 = new G4AssemblyVolume();
		

	G4RotationMatrix *  yRot2 = new G4RotationMatrix;
	yRot2->rotateY(0*deg);

	
	G4ThreeVector pos1 = G4ThreeVector((CrystalInsideLength+err)+SpaceBetweenCrystals/2,(CrystalInsideLength+err)+SpaceBetweenCrystals/2,CrystalLength/2);
	G4ThreeVector pos2 = G4ThreeVector((CrystalInsideLength+err+SpaceBetweenCrystals/2),(CrystalInsideLength+err+SpaceBetweenCrystals/2),
				CrystalLength/2-BackDeadLayer/2+2*err+CrystalLength/2);	
	G4ThreeVector pos3 = G4ThreeVector((CrystalInsideLength+err+SpaceBetweenCrystals/2),(CrystalInsideLength+err+SpaceBetweenCrystals/2),	
				-((WellLength-BackDeadLayer))/2+CrystalLength/2-BackDeadLayer+CrystalLength/2);


	assemblyDetector1->AddPlacedVolume(logicShape2,pos1,yRot2);
	assemblyDetector1->AddPlacedVolume(logicShape1,pos1,yRot2);	
	assemblyDetector1->AddPlacedVolume(logicShape4,pos2,yRot2);	
	assemblyDetector1->AddPlacedVolume(logicShape3,pos3,yRot2);	

	AddActiveVolumeName(name);

	AppsInput* FInput = AppsInput::Instance();
	
	vector<AppsAnalysisObject *> temp;		

	AppsAnalysisObject * analysis = new AppsAnalysisObject();
	analysis->SetObjectName(name);
	analysis->SetDetectorResolution(0.000228,0.000731);		
	analysis->SetObjectHistograms(FInput->GeometryObjects[objectId].GetObjectHistogramType());

	temp.push_back(analysis);

	FInput->GeometryObjects[objectId].AnalysisObjects.push_back(temp);


	

	return assemblyDetector1;
}


void CloverDetector::ConstructHolder(G4String name){

	G4VSolid* Holder1 = new G4Trd("GeTrap1",HolderFrontWidth/2,HolderWidth/2,HolderFrontWidth/2,HolderWidth/2,HolderFrontLength/2);
	G4VSolid* Holder2 = new G4Trd("GeTrap2",(HolderFrontWidth-HolderFrontWallThickness)/2,(HolderWidth-HolderFrontWallThickness)/2,
			(HolderFrontWidth-HolderFrontWallThickness)/2,(HolderWidth-	HolderFrontWallThickness)/2,HolderFrontLength/2);
 	G4VSolid* Holder3 = new G4SubtractionSolid ("GeTrap3", Holder1, Holder2,0,G4ThreeVector(0,0,HolderFrontWallThickness));
  	G4VSolid* Holder4 = new G4Box ("GeBox1",HolderWidth/2,HolderWidth/2,(HolderLength-HolderFrontLength)/2);
	G4VSolid* Holder5 = new G4Box ("GeBox2",(HolderWidth-HolderWallThickness)/2,(HolderWidth-HolderWallThickness)/2,HolderLength);
	G4VSolid* Holder6 = new G4SubtractionSolid ("GeBox3", Holder4, Holder5,0,G4ThreeVector(0,0,HolderFrontWallThickness));
  	G4VSolid* Holder7 = new G4UnionSolid("GeBox3",Holder3,Holder6,0,G4ThreeVector(0,0,(HolderLength)/2));
	G4VSolid* Holder8 = new G4Box ("GeBox1",HolderWidth/2,HolderWidth/2,HolderWallThickness/2);
        G4VSolid* Holder9 = new G4UnionSolid("GeBox3",Holder7,Holder8,0,G4ThreeVector(0,0,HolderLength-HolderFrontLength/2));
	

  	logicShape6 =   new G4LogicalVolume(Holder9, Al_mat, name);	
 
}

void CloverDetector::ConstructShield(){
	
	//Sideshield and its holder (holder first)	
	G4VSolid* SideShielda0 = new G4ExtrudedSolid("SideShielda0",{{0,0},{0,(sideShieldbL-sideShieldaL)/2.},{sideShieldaH,0}},sideShieldbL/2.,{0,0},1.,{0,0},1.);
	
	G4VSolid* SideShieldb0 = new G4Box("SideShieldb0", sideShieldbH/2., (sideShieldbL-sideShieldaL)/4., sideShieldbL/2.);

	G4VSolid* SideShieldb0b = new G4Box("SideShieldb0b", sideShieldbH/2.*1.1, (sideShieldbL-sideShieldaL)/4., sideShieldbL/2.);

	G4RotationMatrix *  yRot1 = new G4RotationMatrix;
	yRot1->rotateX(45*deg);

	G4VSolid* SideShielda00 = new G4SubtractionSolid("SideShielda00", SideShielda0, SideShieldb0, yRot1, G4ThreeVector(0., 0., sideShieldbL/2.+(pow(2.,0.5)-2.)*(sideShieldbL-sideShieldaL)/4.));

	G4VSolid* SideShieldb00 = new G4SubtractionSolid("SideShieldb00", SideShieldb0, SideShieldb0b, yRot1, G4ThreeVector(0., 0., sideShieldbL/2.+(pow(2.,0.5)-1.)*(sideShieldbL-sideShieldaL)/4.));
	
	yRot1->rotateX(90*deg);
	
	G4VSolid* SideShielda000 = new G4SubtractionSolid("SideShielda000", SideShielda00, SideShieldb0, yRot1, G4ThreeVector(0., 0., -sideShieldbL/2.-(pow(2.,0.5)-2.)*(sideShieldbL-sideShieldaL)/4.));
	
	G4VSolid* SideShieldb000 = new G4SubtractionSolid("SideShieldb000", SideShieldb00, SideShieldb0b, yRot1, G4ThreeVector(0., 0., -sideShieldbL/2.-(pow(2.,0.5)-1.)*(sideShieldbL-sideShieldaL)/4.));
	
	G4RotationMatrix *  yRot2 = new G4RotationMatrix;

	G4VSolid* SideShield0 = new G4UnionSolid("SideShield0", SideShielda000, SideShieldb000 , yRot2, G4ThreeVector(-sideShieldbH/2.,(sideShieldbL-sideShieldaL)/4.,0));
	
	//------------
	G4VSolid* SideShieldaa0 = new G4ExtrudedSolid("SideShieldaa0",{{0,0},{0,(sideShieldabL-sideShieldaaL)/2.},{sideShieldaaH,0}},sideShieldabL/2.,{0,0},1.,{0,0},1.);
	
	G4VSolid* SideShieldab0 = new G4Box("SideShieldab0", sideShieldabH/2., (sideShieldabL-sideShieldaaL)/4., sideShieldabL/2.);

	G4VSolid* SideShieldab0b = new G4Box("SideShieldab0b", sideShieldabH/2.*1.1, (sideShieldabL-sideShieldaaL)/4., sideShieldabL/2.);

	yRot1->rotateX(-90*deg);

	G4VSolid* SideShieldaa00 = new G4SubtractionSolid("SideShieldaa00", SideShieldaa0, SideShieldab0, yRot1, G4ThreeVector(sideShieldabH/2.-0.1*mm, 0., sideShieldabL/2.+(pow(2.,0.5)-2.)*(sideShieldabL-sideShieldaaL)/4.));

	G4VSolid* SideShieldab00 = new G4SubtractionSolid("SideShieldab00", SideShieldab0, SideShieldab0b, yRot1, G4ThreeVector(0., 0., sideShieldabL/2.+(pow(2.,0.5)-1.)*(sideShieldabL-sideShieldaaL)/4.));
	
	yRot1->rotateX(90*deg);
	
	G4VSolid* SideShieldaa000 = new G4SubtractionSolid("SideShieldaa000", SideShieldaa00, SideShieldab0, yRot1, G4ThreeVector(sideShieldabH/2.-0.1*mm, 0., -sideShieldabL/2.-(pow(2.,0.5)-2.)*(sideShieldabL-sideShieldaaL)/4.));
	
	G4VSolid* SideShieldab000 = new G4SubtractionSolid("SideShieldab000", SideShieldab00, SideShieldab0b, yRot1, G4ThreeVector(0., 0., -sideShieldabL/2.-(pow(2.,0.5)-1.)*(sideShieldabL-sideShieldaaL)/4.));
	
	SideShield = new G4UnionSolid("SideShield", SideShieldaa000, SideShieldab000 , yRot2, G4ThreeVector(-sideShieldabH/2.,(sideShieldabL-sideShieldaaL)/4.,0));
		
	SideShieldb = new G4SubtractionSolid("SideShieldb", SideShield0, SideShield , yRot2, G4ThreeVector(sideShieldaH-sideShieldaaH-2.*mm, 0.5*mm, 0.));
	
	
	//BackShield and its holder (holder first)
	G4VSolid* BackShield0 = new G4Box("BackShield0", backShieldH/2., backShieldL/2., backShieldL/4.);
	
	G4VSolid* BackShieldT1 = new G4Tubs("BackShieldT1", 0., backShieldR1, backShieldH/2., 0., 360.*deg);
	
	G4VSolid* BackShieldT2 = new G4Tubs("BackShieldT2", 0., backShieldR2, backShieldRH/2., 0., 360.*deg);
	
	G4RotationMatrix *  yRot3 = new G4RotationMatrix;
	yRot3->rotateY(90*deg);

	G4VSolid* BackShield00 = new G4SubtractionSolid("BackShield00", BackShield0, BackShieldT1, yRot3, G4ThreeVector(0., 0., backShieldL/4.));
	
	G4VSolid* BackShield000 = new G4SubtractionSolid("BackShield000", BackShield00, BackShieldT2, yRot3, G4ThreeVector((backShieldRH-backShieldH)/2, 0., backShieldL/4.));

	G4VSolid* BackShielda0 = new G4Box("BackShielda0", backShieldaH/2., backShieldL/2.-backShieldalT, backShieldL/4.-backShieldalT);
	
	G4VSolid* BackShieldaT = new G4Tubs("BackShieldaT", 0., backShieldR1+backShieldalT, backShieldaH/2., 0., 360.*deg);
		
	BackShield = new G4SubtractionSolid("BackShield", BackShielda0, BackShieldaT, yRot3, G4ThreeVector(0., 0., backShieldL/4.));
	
	BackShieldb = new G4SubtractionSolid("BackShieldb", BackShield000, BackShield, yRot2, G4ThreeVector((backShieldH-backShieldaH)/2.-backShieldalT, 0., 0.));
		
	
	//Frontshield and its holder
	frontShieldal = frontShieldcl+2*frontShieldcH*sin(22.5*deg)+(frontShieldaT-frontShieldcT)*cos(22.5*deg);
	frontShieldcL = frontShieldcl+2*frontShieldcH*sin(22.5*deg);
	frontShieldaL = frontShieldal+2*frontShieldaH*sin(22.5*deg);
	FrontShield = new G4Trap("FrontShield", frontShieldaH/2., 0., 0., frontShieldaT/2., frontShieldal/2., frontShieldal/2.-frontShieldaT*sin(67.5*deg), 0., frontShieldaT/2., frontShieldaL/2., frontShieldaL/2.-frontShieldaT*sin(67.5*deg), 0.);
	
	G4VSolid* FrontShieldb0 = new G4Trd("FrontShieldb0", frontShieldbH/2., frontShieldbH/2., frontShieldbL/2., frontShieldbL/2.-frontShieldbT, frontShieldbT/2.);
	
	G4RotationMatrix *  yRot4 = new G4RotationMatrix;
	yRot4->rotateX(-90*deg);
	yRot4->rotateY(90*deg);
	yRot4->rotateX(-22.5*deg);

	FrontShieldb = new G4SubtractionSolid("FrontShieldb", FrontShieldb0, FrontShield, yRot4, G4ThreeVector(frontShieldH-(frontShieldaH+frontShieldcH*2)/2*cos(22.5*deg)+(frontShieldaT-frontShieldcT)/2*sin(22.5*deg)-frontShieldbH/2., 0., -(frontShieldaL+frontShieldal)/4.+frontShieldaT/2.-50*err+frontShieldbL/2.-frontShieldbT/2.));
	
	FrontShieldc = new G4Trap("FrontShieldc", frontShieldcH/2., 0., 0., frontShieldcT/2., frontShieldcl/2., frontShieldcl/2.-frontShieldcT*sin(67.5*deg), 0., frontShieldcT/2., frontShieldcL/2., frontShieldcL/2.-frontShieldcT*sin(67.5*deg), 0.);


	logicShape7 =   new G4LogicalVolume(FrontShieldb, Densimet_mat, "Shield");	

}


void CloverDetector::ConstructBackHolder(){
  
	G4RotationMatrix * zRot1 = new G4RotationMatrix;
 	zRot1->rotateZ(45*deg);

	G4RotationMatrix * zRot2 = new G4RotationMatrix;
 	zRot2->rotateZ(-45*deg);
	
	G4VSolid* HolderBack1 = new G4Box ("HolderBack1",HolderBackWidth/2,HolderBackWidth/2,HolderBackThickness/2);
	G4VSolid* HolderBack2 = new G4Box ("HolderBack2",HolderBackPieceALength/2,HolderBackPieceAWidth/2,HolderBackPieceAThickness/2);
  	G4VSolid* HolderBack3 = new G4UnionSolid("HolderBack3",HolderBack1,HolderBack2,0,G4ThreeVector(HolderBackWidth/2-HolderBackPieceALength/2,0,
				HolderBackThickness/2+HolderBackPieceAThickness/2));
  	G4VSolid* HolderBack4 = new G4UnionSolid("HolderBack4",HolderBack3,HolderBack2,0,G4ThreeVector(-HolderBackWidth/2+HolderBackPieceALength/2,0,
				HolderBackThickness/2+HolderBackPieceAThickness/2));
	G4VSolid* HolderBack5 = new G4Box ("HolderBack5",HolderBackPieceAWidth/2,HolderBackPieceALength/2,HolderBackPieceAThickness/2);  	
	G4VSolid* HolderBack6 = new G4UnionSolid("HolderBack6",HolderBack4,HolderBack5,0,G4ThreeVector(0,-HolderBackWidth/2+HolderBackPieceALength/2,
				HolderBackThickness/2+HolderBackPieceAThickness/2));
	G4VSolid* HolderBack7 = new G4UnionSolid("HolderBack7",HolderBack6,HolderBack5,0,G4ThreeVector(0,+HolderBackWidth/2-HolderBackPieceALength/2,
				HolderBackThickness/2+HolderBackPieceAThickness/2));
	G4VSolid* HolderBack8 = new G4Box ("HolderBack8",(HolderBackWidth-2*HolderBackPieceALength)/2,HolderBackPieceBWidth/2,HolderBackPieceAThickness/2);  	
	G4VSolid* HolderBack9 = new G4UnionSolid("HolderBack9",HolderBack7,HolderBack8,0,G4ThreeVector(0,0,HolderBackPieceAThickness/2));
	G4VSolid* HolderBack10 = new G4Box ("HolderBack10",HolderBackPieceBWidth/2,(HolderBackWidth-2*HolderBackPieceALength)/2,HolderBackPieceAThickness/2);  
	G4VSolid* HolderBack11 = new G4UnionSolid("HolderBack11",HolderBack9,HolderBack10,0,G4ThreeVector(0,0,HolderBackPieceAThickness/2));
	G4VSolid* HolderBack12 = new G4Tubs("HolderBack12",0,HolderBackPieceCRadius,HolderBackPieceCThickness/2, 0.*deg, 360.*deg);
	G4VSolid* HolderBack13 = new G4UnionSolid("HolderBack13",HolderBack11,HolderBack12,0,G4ThreeVector(0,0,HolderBackPieceCThickness/2));
	G4VSolid* HolderBack14 = new G4Box ("HolderBack14",HolderBackWidth,HolderBackWidth/3.5,HolderBackWidth/2);
	G4VSolid* HolderBack15 = new G4Box ("HolderBack15",HolderBackWidth/2,HolderBackWidth/2,HolderBackWidth/4);
 	G4VSolid* HolderBack16 = new G4SubtractionSolid ("HolderBack16", HolderBack15, HolderBack14,zRot1,G4ThreeVector(0,0,1*cm));
 	G4VSolid* HolderBack17 = new G4SubtractionSolid ("HolderBack17", HolderBack16, HolderBack14,zRot2,G4ThreeVector(0,0,1*cm));
	G4VSolid* HolderBack18 = new G4UnionSolid("HolderBack18",HolderBack13,HolderBack17,0,G4ThreeVector(0,0,-HolderBackWidth/4-HolderBackThickness/2));

  	logicShape5 =   new G4LogicalVolume(HolderBack18, Al_mat, "HolderBack");	
	G4VisAttributes* vis = new G4VisAttributes(false);
	logicShape5->SetVisAttributes(vis);
}


void CloverDetector::ConstructCrystalSegments(){


	G4VSolid* SegBox1 = new G4Box ("GeBox1",length,length,length*2);

	double dev = -1.6*mm;
	double dev2 = 0.6*mm;



	G4VSolid* FrontSegments = new G4SubtractionSolid ("FrontSegments", GeCrystal7, SegBox1,0,G4ThreeVector(0,0,length*2-(CrystalLength/2-CrystalFrontLength)));
	G4VSolid* BackSegments = new G4SubtractionSolid ("BackSegments", GeCrystal7, SegBox1,0,G4ThreeVector(0,0,(-CrystalLength/2.-length*2+CrystalFrontLength)));

	G4VSolid* Seg11 = new G4SubtractionSolid ("Seg11", FrontSegments, SegBox1,0,G4ThreeVector(length+dev,0,0));
	Seg1 = new G4SubtractionSolid ("Seg1", Seg11, SegBox1,0,G4ThreeVector(0,length+dev,0));

	G4VSolid* Seg21 = new G4SubtractionSolid ("Seg21", FrontSegments, SegBox1,0,G4ThreeVector(length+dev,0,0));
	Seg2 = new G4SubtractionSolid ("Seg2", Seg21, SegBox1,0,G4ThreeVector(0,-length+dev,0));

	G4VSolid* Seg31 = new G4SubtractionSolid ("Seg31", FrontSegments, SegBox1,0,G4ThreeVector(-length+dev,0,0));
	Seg3 = new G4SubtractionSolid ("Seg3", Seg31, SegBox1,0,G4ThreeVector(0,-length+err+dev,0));
	
	G4VSolid* Seg41 = new G4SubtractionSolid ("Seg41", FrontSegments, SegBox1,0,G4ThreeVector(-length+dev,0,0));
	Seg4 = new G4SubtractionSolid ("Seg4", Seg41, SegBox1,0,G4ThreeVector(0,length+dev,0));

	G4VSolid* Seg51 = new G4SubtractionSolid ("Seg51", BackSegments, SegBox1,0,G4ThreeVector(0,length+dev2,0));
	Seg5 = new G4SubtractionSolid ("Seg5", Seg51, SegBox1,0,G4ThreeVector(length+dev2,0,0));

	G4VSolid* Seg61 = new G4SubtractionSolid ("Seg61", BackSegments, SegBox1,0,G4ThreeVector(length+dev2,0,0));
	Seg6 = new G4SubtractionSolid ("Seg6", Seg61, SegBox1,0,G4ThreeVector(0,-length+dev2,0));

	G4VSolid* Seg71 = new G4SubtractionSolid ("Seg71", BackSegments, SegBox1,0,G4ThreeVector(-length+dev2,0,0));
	Seg7 = new G4SubtractionSolid ("Seg7", Seg71, SegBox1,0,G4ThreeVector(0,-length+dev2,0));
	
	G4VSolid* Seg81 = new G4SubtractionSolid ("Seg81", BackSegments, SegBox1,0,G4ThreeVector(-length+dev2,0,0));
	Seg8 = new G4SubtractionSolid ("Seg8", Seg81, SegBox1,0,G4ThreeVector(0,length+dev2,0));

}



void CloverDetector::ConstructInsideDeadLayer(){

	InsideCylinder1 = new G4Tubs ("InsideCylinder1",WellRadius,WellRadius+WellLateralDeadLayer,(WellLength-BackDeadLayer)/2,0.*deg,360.*deg);

}

void CloverDetector::ConstructDetectorCrystal(){
	
	G4RotationMatrix*  zRot1 = new G4RotationMatrix;
  	zRot1->rotateZ(90*deg);	

	G4VSolid* InsideCylinder5 = new G4Tubs ("InsideCylinder5",0,WellRadius+WellLateralDeadLayer,(WellLength)/2,0.*deg,360.*deg);
  	GeCrystal7 = new G4SubtractionSolid ("GeCrystal7",GeCrystal6,InsideCylinder5,zRot1,G4ThreeVector(0,0,-((WellLength))/2+CrystalLength/2)); 

}


void CloverDetector::ConstructDeadLayer(){



	G4RotationMatrix* xRot1 = new G4RotationMatrix;
	xRot1->rotateX(CrystalCutAngle);

	G4RotationMatrix* yRot1 = new G4RotationMatrix;
  	yRot1->rotateY(-CrystalCutAngle);   


	G4VSolid*  GeBox1 = new G4Box ("GeBox1",length,length,length*2);
  	G4VSolid* DeadGeCylinder1 = new G4Tubs ("DeadGeCylinder1",0*mm, CrystalRadius, CrystalLength/2, 0.*deg, 360.*deg);
  	G4VSolid* DeadGeCrystal1 = new G4SubtractionSolid ("DeadGeCrystal1", DeadGeCylinder1, GeBox1,0,G4ThreeVector(0,length+CrystalOutsideLength,0));
  	G4VSolid* DeadGeCrystal2 = new G4SubtractionSolid ("DeadGeCrystal2", DeadGeCrystal1, GeBox1,0,G4ThreeVector(length+CrystalOutsideLength,0,0));
  	G4VSolid* DeadGeCrystal3 = new G4SubtractionSolid ("DeadGeCrystal3", DeadGeCrystal2, GeBox1,0,G4ThreeVector(-(length+CrystalInsideLength),0,0));
  	G4VSolid* DeadGeCrystal4 = new G4SubtractionSolid ("DeadGeCrystal4", DeadGeCrystal3, GeBox1,0,G4ThreeVector(0,-(length+CrystalInsideLength),0));
  	G4VSolid* DeadGeCrystal5 = new G4SubtractionSolid ("DeadGeCrystal5", DeadGeCrystal4, GeBox1,xRot1,G4ThreeVector(0,CrystalOutsideLength+length,
					-(CrystalLength/2-CrystalCutDistance)));
  	G4VSolid* DeadGeCrystal6 = new G4SubtractionSolid ("DeadGeCrystal6", DeadGeCrystal5, GeBox1,yRot1,G4ThreeVector(CrystalOutsideLength+length,0,
					-(CrystalLength/2-CrystalCutDistance)));
  	G4VSolid* DeadGeCrystal7 = new G4SubtractionSolid ("DeadGeCrystal7", DeadGeCrystal6, GeCrystal6,0,G4ThreeVector(0,0,BackDeadLayer+err));  
  	DeadGeCrystal8 = new G4SubtractionSolid ("DeadGeCrystal8", DeadGeCrystal7, GeCrystal6,0,G4ThreeVector(0,0,0));    

}



void CloverDetector::ConstructBackDeadLayer(){

	G4VSolid*  GeBox1 = new G4Box ("GeBox1",length,length,length*2);
  	G4VSolid* BackDeadLayer1 = new G4Tubs ("BackDeadLayer1",WellRadius+err,CrystalRadius-LateralDeadLayer-err,BackDeadLayer/2, 0.*deg, 360.*deg);
  	G4VSolid* BackDeadLayer2 = new G4SubtractionSolid ("GeCrystal1", BackDeadLayer1, GeBox1,0,G4ThreeVector(0,length+CrystalOutsideLength-LateralDeadLayer,0));
  	G4VSolid* BackDeadLayer3 = new G4SubtractionSolid ("GeCrystal2", BackDeadLayer2, GeBox1,0,G4ThreeVector(length+CrystalOutsideLength-LateralDeadLayer,0,0));
  	G4VSolid* BackDeadLayer4 = new G4SubtractionSolid ("GeCrystal3", BackDeadLayer3, GeBox1,0,G4ThreeVector(-(length+CrystalInsideLength-LateralDeadLayer),0,0));
  	BackDeadLayer5 = new G4SubtractionSolid ("GeCrystal4", BackDeadLayer4, GeBox1,0,G4ThreeVector(0,-(length+CrystalInsideLength-LateralDeadLayer),0));
  
}


void CloverDetector::ConstructCrystal(){ 


	G4RotationMatrix* xRot1 = new G4RotationMatrix;
	xRot1->rotateX(CrystalCutAngle);

	G4RotationMatrix* yRot1 = new G4RotationMatrix;
  	yRot1->rotateY(-CrystalCutAngle);   

	G4VSolid*  GeCylinder1 = new G4Tubs ("GeCylinder1",0*mm, CrystalRadius-LateralDeadLayer, CrystalLength/2-FrontDeadLayer/2-BackDeadLayer/2, 0.*deg, 360.*deg);
	G4VSolid*  GeBox1 = new G4Box ("GeBox1",length,length,length*2);
	G4VSolid*  GeCrystal1 = new G4SubtractionSolid ("GeCrystal1", GeCylinder1, GeBox1,0,G4ThreeVector(0,length+CrystalOutsideLength-LateralDeadLayer,0));
	G4VSolid*  GeCrystal2 = new G4SubtractionSolid ("GeCrystal2", GeCrystal1, GeBox1,0,G4ThreeVector(length+CrystalOutsideLength-LateralDeadLayer,0,0));
	G4VSolid*  GeCrystal3 = new G4SubtractionSolid ("GeCrystal3", GeCrystal2, GeBox1,0,G4ThreeVector(-(length+CrystalInsideLength-LateralDeadLayer),0,0));
	G4VSolid*  GeCrystal4 = new G4SubtractionSolid ("GeCrystal4", GeCrystal3, GeBox1,0,G4ThreeVector(0,-(length+CrystalInsideLength-LateralDeadLayer),0));
	G4VSolid*  GeCrystal5 = new G4SubtractionSolid ("GeCrystal5", GeCrystal4, GeBox1,xRot1,G4ThreeVector(0,CrystalOutsideLength+length-LateralDeadLayer,
					-(CrystalLength/2-CrystalCutDistance)));
	GeCrystal6 = new G4SubtractionSolid ("GeCrystal6", GeCrystal5, GeBox1,yRot1,G4ThreeVector(CrystalOutsideLength+length-LateralDeadLayer,0,
					-(CrystalLength/2-CrystalCutDistance)));


}

void CloverDetector::ConstructMaterials(){

	G4NistManager* nist = G4NistManager::Instance();

	Ge_mat = nist->FindOrBuildMaterial("G4_Ge");
	Al_mat = nist->FindOrBuildMaterial("G4_Al");
	BGO_mat = nist->FindOrBuildMaterial("G4_BGO");
	CsI_mat = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
	Pb_mat = nist->FindOrBuildMaterial("G4_Pb");
	Air_mat = nist->FindOrBuildMaterial("G4_AIR");

	  // StainlessSteel definition
	G4double density;
   	G4int ncomponents;
   	G4double fractionmass;
  	G4Element* C  = nist->FindOrBuildElement("C");
   	G4Element* Si = nist->FindOrBuildElement("Si");
   	G4Element* Cr = nist->FindOrBuildElement("Cr");
   	G4Element* Mn = nist->FindOrBuildElement("Mn");
   	G4Element* Fe = nist->FindOrBuildElement("Fe");
   	G4Element* Ni = nist->FindOrBuildElement("Ni");
   	StainlessS_mat = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
   	StainlessS_mat->AddElement(C, fractionmass=0.001);
   	StainlessS_mat->AddElement(Si, fractionmass=0.007);
   	StainlessS_mat->AddElement(Cr, fractionmass=0.18);
   	StainlessS_mat->AddElement(Mn, fractionmass=0.01);
   	StainlessS_mat->AddElement(Fe, fractionmass=0.712);
   	StainlessS_mat->AddElement(Ni, fractionmass=0.09);



 	 // Densimet definition
   	G4int mcomponents;
  	G4Element* W  = nist->FindOrBuildElement("W");
   	Densimet_mat = new G4Material("Densimet", density= 18.55*g/cm3, mcomponents=3);
   	Densimet_mat->AddElement(W, fractionmass=0.97);   
   	Densimet_mat->AddElement(Fe, fractionmass=0.015);
   	Densimet_mat->AddElement(Ni, fractionmass=0.015);



}

void CloverDetector::ReadInputParameters(){


	G4String z; 
	G4double a[40];

	ifstream myfile ("CloverParameters.txt");
	if(myfile.is_open()){
		getline(myfile, z);
		getline(myfile, z);
   		getline(myfile, z);
   
		for(G4int i=0;i<40;i++){
   
			myfile >> a[i];
  	 		getline(myfile, z);
  			getline(myfile, z); 

		}
	} else {
   		cout << "Unable to open CloverParameters.txt" << G4endl;
	}  

	CrystalRadius = a[0]*cm; 
	CrystalLength = a[1]*cm;  
	CrystalInsideLength=a[2]*cm; 
	CrystalOutsideLength=a[3]*cm; 
	FrontDeadLayer=a[4]*cm;  
	BackDeadLayer=a[5]*cm; 
	LateralDeadLayer=a[6]*cm; 
	WellRoundRadius=a[7]*cm;
	WellRadius=a[8]*cm;  	
	WellBottomDeadLayer=a[9]*cm;
	WellLateralDeadLayer=a[10]*cm;  
	WellLength=a[11]*cm;
	HolderLength=a[12]*cm; 
	HolderWidth=a[13]*cm; 
	HolderRoundRadius=a[14]*cm;  
	HolderFrontRoundRadius=a[15]*cm; 
	HolderFrontLength=a[16]*cm;  
	HolderFrontWidth=a[17]*cm; 
	HolderWallThickness=a[18]*cm;  
	HolderFrontWallThickness= a[19]*cm;
	SpaceBetweenCrystals=a[20]*cm; 
	HolderBackWidth=a[21]*cm; 
	HolderBackThickness=a[22]*cm; 
	HolderBackPieceALength=a[23]*cm;  
	HolderBackPieceAWidth=a[24]*cm; 
	HolderBackPieceAThickness= a[25]*cm; 
	HolderBackPieceBWidth=a[26]*cm;
	HolderBackPieceCRadius=a[27]*cm;
	HolderBackPieceCThickness= a[28]*cm; 
	HolderBackHoleRadius=a[29]*cm; 
	HolderBackHolePosition=a[30]*cm; 
	HolderBackPieceLength2=a[31]*cm; 
	InsideAluminiumCylinderLength=a[32]*cm; 
	CrystalCutAngle= a[33]*deg;
	HolderAngle=a[34]*deg; 
	HolderToCrystalDistance=a[35]*cm;
	CrystalFrontLength= a[36]*cm; 
	FrontToBackHolderDistance = a[37]*cm;
	CrystalCutDistance=a[38]*cm;
	CrystalCutDeadLayer=a[39]*cm;
	length=5*cm;
	err=0.001*cm;

}

void CloverDetector::SetFWHMParameters(G4double slope, G4double intercept){
	
	FWHMParameters.push_back(slope);
	FWHMParameters.push_back(intercept);

}	


