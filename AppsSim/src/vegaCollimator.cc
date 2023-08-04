#include "vegaCollimator.hh"
#include "AppsInput.hh"

vegaCollimator::vegaCollimator(G4int id) : objectId(id)
{

    ConstructMaterials();
    AppsInput *FInput = AppsInput::Instance();

    auto object = FInput->GetGeometryObjects()[id];

    ConstructCollimator(object.GetShapeParameters()[0]);
}

void vegaCollimator::ConstructCollimator(int HoleOnBeam)
{

    collimatorAssembly = new G4AssemblyVolume();

    // the optical cavity window:

    G4double Wind_locIP = 5.10 * m, Wind_length = 1 * mm, Wind_diam = 1 * cm;
    G4Tubs *Wind_shape = new G4Tubs("OPWindow", 0, Wind_diam / 2, Wind_length, 0 * deg, 360 * deg);
    G4ThreeVector Wind_cent = G4ThreeVector(0, 0, Wind_locIP + Wind_length / 2);
    G4LogicalVolume *logicWind = new G4LogicalVolume(Wind_shape, al_mat, "WindowVolume");
    collimatorAssembly->AddPlacedVolume(logicWind, Wind_cent, 0);

    // the beam stopper:
    G4double Stop_locIP = 7.68 * m,
             Stop_length = 20 * cm, Stop_height = 15 * cm, Stop_width = 24.4 * cm, Stop_hole = 5 * mm;
    G4Box *Stop_box = new G4Box("StopBox", Stop_width / 2, Stop_height / 2, Stop_length / 2);
    G4Tubs *Stop_cut = new G4Tubs("StopHole", 0, Stop_hole / 2, Stop_length, 0 * deg, 360 * deg);
    G4VSolid *Stop_shape = new G4SubtractionSolid("StopShape", Stop_box, Stop_cut);
    G4ThreeVector Stop_cent = G4ThreeVector(0, 0, Stop_locIP + Stop_length / 2);
    G4LogicalVolume *logicStop = new G4LogicalVolume(Stop_shape, w_mat, "StopperVolume");
    collimatorAssembly->AddPlacedVolume(logicStop, Stop_cent, 0);

    // 1st air gap:
    G4double Air1_locIP = 7.88 * m, Air1_length = 14 * cm;
    G4Box *Air1_shape = new G4Box("Air1Shape", 1 * m, 1 * m, Air1_length / 2);
    G4ThreeVector Air1_cent = G4ThreeVector(0, 0, Air1_locIP + Air1_length / 2);
    G4LogicalVolume *logicAir1 = new G4LogicalVolume(Air1_shape, air_mat, "Air1Volume");
    collimatorAssembly->AddPlacedVolume(logicAir1, Air1_cent, 0);

    // the neutron shield (inside wall):
    G4double NeSh_locIP = 8.02 * m, NeSh_length = 20 * cm, NeSh_diam = 14 * cm, NeSh_hole = 25 * mm;
    G4Tubs *NeSh_shape = new G4Tubs("NeShTube", NeSh_hole / 2, NeSh_diam / 2, NeSh_length / 2, 0. * deg, 360. * deg);
    G4ThreeVector NeSh_cent = G4ThreeVector(0, 0, NeSh_locIP + NeSh_length / 2);
    G4LogicalVolume *logicNeSh = new G4LogicalVolume(NeSh_shape, pe_mat, "NeShVolume");
    collimatorAssembly->AddPlacedVolume(logicNeSh, NeSh_cent, 0);

    // 2nd air gap (inside wall):
    G4double Air2_locIP = 8.22 * m, Air2_length = 70 * cm, Air2_diam = 14 * cm;
    G4Tubs *Air2_shape = new G4Tubs("Air2Tube", 0, Air2_diam / 2, Air2_length / 2, 0. * deg, 360. * deg);
    G4ThreeVector Air2_cent = G4ThreeVector(0, 0, Air2_locIP + Air2_length / 2);
    G4LogicalVolume *logicAir2 = new G4LogicalVolume(Air2_shape, air_mat, "Air2Volume");
    collimatorAssembly->AddPlacedVolume(logicAir2, Air2_cent, 0);

    // the pre-collimator (inside wall):
    G4double PreC_locIP = 8.92 * m, PreC_length = 10 * cm, PreC_diam = 14 * cm, PreC_hole = 10 * mm;
    G4Tubs *PreC_shape = new G4Tubs("PreCTube", PreC_hole / 2, PreC_diam / 2, PreC_length / 2, 0. * deg, 360. * deg);
    G4ThreeVector PreC_cent = G4ThreeVector(0, 0, PreC_locIP + PreC_length / 2);
    G4LogicalVolume *logicPreC = new G4LogicalVolume(PreC_shape, pb_mat, "PreCVolume");
    collimatorAssembly->AddPlacedVolume(logicPreC, PreC_cent, 0);

    // the wall:
    G4double Wall_locIP = NeSh_locIP;
    G4double Wall_length = NeSh_length + Air2_length + PreC_length;
    G4Box *Wall_box = new G4Box("WallBox", 2 * m, 2 * m, Wall_length / 2);
    if (abs(PreC_diam - Air2_diam) / mm > 1.e-03 * mm || abs(PreC_diam - NeSh_diam) / mm > 1.e-03 * mm)
    {
        G4cout << "CollimationSystem: non-uniform wall penetrations. Wall geometry needs to be rewritten!" << G4endl;
        exit(0);
    }
    G4Tubs *Cut_cyl = new G4Tubs("CutCyl", 0., PreC_diam / 2, Wall_length / 2, 0. * deg, 360. * deg);
    G4VSolid *Wall_shape = new G4SubtractionSolid("WallShape", Wall_box, Cut_cyl);
    G4LogicalVolume *logicWall = new G4LogicalVolume(Wall_shape, conc_mat, "WallVolume");
    G4ThreeVector Wall_cent = G4ThreeVector(0, 0, Wall_locIP + Wall_length / 2);
    collimatorAssembly->AddPlacedVolume(logicWall, Wall_cent, 0);

    // 3rd air gap:
    G4double Air3_locIP = 9.02 * m, Air3_length = 7.9 * cm; // bit shorter to make room for a test layer
    G4Box *Air3_shape = new G4Box("Air3Shape", 1 * m, 1 * m, Air3_length / 2);
    G4ThreeVector Air3_cent = G4ThreeVector(0, 0, Air3_locIP + Air3_length / 2);
    G4LogicalVolume *logicAir3 = new G4LogicalVolume(Air3_shape, air_mat, "Air3Volume");
    collimatorAssembly->AddPlacedVolume(logicAir3, Air3_cent, 0);

    // THE GAMMA BEAM COLLIMATOR:
    G4double Coll_locIP = 9.10 * m, Plate_width = 11 * cm, Plate_height = 2 * cm, Plate_thick = 1 * cm, Plate_space = 1 * cm;
    const G4int nPlates = 20;
    const G4int nHoles = 7;
    G4double Hole_cent[nHoles] = {20 * mm, 30 * mm, 40 * mm, 50 * mm, 60 * mm, 70 * mm, 85 * mm}; // hole center starting from left
    G4double Hole_diam[nHoles] = {0.5 * mm, 0.7 * mm, 0.9 * mm, 1.1 * mm, 1.4 * mm, 3 * mm, 10 * mm};
    HoleOnBeam -= 1; // turn it into index
    if (HoleOnBeam < 0 || HoleOnBeam >= nHoles)
    {
        G4cout << "CollimationSystem Error - Wrong hole:" << HoleOnBeam << G4endl;
        exit(-1);
    }
    else
    {
        G4cout << "DetectorConstruction: collimator hole diameter = " << Hole_diam[HoleOnBeam] / mm << "mm" << G4endl;
    }
    G4double alignHole = Plate_width / 2 - Hole_cent[HoleOnBeam];
    G4Box *Plate_box = new G4Box("PlateBox", Plate_width / 2, Plate_height / 2, Plate_thick / 2);
    G4Tubs *Plate_hole[nHoles];
    G4String ObjName;
    G4MultiUnion *HolesUnion = new G4MultiUnion("Holes_Union"); // collect all holes to be cut out from the plate
    for (G4int ih = 0; ih < nHoles; ih++)
    {
        ObjName = "Hole";
        ObjName += std::to_string(ih);
        Plate_hole[ih] = new G4Tubs(ObjName, 0, Hole_diam[ih] / 2, Plate_thick / 2, 0. * deg, 360. * deg);
        G4Transform3D HolePlace = G4Transform3D(G4RotationMatrix(), G4ThreeVector(-Plate_width / 2 + Hole_cent[ih], 0, 0));
        HolesUnion->AddNode(*Plate_hole[ih], HolePlace);
    }
    HolesUnion->Voxelize();
    G4VSolid *Plate_shape = new G4SubtractionSolid("PlateShape", Plate_box, HolesUnion);
    G4Box *Space_box = new G4Box("SpaceBox", Plate_width / 2, Plate_height / 2, Plate_thick / 2);
    G4VSolid *Plate[nPlates];
    G4VSolid *Space[nPlates];
    for (G4int ip = 0; ip < nPlates; ip++)
    {
        Plate[ip] = (G4VSolid *)Plate_shape->Clone();
        ObjName = "Plate";
        ObjName += std::to_string(ip);
        G4LogicalVolume *logicPlate = new G4LogicalVolume(Plate[ip], w_mat, ObjName);
        G4double Plate_cent = Coll_locIP + ip * (Plate_thick + Plate_space) + Plate_thick / 2;
        ObjName = "CollPlate";
        ObjName += std::to_string(ip);

        G4ThreeVector platePos = G4ThreeVector(alignHole, 0, Plate_cent);
        collimatorAssembly->AddPlacedVolume(logicPlate, platePos, 0);

        Space[ip] = (G4VSolid *)Space_box->Clone();
        ObjName = "Space";
        ObjName += std::to_string(ip);
        G4LogicalVolume *logicSpace = new G4LogicalVolume(Space[ip], air_mat, ObjName);
        G4double Space_cent = Coll_locIP + ip * (Plate_thick + Plate_space) + Plate_thick + Plate_space / 2;
        ObjName = "CollSpace";
        ObjName += std::to_string(ip);

        G4ThreeVector spacePos = G4ThreeVector(alignHole, 0, Space_cent);
        collimatorAssembly->AddPlacedVolume(logicSpace, spacePos, 0);
    }
}

void vegaCollimator::ConstructMaterials()
{
    G4NistManager *nist = G4NistManager::Instance();

    w_mat = nist->FindOrBuildMaterial("G4_W");
    al_mat = nist->FindOrBuildMaterial("G4_Al");
    air_mat = nist->FindOrBuildMaterial("G4_AIR");
    pb_mat = nist->FindOrBuildMaterial("G4_Pb");
    conc_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    pe_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE"); // ADD BORON!!!
}
