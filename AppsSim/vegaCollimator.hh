#ifndef VEGACOLLIMATOR_H
#define VEGACOLIMATER_H

#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "AppsGeometryObject.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4MultiUnion.hh"

class vegaCollimator
{

public:
    vegaCollimator(G4int id);
    G4AssemblyVolume *GetAssembly() { return collimatorAssembly; }
    void ConstructCollimator(int HoleOnBeam);

private:
    void ConstructMaterials();

    G4int objectId;

    G4Material *w_mat;
    G4Material *al_mat;
    G4Material *air_mat;
    G4Material *pb_mat;
    G4Material *conc_mat;
    G4Material *pe_mat;

    G4AssemblyVolume *collimatorAssembly;
};

#endif
