#ifndef GammaGenerator_hh
#define GammaGenerator_hh 1

#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4LorentzVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4String.hh"

#include <iostream>
#include <fstream>
using namespace std;

class MCGamma {
public:
  MCGamma();
  virtual ~MCGamma();
  
  G4double energy, time;
  G4ThreeVector position, polarization;
  G4ParticleMomentum direction;
};

class GammaGenerator {
public:
  GammaGenerator(CLHEP::HepRandomEngine*);
  virtual ~GammaGenerator();

  MCGamma* Gamma();   // fills the fGamma object with new data

private:
  CLHEP::HepRandomEngine* randGen;
  G4LorentzVector Electron4Mom, Photon4Mom, Gamma4Mom; 
  G4ThreeVector Electron3Pos, Photon3Pol, GammaStokes, Electron3Boost;
  MCGamma* fGamma;
  G4RandFlat*  flatRand;
  G4RandExponential*  expRand;
  G4RandGauss* gaussRand;
  G4int nGamma; G4double XSum; G4double mec2;
  G4double eEner, deEner, eSigZ, eEmit, eAlpha, eBeta,  eNB, eChar;   // electron beam parameters at IP
  G4double lLamb, dlEner, lRayl, lSigT, lPowr, lPolD, lPolA;                     // laser beam parameters at IP
  G4double elCross, elFreq;                                                                        // electron-laser interaction parameters
  G4bool useRecoil, useNonLin;                                                                 // switches for various effects
  G4bool printOut;
  
  void GenerateElectron();   // generate electron in the lab frame
  void NonLinEffect();           // apply nonlinear effects to electron
  void GeneratePhoton();     // generate photon in the lab frame
  void GoToRestFrame();      // move to electron rest frame with aligned photon
  void GenerateGamma();    // generate gamma
  void ComputeIntensity();  // compute and print gamma source intensity
  void Print4Vector(G4String, G4LorentzVector);
  void PrintParameters();
};

#endif
