//*******************************************************************************
//GEANT4 physics class-­‐-­‐> NRFCrossSectionData -­‐-­‐ header file	
// Created: N. Kikuzawa, Japan Atomic Energy Agency, 24-­‐NOV-­‐06	
// The last update: N. Kikuzawa, JAEA,	 24-­‐NOV-­‐06	
//  Modified by Hani Negm to include the angular distribution	
//  for  the  polarization  effect.	 2012 -­2014	
//********************************************************************************

#ifndef	 NRFCrossSectionData_h  
#define NRFCrossSectionData_h 1	
#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"	  
#include "G4Element.hh"	
#include "G4Isotope.hh"
  
class  MyLevelData  {	
  public:	
 	G4double Excite;	
 	G4double   Spin;
 	G4int Parity;	
 	G4double Width;	
 	G4double CrossSection;	
  	G4double F2;	
  	G4double BranchingRatio;
  	G4double ExcitedLevel;
};

class NRFCrossSectionData : public  G4VCrossSectionDataSet	
{

  public:

	NRFCrossSectionData(G4Isotope* anIsotope);
   	~NRFCrossSectionData();	
   	
	G4bool IsApplicable(const G4DynamicParticle* aParticle, const G4Isotope* anIsotope);
	void Init(G4Isotope* anIsotope);

	G4double GetCrossSection(const G4DynamicParticle* aDynamicGamma);

  	G4ThreeVector GetMomentumDirection(const G4DynamicParticle *aParticle);
	G4ThreeVector ComputeMomentumDirection(const G4DynamicParticle* aParticle);	

 	void DumpPhysicsTable(const G4ParticleDefinition&);
	void PrintLevelData();	


	G4double GetZ() {return theZ;}	
  	G4double GetN() {return theN;}
	G4String GetName() {return theName;}
  	G4double GetLastSpin() {return theLastSpin;}
  	G4int GetLastParity() {return theLastParity;}
  	G4double GetLastBranching(){ return theLastBranchingRatio;}
	

	G4double theLastBranchingRatio;
  	G4double theFirstExcitedLevel;
  	typedef std::vector<MyLevelData*> LevelData;

  	LevelData* theLevelData;
	


  private:	
  
	G4String theName;	
  	G4double theN;	
  	G4double theZ;	
  	G4double theLastCrossSection;	
  	G4double theLastSpin;	
  	G4double theGroundStateSpin;	
  	G4int theLastParity;	
 	G4int theGroundStateParity;	
  	G4ThreeVector theLastMomentumDirection;	
  	G4bool hasCrossSectionData;	
  	
	

};

#endif	
  

