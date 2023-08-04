//******************************************************************************
// GEANT4 physics class-­‐-­‐>NRFCrossSectionDataStore -­‐-­‐ header file	
// N Kikuzawa, JAEA, 29-­‐JUN-­‐09	
//  Class Description
// This is the class to which to register data-­‐sets. You can get the instance	
// from energy hadronic process, and use its 'AddDataSet(...)' method to tailor	
// the cross-­‐sections for your application.	
//******************************************************************************** 

#ifndef	 NRFCrossSectionDataStore_h	
#define	 NRFCrossSectionDataStore_h 1	
#include "G4ParticleDefinition.hh"	
#include "G4DynamicParticle.hh"	
#include "G4Element.hh"	 
#include "NRFCrossSectionData.hh"	
  
//class NRFCrossSectionData;

class  NRFCrossSectionDataStore	
 {
  
  public:

	NRFCrossSectionDataStore();
	~NRFCrossSectionDataStore();	
   	
	G4double GetCrossSection(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope);

  	G4double GetBranching(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope, G4double aTemperature=0.);
	G4double GetExcitation(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotopet, G4double aTemperature=0.);

	G4bool IsApplicable(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope);


   	void AddDataSet(NRFCrossSectionData*);		
 	void BuildPhysicsTable(const G4ParticleDefinition&);	
  	void DumpPhysicsTable(const G4ParticleDefinition&);	
   	G4double GetLastSpin(){ return theLastSpin;};	

  	G4int GetLastParity() { return theLastParity;};	
  	G4ThreeVector GetMomentumDirection( const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope);
  
	G4double FirstExcitedLevel;
  	G4double LastBranchingRatio;

  private:
	
  	std::vector<NRFCrossSectionData*> DataSetList;		
  	G4int NDataSetList;	
 	G4int verboseLevel;	
  
 	G4double theLastCrossSection;	
 	G4ThreeVector theLastMomentumDirection;		
  	G4double theLastSpin;
  	G4int theLastParity;	


  };
#endif
//**************************************************************************************


  









  




