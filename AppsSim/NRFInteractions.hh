//***********************************************************************
//GEANT4 Physics->NRFInteractions class: header file
//***********************************************************************

#ifndef  NRFInteractions_h	
  
#define	  NRFInteractions_h 1	
  
#include "G4ios.hh"	
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"	
#include "G4Element.hh"	
#include "G4Gamma.hh"	
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Step.hh"
#include "NRFCrossSectionDataStore.hh"
#include "NRFCrossSectionData.hh"

//----------------------------------------------------------

//class NRFCrossSectionData

class NRFInteractions : public G4VDiscreteProcess	
{
  public: 
	
	NRFInteractions(const G4String&	processName ="NRF", G4ProcessType type = fPhotolepton_hadron);
	~NRFInteractions();

 	G4bool IsApplicable(const G4ParticleDefinition&);	//true for Gamma only.	
  

 	void BuildPhysicsTable(const G4ParticleDefinition&);	// here dummy, the total cross	 section parameterization is used rather	
 						 		// than tables, just calling PrintInfoDefinition	

  
  	void PrintInfoDefinition();				// Print few lines of information about the process: validity range,	
  								// origine ..etc..	
  								//  Invoked  by BuildThePhysicsTable().	
   	
  	void SetCrossSecFactor(G4double fac);			// Set the factor to artificially increase the crossSection (default 1)	
   	
  	G4double GetCrossSecFactor() {return CrossSecFactor;}   //Get the factor to artificially increase the cross section	
  	
	

	G4double CrossSectionPerVolume(const G4DynamicParticle* aDynamicGamma, const G4Material* aMaterial);
  
	
  	G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);	
  								//Computes and returns the final state of the	 process (at end of step),	
  								//returned as a ParticleChange object.	
  								//This function overloads a virtual function of the base class
  								// It is invoked by the ProcessManager of the Particle.	
  
	
	G4double ComputeCrossSectionPerAtom(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope);
	const G4Isotope * SelectRandomIsotope(const G4DynamicParticle* aDynamicGamma, const G4Material* material);

  	G4double ComputeMeanFreePath (const G4DynamicParticle* aDynamicGamma, const G4Material* aMaterial);	
	G4double GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*);
 
  private:
		
	NRFInteractions& operator=(const NRFInteractions &right);	
  	NRFInteractions(const	 NRFInteractions&);	

  	G4double LowestEnergyLimit;  				// low energy limit of the tables	
 	G4double HighestEnergyLimit;				// high energy limit of the tables	
  	G4double fminimalEnergy; 				// minimalEnergy of produced particles	
  	G4double CrossSecFactor;	 			// factor to artificially increase the cross section	
  
  	NRFCrossSectionDataStore* theCrossSectionDataStore;	

	std::vector<G4double> xsec;
	


};

#endif  


  	
  	
  	
  	


  	
  	
  	
  	
  	
  	
  	
  	
  


  	
  	


  	
  

























  




