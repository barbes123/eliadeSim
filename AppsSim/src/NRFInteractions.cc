//********************************************************************
// GEANT4	 Physics-­‐-­‐>NRFInteractions	 class:	 source	 file
//Class	 Description	
//This  class describes the	photon	interaction with	 
// NRF interaction with  its cross  section or the  atomic interactions with its cross	sections	  
//As well the nuclear	recoil	after	the	NRF	interaction	 -­‐	Hani	Negm	 2012-­‐2014	 
// ********************************************************************	
	
#include "../include/NRFInteractions.hh"
#include "G4EnergyLossTables.hh"	
#include "G4UnitsTable.hh"	
#include "G4SystemOfUnits.hh"	
#include "G4ios.hh"	
#include "G4ElementTable.hh"
  
using namespace std;

NRFInteractions::NRFInteractions(const G4String&  processName, G4ProcessType type)
: G4VDiscreteProcess (processName, type), CrossSecFactor(1.)
{	  
   	
  	theCrossSectionDataStore = new NRFCrossSectionDataStore();

  	G4int numEle = G4Element::GetNumberOfElements();		

	G4ElementTable * elements = G4Element::GetElementTable();
	
  
	for (G4int i = 0; i< numEle; i++){	
		
		G4IsotopeVector * isotopes = (*elements)[i]->GetIsotopeVector();
		G4int numIso = (*elements)[i]->GetNumberOfIsotopes();

		for(G4int o = 0; o < numIso; o++){
		
			NRFCrossSectionData  * aGammaData = new	NRFCrossSectionData((*isotopes)[o]);	
		  	theCrossSectionDataStore->AddDataSet(aGammaData);

		}

  	}


	

  
}	
  
NRFInteractions::~NRFInteractions(){
	delete theCrossSectionDataStore;
}

G4bool NRFInteractions::IsApplicable(const G4ParticleDefinition& particle){

	  return (&particle == G4Gamma::Gamma());

}

void NRFInteractions::BuildPhysicsTable(const G4ParticleDefinition&){	
  
	PrintInfoDefinition();

}


	 
G4double NRFInteractions::CrossSectionPerVolume(const G4DynamicParticle* aDynamicGamma,const G4Material* aMaterial){

	G4double cross = 0.0;

  	const G4ElementVector* theElementVector = aMaterial->GetElementVector();
	const G4double* theAtomNumDensityVector = aMaterial->GetVecNbOfAtomsPerVolume();
  	G4int nelm = aMaterial->GetNumberOfElements(); 
	
	xsec.clear();

	for(int i = 0; i < nelm; i++){
	
		const G4IsotopeVector* theIsotopeVector  = (*theElementVector)[i]->GetIsotopeVector();
		G4int niso = (*theElementVector)[i]->GetNumberOfIsotopes();
		const G4double * FractionVector = (*theElementVector)[i]->GetRelativeAbundanceVector();

		for(int o = 0; o < niso; o++){
 		
			cross += theAtomNumDensityVector[i]
			*ComputeCrossSectionPerAtom(aDynamicGamma, (*theIsotopeVector)[o])*FractionVector[o];
			xsec.push_back(cross);

		}
	
	}
	
	return cross;


}

	
G4double NRFInteractions::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*){ 

	const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();	 	
  	G4Material* aMaterial = aTrack.GetMaterial();	
 	G4double MeanFreePath = ComputeMeanFreePath(aDynamicGamma, aMaterial);	
	
  	return MeanFreePath;	 
}

// Computes and returns the photon mean free path in GEANT4 internal units	
G4double NRFInteractions::ComputeMeanFreePath(const G4DynamicParticle* aDynamicGamma,const G4Material* aMaterial){

	G4double cross = CrossSectionPerVolume(aDynamicGamma, aMaterial);
  	G4double MeanFreePath = DBL_MAX;

	if(cross > DBL_MIN){
		MeanFreePath = 1./cross;	
	} 

 	return	MeanFreePath;	
}

  
//Calculates the microscopic cross section in GEANT4 internal units.	
// Total cross section parametrisation from H.Burkhardt	
// It gives a good description at any energy (from 0 to 10**21  eV)


G4double NRFInteractions::ComputeCrossSectionPerAtom(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope){  
 	

	G4double crossSection = theCrossSectionDataStore->GetCrossSection(aDynamicGamma,anIsotope);
	crossSection *= CrossSecFactor;	//increase  the	CrossSection by	 (by default 1)	


	return crossSection;
}


//  Set	 the factor to artificially increase the cross section	
void NRFInteractions::SetCrossSecFactor(G4double fac){
 	
	CrossSecFactor = fac;

 	G4cout << "The cross section for NRF Interaction is artificially increased by the CrossSecFactor= " << CrossSecFactor << G4endl;	
}


G4VParticleChange* NRFInteractions::PostStepDoIt(const G4Track& aTrack, const	G4Step& aStep){
 	

	aParticleChange.Initialize(aTrack);		
  	G4Material* aMaterial = aTrack.GetMaterial();	
 
	// current Gamma energy and direction
 	const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();	
  	G4double Egam = aDynamicGamma->GetKineticEnergy();
  	G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();

	//select randomly one element constituting the material	

	const G4Isotope* anIsotope = SelectRandomIsotope(aDynamicGamma, aMaterial);

    	if(!theCrossSectionDataStore->IsApplicable(aDynamicGamma, anIsotope)){
  	 	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);	
  	}
 
	G4ParticleMomentum aMomentumDirection =	theCrossSectionDataStore->GetMomentumDirection(aDynamicGamma, anIsotope);

	//Calculate the	 recoil energy of interested isotope	
 
  	G4double AtomicMass = anIsotope->GetA()/CLHEP::g/CLHEP::mole;
  	G4double AtomicMassUnit = 931.494 *MeV;
  	G4double energy = aDynamicGamma->GetKineticEnergy()/CLHEP::MeV;
  	G4double E_Recoil = energy*energy/(2. *AtomicMass * AtomicMassUnit);
  	G4double v = G4UniformRand();
  	G4double beta = theCrossSectionDataStore->GetBranching(aDynamicGamma, anIsotope);
  	G4double excitation = theCrossSectionDataStore->GetExcitation(aDynamicGamma, anIsotope);


  	if (v<=beta) {

	  	G4double NRF_Energy = energy - E_Recoil ;
	  	Egam = NRF_Energy;

 	} else {

  		G4double NRF_Energy = energy - E_Recoil - excitation *MeV;
		Egam = NRF_Energy;
	}

 
	//NRF gamma ray properties:	
  	G4DynamicParticle* aParticle = new G4DynamicParticle(G4Gamma::Gamma(), aMomentumDirection, Egam);	
	

 	aParticleChange.SetNumberOfSecondaries(1);	
 	aParticleChange.AddSecondary(aParticle);	

  	aParticleChange.ProposeMomentumDirection(0., 0., 0.);	
  	aParticleChange.ProposeEnergy(0.);	
  	aParticleChange.ProposeTrackStatus(fStopAndKill);	




    	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}




const G4Isotope * NRFInteractions::SelectRandomIsotope(const G4DynamicParticle* aDynamicGamma, const G4Material* aMaterial){
  
  	const G4int numElem = aMaterial->GetNumberOfElements();
 	const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  	
	if(numElem == 1){

		G4int numIso = (*theElementVector)[0]->GetNumberOfIsotopes();

		if(numIso == 1){
			return (*theElementVector)[0]->GetIsotope(0);	
		} else {
		  	
  			G4double x = G4UniformRand()*CrossSectionPerVolume(aDynamicGamma, aMaterial);
			for(G4int i = 0; i < numIso; i++){
				if(x <= xsec[i]){
					return (*theElementVector)[0]->GetIsotope(i);
				}	
			}

		}

	} else {

		G4double x = G4UniformRand()*CrossSectionPerVolume(aDynamicGamma, aMaterial);
		G4int cnt = 0;		

		for(G4int i = 0; i < numElem; i++){
			
			G4int numIso = (*theElementVector)[i]->GetNumberOfIsotopes();

			for(G4int o = 0; o < numIso; o++){
				if(x <= xsec[cnt]){
					return (*theElementVector)[i]->GetIsotope(o);
				}
				cnt++;
			}	
		}
	}

	
	return NULL;	
}

void NRFInteractions::PrintInfoDefinition(){

 	G4cout <<G4endl<< GetProcessName() << ": "<< "gamma-­‐>gamma process"	   /* << GetProcessSubType()	 */<<G4endl;	
  
}















  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  	
  
  


