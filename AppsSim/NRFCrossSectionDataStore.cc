//**********************************************************************************
// GEANT4	 physics class-­‐-­‐>NRFCrossSectionDataStore: source file	
// N. Kikuzawa,	 05-­‐Feb-­‐07	
//****************************************************************  

#include "../include/NRFCrossSectionDataStore.hh"
#include "G4HadronicException.hh"	
#include <iostream>

using namespace std;

NRFCrossSectionDataStore::NRFCrossSectionDataStore(): NDataSetList(0), verboseLevel(0){

}


NRFCrossSectionDataStore::~NRFCrossSectionDataStore(){

}	

G4bool NRFCrossSectionDataStore::IsApplicable(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope){

	for (G4int i = NDataSetList-1; i >= 0; i--){
  		if (DataSetList[i]->IsApplicable(aDynamicGamma, anIsotope)){
  			return true;	
		}
	}

return false;	
}

G4double NRFCrossSectionDataStore::GetCrossSection(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope){

	if (NDataSetList == 0){
  		throw	 G4HadronicException(__FILE__, __LINE__,"NRFCrossSectionDataStore: no data sets registered");	 
  		return DBL_MIN;
	}

	for (G4int i = NDataSetList-1; i >= 0; i--) {
  		if (DataSetList[i]->IsApplicable(aDynamicGamma, anIsotope)){
  			G4double crossSection = DataSetList[i]->GetCrossSection(aDynamicGamma);
 			return crossSection;	
		}
	}

	throw G4HadronicException(__FILE__, __LINE__,"NRFCrossSectionDataStore: no applicable data set found for particle/element");	

	return DBL_MIN;	
}

G4ThreeVector  NRFCrossSectionDataStore::GetMomentumDirection(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope){
  
	G4double AtomicZ = anIsotope->GetZ();
  	G4double AtomicN = anIsotope->GetN();
 	G4ThreeVector aMomentumDirection = aDynamicGamma->GetMomentumDirection();
	
  	for (G4int i = NDataSetList-1; i>= 0; i--){
 		if (DataSetList[i]->GetZ() == AtomicZ && DataSetList[i]->GetN() == AtomicN){
  			aMomentumDirection = DataSetList[i]->GetMomentumDirection(aDynamicGamma);
  			return aMomentumDirection;
			}
	}

  	throw G4HadronicException(__FILE__, __LINE__,"NRFCrossSectionDataStore: no applicable data set found for particle/element");	
  	return  aMomentumDirection;	
}

void NRFCrossSectionDataStore::AddDataSet(NRFCrossSectionData* aDataSet){
	DataSetList.push_back( aDataSet );		
  	NDataSetList++;	
}


G4double NRFCrossSectionDataStore::GetBranching(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope, G4double){
	
	

	for (G4int i = NDataSetList-1; i >= 0; i--) {
		
		if (DataSetList[i]->IsApplicable(aDynamicGamma, anIsotope)){
	   		LastBranchingRatio = DataSetList[i]->theLastBranchingRatio;
			//FirstExcitedLevel = DataSetList[i]->theFirstExcitedLevel;
	  		return LastBranchingRatio;
		}

	}

	return 0;
}


G4double NRFCrossSectionDataStore::GetExcitation(const G4DynamicParticle* aDynamicGamma, const G4Isotope* anIsotope, G4double){

	for (G4int i = NDataSetList-1; i >= 0; i--) {
		if (DataSetList[i]->IsApplicable(aDynamicGamma, anIsotope)){
 	   		FirstExcitedLevel = DataSetList[i]->theFirstExcitedLevel;
 	  		return FirstExcitedLevel;
 		}

 	}

	return 0;
}



void NRFCrossSectionDataStore::BuildPhysicsTable(const G4ParticleDefinition& aParticleType){

 	if (NDataSetList > 0) {	
  		for (G4int i = 0; i < NDataSetList; i++) {	
  			DataSetList[i]->BuildPhysicsTable(aParticleType);
  		}
	}
}


void NRFCrossSectionDataStore::DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
	if (NDataSetList == 0) {	
		G4cout << "WARNING -­‐ NRFCrossSectionDataStore::DumpPhysicsTable: no data sets registered" << G4endl;	
		return;	  
	}

	for (G4int i = NDataSetList-1; i >= 0; i--) {
		DataSetList[i]->DumpPhysicsTable(aParticleType);
  	}
}

