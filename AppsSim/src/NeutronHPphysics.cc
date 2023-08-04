//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr04/src/NeutronHPphysics.cc
/// \brief Implementation of the NeutronHPphysics class
//
// $Id: NeutronHPphysics.cc 66587 2012-12-21 11:06:44Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NeutronHPphysics.hh"

//#include "NeutronHPMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Processes

#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPInelastic.hh"

#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPCapture.hh"

#include "G4NeutronFissionProcess.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPFission.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::NeutronHPphysics(const G4String& name)
:  G4VPhysicsConstructor(name)
{
//  fNeutronMessenger = new NeutronHPMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::~NeutronHPphysics()
{
//  delete fNeutronMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHPphysics::ConstructProcess()
{
   G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();
   
   // process: elastic
   //
   G4HadronElasticProcess* process1 = new G4HadronElasticProcess();
   pManager->AddDiscreteProcess(process1);   
   //
   // cross section data set
   G4NeutronHPElasticData* dataSet1a = new G4NeutronHPElasticData();
   G4NeutronHPThermalScatteringData* dataSet1b 
                               = new G4NeutronHPThermalScatteringData();
   process1->AddDataSet(dataSet1a);                               
   process1->AddDataSet(dataSet1b);
   //
   // models
   G4NeutronHPElastic*           model1a = new G4NeutronHPElastic();
   G4NeutronHPThermalScattering* model1b = new G4NeutronHPThermalScattering();
   model1a->SetMinEnergy(4*eV);
   process1->RegisterMe(model1a);
	 model1b->SetMaxEnergy(4*eV);
   process1->RegisterMe(model1b);
   
   // process: inelastic
   //
   G4HadronInelasticProcess* process2 = new G4HadronInelasticProcess("neutron-Inelastic",G4Neutron::Definition());
   pManager->AddDiscreteProcess(process2);   
   //
   // cross section data set
   G4NeutronHPInelasticData* dataSet2 = new G4NeutronHPInelasticData();
   process2->AddDataSet(dataSet2);                               
   //
   // models
   G4NeutronHPInelastic* model2 = new G4NeutronHPInelastic();
   process2->RegisterMe(model2);    

   // process: nCapture   
   //
   G4NeutronCaptureProcess* process3 = new G4NeutronCaptureProcess();
   pManager->AddDiscreteProcess(process3);    
   //
   // cross section data set
   G4NeutronHPCaptureData* dataSet3 = new G4NeutronHPCaptureData();
   process3->AddDataSet(dataSet3);                               
   //
   // models
   G4NeutronHPCapture* model3 = new G4NeutronHPCapture();
   process3->RegisterMe(model3);
   
   // process: nFission   
   //
   G4NeutronFissionProcess* process4 = new G4NeutronFissionProcess();
   pManager->AddDiscreteProcess(process4);    
   //
   // cross section data set
   G4NeutronHPFissionData* dataSet4 = new G4NeutronHPFissionData();
   process4->AddDataSet(dataSet4);                               
   //
   // models
   G4NeutronHPFission* model4 = new G4NeutronHPFission();
   process4->RegisterMe(model4);       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
