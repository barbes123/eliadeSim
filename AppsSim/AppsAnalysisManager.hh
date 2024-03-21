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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it
//
 
#ifndef APPSANALYSISMANAGER_HH
#define APPSANALYSISMANAGER_HH 

#include "globals.hh"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "RtypesCore.h"
#include <G4Types.hh>
#include "G4ThreadLocalSingleton.hh"
#include "G4Run.hh"


// Define the total number of columns in the ntuple
const G4int MaxNtCol = 9;

class AppsAnalysisManager
{ 
	friend class G4ThreadLocalSingleton< AppsAnalysisManager >;
public:
	
	virtual ~AppsAnalysisManager();
  
	static AppsAnalysisManager *Instance();
	
	void Write();
	void CloseFile();

	void CreateELIADETree();
	void FillELIADETreeColumn(G4int id, G4double value);
	void AddELIADETreeRow();
	
protected:
	AppsAnalysisManager();

private:
	static G4ThreadLocal AppsAnalysisManager* fgInstance;
	TFile* fFile{ nullptr };
	TTree* ftree{ nullptr };
	UChar_t Mod{' '};
	UChar_t Ch{' '};
	ULong64_t TimeStamp{0};
	Double_t FineTS{0.};
	UShort_t ChargeLong{0};
	UShort_t ChargeShort{0};
	UInt_t RecordLength{0};
	bool Signal{false};

};

#endif




