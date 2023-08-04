 #include "globals.hh"
 #include "G4VEmProcess.hh"
 #include "G4Gamma.hh"
 

 class G4ParticleDefinition;
 class G4VEmModel;
 class G4MaterialCutsCouple;
 class G4DynamicParticle;
 class G4PhotoElastic : public G4VEmProcess
 
 {
 public:
 
   G4PhotoElastic(const G4String& processName ="PhotoElastic",
                       G4ProcessType type = fElectromagnetic);
 
   virtual ~G4PhotoElastic();
 
   // true for Gamma only.  
   virtual G4bool IsApplicable(const G4ParticleDefinition&);
   virtual void PrintInfo();
 
 protected:
 
   virtual void InitialiseProcess(const G4ParticleDefinition*);
 
 private:
      
   G4bool       isInitialised;
 };
