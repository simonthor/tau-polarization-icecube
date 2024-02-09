/// \file PhysicsList.cc
/// \brief Implementation of the B1::PhysicsList class

#include "PhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4EmExtraPhysics.hh"
// TODO try to copy the relevant files related to G4NeutrinoPhysics and compile together with the rest
// #include "G4NeutrinoPhysics.hh"
// #include "G4HadronPhysicsNuBeam.hh"

namespace B1
{

PhysicsList::PhysicsList()
{
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
  RegisterPhysics(new G4EmStandardPhysics());
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  G4EmExtraPhysics* nuphys = new G4EmExtraPhysics();
  nuphys->NeutrinoActivated(true);
  nuphys->NuETotXscActivated(true);
  nuphys->SetNuDetectorName("DetectorLog");
  nuphys->SetNuNucleusBias(1e34);
  RegisterPhysics(nuphys);
  // RegisterPhysics(new G4HadronPhysicsNuBeam());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
