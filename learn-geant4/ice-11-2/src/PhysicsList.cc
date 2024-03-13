/// \file PhysicsList.cc
/// \brief Implementation of the B1::PhysicsList class

#include "PhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
// #include "G4RadioactiveDecayPhysics.hh"
// #include "G4HadronPhysicsFTFP_BERT.hh"
// #include "G4EmExtraPhysics.hh"
// #include "G4NeutrinoPhysics.hh"

namespace B1
{

PhysicsList::PhysicsList()
{
  RegisterPhysics(new G4DecayPhysics());
  // RegisterPhysics(new G4RadioactiveDecayPhysics());
  RegisterPhysics(new G4EmStandardPhysics());
  // RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  // RegisterPhysics(new G4NeutrinoPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
