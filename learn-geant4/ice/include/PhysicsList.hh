/// \file PhysicsList.hh
/// \brief Definition of the B3::PhysicsList class

#ifndef B1PhysicsList_h
#define B1PhysicsList_h 1

#include "G4VModularPhysicsList.hh"

namespace B1
{

/// Modular physics list
///
/// It includes the folowing physics builders
/// - G4DecayPhysics
/// - G4RadioactiveDecayPhysics
/// - G4EmStandardPhysics

class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList();
  ~PhysicsList() override;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
