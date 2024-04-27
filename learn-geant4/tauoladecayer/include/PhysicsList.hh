/// \file PhysicsList.hh
/// \brief Definition of the B3::PhysicsList class

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"


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

#endif
