/// \file SensitiveDetector.hh
/// \brief Definition of the B1::SensitiveDetector class

#ifndef B1SensitiveDetector_h
#define B1SensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VTouchable.hh"


namespace B1
{

class SensitiveDetector : public G4VSensitiveDetector
{
public:
    SensitiveDetector(G4String);
    ~SensitiveDetector();
    
private:
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
};
}
#endif

