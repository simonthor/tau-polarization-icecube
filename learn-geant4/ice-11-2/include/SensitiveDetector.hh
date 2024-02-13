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
    // Declare Initialize, which will be called at the beginning of each event
    void Initialize(G4HCofThisEvent * HCE);
    // Declare EndOfEvent, which will be called at the end of each event
    void EndOfEvent(G4HCofThisEvent * HCE);
    
private:
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
    // Create the file to write to
    std::ofstream outputFile;
};
}
#endif

