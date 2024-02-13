#include "SensitiveDetector.hh"

namespace B1
{

SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveDetector::~SensitiveDetector()
{}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{   
    // aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    // Continue if the particle is not a tau neutrino
    // if (aStep->GetTrack()->GetDefinition()->GetPDGEncoding() != 16) return true;
    // Get energy of particle and write to file
    const G4VTouchable *detector = aStep->GetPreStepPoint()->GetTouchable();
    // Write 4-momentum and PDG to file
    #ifdef G4MULTITHREADED
    static G4Mutex mutex = G4MUTEX_INITIALIZER;
    G4AutoLock al(&mutex);
    #endif
    static std::ofstream outputFile("hits.csv");
    static bool first = true;
    if (first) {
        first = false;
        outputFile << "eventNumber,PDGID,Energy(GeV),Px(GeV),Py(GeV),Pz(GeV)" << std::endl;
    };
    outputFile << detector->GetCopyNumber() << ","
        << aStep->GetTrack()->GetDefinition()->GetPDGEncoding() << ","
        << aStep->GetPreStepPoint()->GetTotalEnergy() / GeV  << ","
        << aStep->GetPreStepPoint()->GetMomentum().x() / GeV << ","
        << aStep->GetPreStepPoint()->GetMomentum().y() / GeV << ","
        << aStep->GetPreStepPoint()->GetMomentum().z() / GeV
        << "\n";
    
    return true;
}
}