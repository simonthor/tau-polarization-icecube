#include "SensitiveDetector.hh"

namespace B1
{

SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveDetector::~SensitiveDetector()
{
    // delete &outputFile;
}

void SensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
    int event_number = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    // Create a string for the file name 
    std::string filename = "output" + std::to_string(event_number) + ".csv";
    // Open file to write to
    outputFile.open(filename);
    outputFile << "eventNumber,PDGID,Energy(GeV),Px(GeV),Py(GeV),Pz(GeV)" << std::endl;
    G4cout << "Output file: " << filename << G4endl;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{
    // Close the file
    outputFile.close();
}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{   
    // aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    // Continue if the particle is not a tau neutrino
    if (aStep->GetTrack()->GetDefinition()->GetPDGEncoding() != 16) return true;
    // Get energy of particle and write to file
    const G4VTouchable *detector = aStep->GetPreStepPoint()->GetTouchable();
    // Write 4-momentum and PDG to file

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