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
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include <fstream>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{
  // Read settings.yaml, find the energy and input it into the file name
  std::ifstream settingsFile("/home/simon/Code/icecube/icecube-tau-decays/settings.yaml");
  // Iterate over each line and if the first characters are "energy", split the line at the space and take the second element
  for (std::string line; std::getline(settingsFile, line); ) {
    if (line.substr(0, 8) == "energy: ") {
      out_filename = "/home/simon/Code/icecube/data/geant4_output_e" + line.substr(8) + ".csv";
      energy = std::stod(line.substr(8));
      tau_out_filename = "/home/simon/Code/icecube/data/geant4_tau_output_e" + line.substr(8) + ".csv";
      // G4cout << "Output file: " << out_filename << G4endl;
      break;
    }
    if (line == "g4_tau_info: on") {
      write_tau_info = true;
    }
  }
  settingsFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  int pdgID = step->GetTrack()->GetDefinition()->GetPDGEncoding();
  // Check that the particle with an ID of 1 always is a tau lepton
  if ((step->GetTrack()->GetTrackID() == 1) && (pdgID != 15)) {
    G4cout << "Parent PDG ID: " << pdgID << G4endl;
    return;
  }

  // If the parent of the track is a tau lepton and the track is not a tau lepton, print out the track's information
  if (step->GetTrack()->GetParentID() == 1 
      && step->GetTrack()->GetCurrentStepNumber() == 1 // There are multiple steps for each track/particle, so we only pick the first step
      && pdgID != 15) {
    
    auto momentum = step->GetTrack()->GetMomentum();
    // TODO create a RunAction that clears the file and adds a header
    fileMutex.lock();
    // Open the file in append mode
    std::ofstream outputFile(out_filename, std::ios_base::app);
    if (outputFile.is_open()) {
        // Write event number to file
        outputFile << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << ","
          // Write PDG ID
          << pdgID << ","
          // Write 4-momentum
          << step->GetTrack()->GetTotalEnergy() / GeV << "," << momentum.x()/GeV << "," << momentum.y()/GeV << "," << momentum.z()/GeV << "\n";
        
        // Close the file
        outputFile.close();
    } else {
        // Handle file open error
        G4cerr << "Error: Unable to open file for writing!" << G4endl;
    }
    // Unlock the mutex after writing to the file
    fileMutex.unlock();
    
    // Kill a particle to stop it from propagating further and/or decaying
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    
    // G4cout << "Track ID: " << step->GetTrack()->GetTrackID() << G4endl;
    // G4cout << "PDG ID: " << step->GetTrack()->GetDefinition()->GetPDGEncoding() << G4endl;
    // G4cout << "Energy: " << step->GetTrack()->GetTotalEnergy() / GeV << G4endl;
    // G4cout << "Momentum: " << step->GetTrack()->GetMomentum() / GeV << G4endl;
  }

  if (!write_tau_info) {
    return;
  }
  
  // If the particle is a tau and it is just created, store its 4-momentum and position
  if ((step->GetTrack()->GetParentID() == 0)  &&
       (step->GetTrack()->GetCurrentStepNumber() == 1)) {
    if (pdgID != 15) {
      G4cout << "Warning: the initial particle is not a tau lepton\n";
    }
    auto momentum = step->GetTrack()->GetMomentum();
    auto position = step->GetTrack()->GetPosition();
    fileMutex.lock();
    std::ofstream outputFile(tau_out_filename, std::ios_base::app);
    if (outputFile.is_open()) {
      outputFile << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << ","
        << pdgID << ","
        << step->GetTrack()->GetTotalEnergy() / GeV << "," << momentum.x()/GeV << "," << momentum.y()/GeV << "," << momentum.z()/GeV << ","
        << position.x()/m << "," << position.y()/m << "," << position.z()/m << ",0\n"; // The final 0 indicates that this is the first step in the tau's life
      outputFile.close();
    } else {
      G4cerr << "Error: Unable to open file for writing!" << G4endl;
    }
    fileMutex.unlock();
  }
  // If the tau is about to decay, store its 4-momentum and position 
  // If-statement taken from Hadr10 example: https://github.com/Geant4/geant4/blob/master/examples/extended/hadronic/Hadr10/src/SteppingAction.cc
  if ((step->GetTrack()->GetParentID() == 0)  &&
       (step->GetPostStepPoint()->GetProcessDefinedStep() != nullptr)  &&
       (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().find( "Decay" ) != std::string::npos) ) {
    if (pdgID != 15) {
      G4cout << "Warning: the initial particle is not a tau lepton\n";
    }
    auto momentum = step->GetTrack()->GetMomentum();
    auto position = step->GetTrack()->GetPosition();
    fileMutex.lock();
    std::ofstream outputFile(tau_out_filename, std::ios_base::app);
    if (outputFile.is_open()) {
      outputFile << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << ","
        << pdgID << ","
        << step->GetTrack()->GetTotalEnergy() / GeV << "," << momentum.x()/GeV << "," << momentum.y()/GeV << "," << momentum.z()/GeV << ","
        << position.x()/m << "," << position.y()/m << "," << position.z()/m << ",1\n"; // The final 1 indicates that this is the last step in the tau's life
      outputFile.close();
    } else {
      G4cerr << "Error: Unable to open file for writing!" << G4endl;
    }
    fileMutex.unlock();
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}