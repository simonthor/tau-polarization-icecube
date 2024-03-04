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
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
  : G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // Iterate over all particles in the event, not just the primary vertex. 
  // Find the particles that have a mother particle of a tau lepton, but are not a tau lepton themselves
  // and print out their information (PDG ID, E, px, py, pz)
  // G4TrajectoryContainer* trajectoryContainer = event->GetStepCon();
  // if (!trajectoryContainer) {
  //     G4cerr << "Error: No trajectory container found for the event!" << G4endl;
  //     return;
  // }

  // for (size_t i = 0; i < trajectoryContainer->size(); ++i) {
  //   G4VTrajectory* trajectory = (*trajectoryContainer)[i];
  //   if (trajectory) {
  //       const G4Track* track = trajectory->GetPrimaryTrack();
  //       const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();

  //       // Check if the parent particle is a tau lepton
  //       if (track->GetParentID() == 1) { // Assuming tau lepton's parent ID is 1
  //           const G4ParticleDefinition* parentDefinition = track->GetParent()->GetDefinition();
  //           if (parentDefinition->GetPDGEncoding() == 15) { // PDG ID for tau lepton
  //               // Exclude particles that are tau leptons themselves
  //               if (particleDefinition->GetPDGEncoding() != 15) {
  //                   // Print information about daughter particles
  //                   G4int pdgID = particleDefinition->GetPDGEncoding();
  //                   G4double energy = track->GetKineticEnergy();
  //                   G4double time = track->GetGlobalTime();
  //                   G4cout << "Particle PDG ID: " << pdgID << G4endl;
  //                   G4cout << "Particle Kinetic Energy: " << energy << " MeV" << G4endl;
  //                   G4cout << "Particle Global Time: " << time << " ns" << G4endl;
  //                   G4cout << G4endl;
  //               }
  //           }
  //       }
  //   } else {
  //       G4cerr << "Error: No trajectory found at index " << i << G4endl;
  //   }
  // }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

}