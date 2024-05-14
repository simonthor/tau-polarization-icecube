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

#include "TrackingAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"


TrackingAction::TrackingAction(G4int pdg_id)
  : G4UserTrackingAction()
{
  fpdg_id = pdg_id;
}


void TrackingAction::PreUserTrackingAction(const G4Track *track) {
    int pdgID = track->GetDefinition()->GetPDGEncoding();
  // Check that the particle with an ID of 1 always is the beam particle
  if ((track->GetTrackID() == 1) && (pdgID != fpdg_id)) {
    G4cout << "Parent PDG ID: " << pdgID << G4endl;
    return;
  }

  // If the particle is the beam particle and it is just created, store its 4-momentum and position
  if ((track->GetParentID() == 0)) {
    if (pdgID != fpdg_id) {
      G4cerr << "Warning: the initial particle is not the beam particle\n";
    }

    auto momentum = track->GetMomentum();
    auto position = track->GetPosition();
    G4cout << "Event " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << ", "
          << "Beam particle: "
          // Write PDG ID
          << pdgID << ", p = "
          // Write 4-momentum
          << track->GetTotalEnergy() / GeV << "," << momentum.x()/GeV << "," << momentum.y()/GeV << "," << momentum.z()/GeV
          << ", x = " << position.x()/m << "," << position.y()/m << "," << position.z()/m << "\n";
  }

  // If the parent of the track is a tau lepton and the track is not a tau lepton, print out the track's information
  if (track->GetParentID() == 1 
      // && track->GetCreatorProcess()->GetProcessType() == fDecay
    //   && track->GetCurrentStepNumber() == 1 // There are multiple steps for each track/particle, so we only pick the first step
      && pdgID != fpdg_id) {
    
    auto momentum = track->GetMomentum();
    // Open the file in append mode
        // Write event number to file
        G4cout << "Event " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << ", "
          << "Creator process type: " << track->GetCreatorProcess()->GetProcessName() << ", "
          << "Child to beam particle: "
          // Write PDG ID
          << pdgID << ", p = "
          // Write 4-momentum
          << track->GetTotalEnergy() / GeV << "," << momentum.x()/GeV << "," << momentum.y()/GeV << "," << momentum.z()/GeV << "\n";
  } 
}
