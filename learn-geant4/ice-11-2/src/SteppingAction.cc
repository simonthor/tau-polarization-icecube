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

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // If the parent of the track is a tau lepton and the track is not a tau lepton, print out the track's information
  if (step->GetTrack()->GetParentID() == 1 
      && step->GetTrack()->GetCurrentStepNumber() == 1 // There are multiple steps for each track/particle, so we only pick the first step
      && step->GetTrack()->GetDefinition()->GetPDGEncoding() != 15) {
    G4cout << "Track ID: " << step->GetTrack()->GetTrackID() << G4endl;
    // G4cout << "Track step number: " << step->GetTrack()->GetCurrentStepNumber() << G4endl;
    G4cout << "PDG ID: " << step->GetTrack()->GetDefinition()->GetPDGEncoding() << G4endl;
    G4cout << "Energy: " << step->GetTrack()->GetKineticEnergy() << G4endl;
    G4cout << "Momentum: " << step->GetTrack()->GetMomentum() << G4endl;
  }
  // if ( theStep->GetTrack()->GetParentID() == 0  &&
  //     theStep->GetTrack()->GetCurrentStepNumber() == 1 ) {
  // }
  // if ( theStep->GetTrack()->GetParentID() == 0  &&
  //       theStep->GetPostStepPoint()->GetProcessDefinedStep() != nullptr  &&
  //       theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().find( "Decay" )
  //       != std::string::npos ) {
  // }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}