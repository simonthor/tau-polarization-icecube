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
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>
#include <sstream>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Detector volume
  // from G4LogicalVolumeStore.

  // G4double envSizeXY = 0;
  // G4double envSizeZ = 0;

  // if (!fDetectorBox)
  // {
  //   G4LogicalVolume* envLV
  //     = G4LogicalVolumeStore::GetInstance()->GetVolume("DetectorLog");
  //   if ( envLV ) fDetectorBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  // }

  // if ( fDetectorBox ) {
  //   envSizeXY = fDetectorBox->GetXHalfLength()*2.;
  //   envSizeZ = fDetectorBox->GetZHalfLength()*2.;
  // }
  // else  {
  //   G4ExceptionDescription msg;
  //   msg << "Detector volume of box shape not found.\n";
  //   msg << "Perhaps you have changed geometry.\n";
  //   msg << "The gun will be place at the center.";
  //   G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
  //    "MyCode0002",JustWarning,msg);
  // }

  // G4double x0 = (G4UniformRand()-0.5) * envSizeXY/2;
  // G4double y0 = (G4UniformRand()-0.5) * envSizeXY/2;
  // G4double z0 = -0.55 * envSizeZ;
  // Open the file to read the tau leptons from
  // G4cout << "Event ID: " << anEvent->GetEventID() << G4endl;
  std::ifstream file;
  file.open(ffilename);

  // Read the line that corresponds to the event number (+1 because the first line is the header)
  std::string line;
  for (int i = 0; i < anEvent->GetEventID() + 2; ++i) {
    std::getline(file, line);
  }

  // G4cout << "Reading file" << line << G4endl;
  // Read the PDG ID and 4-momentum from the line, by splitting the line at each comma
  // The columns are (starting at column 0): column 1: PDG ID, column 2: E, column 3: px, column 4: py, column 5: pz
  std::string pdg_id_str;
  std::string E_str;
  std::string px_str;
  std::string py_str;
  std::string pz_str;
  std::stringstream ss(line);
  std::getline(ss, pdg_id_str, ','); // Skip the first number, which is the event number
  std::getline(ss, pdg_id_str, ',');
  std::getline(ss, E_str, ',');
  std::getline(ss, px_str, ',');
  std::getline(ss, py_str, ',');
  std::getline(ss, pz_str, ',');
  // Convert the strings to the correct data types
  int pdg_id = std::stoi(pdg_id_str);
  double E = std::stod(E_str) * GeV;
  double px = std::stod(px_str) * GeV;
  double py = std::stod(py_str) * GeV;
  double pz = std::stod(pz_str) * GeV;
  // Close the file
  file.close();
  // G4cout << "PDG ID: " << pdg_id << " E: " << E << " px: " << px << " py: " << py << " pz: " << pz << " MeV" << G4endl;

  // Set the particle gun position to (0,0,0)
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));

  // Get the particle with the correct PDG ID
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(15);
  fParticleGun->SetParticleDefinition(particle);
  
  // Set momentum. The energy is then inferred from the momentum and mass
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
  fParticleGun->SetParticleEnergy(E);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

