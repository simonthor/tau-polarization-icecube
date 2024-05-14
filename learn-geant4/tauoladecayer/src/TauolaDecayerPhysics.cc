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

#include "TauolaDecayerPhysics.hh"
#include "TauolaDecayer.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Decay.hh"
#include "G4DecayTable.hh"

// factory
//
#include "G4PhysicsConstructorFactory.hh"
//
// register it with contructor factory
//
G4_DECLARE_PHYSCONSTR_FACTORY(TauolaDecayerPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TauolaDecayerPhysics::TauolaDecayerPhysics(G4int)
    : G4VPhysicsConstructor("TauolaDecayerPhysics")
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TauolaDecayerPhysics::~TauolaDecayerPhysics() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TauolaDecayerPhysics::ConstructParticle()
{
   // Nothing needs to be done here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TauolaDecayerPhysics::ConstructProcess()
{
    // Setting external decayer to tau leptons.
    // G4Decay will use the external decayer if G4Decay process is
    // assigned to an unstable particle and that particle does not
    // have its decay table.

    // Loop over all particles instantiated and remove already-assigned
    // decay table for tau's so that they will decay through
    // the external decayer (Tauola).

    // NOTE: The extDecayer will be deleted in G4Decay destructor

    // Define a new G4Decay object with an external decayer that will be assigned to the tau leptons
    TauolaDecayer* extDecayer = new TauolaDecayer();
    G4Decay* theDecayProcess = new G4Decay();
    theDecayProcess->SetExtDecayer(extDecayer);

    // Iterate over all particles defined in Geant4
    auto particleIterator=GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)())
    {    
        G4ParticleDefinition* particle = particleIterator->value();
        
        // If the particle is a tau lepton, remove its decay table and original decay process
        if ( std::abs(particle->GetPDGEncoding()) == 15 )
        {
            G4cout << "found particle " << particle->GetPDGEncoding() << "\n";
            particle->SetDecayTable(0);
            
            // Iterate over all processes via the process manager to find the decay process and remove it
            // The decay process is a shared object by all unstable particles without a defined decay table
            // Setting the external decayer on that object will therefore set the external decayer for all
            // unstable particles without decay tables, which is not what we want.
            G4ProcessManager *pmanager = particle->GetProcessManager();
            G4ProcessVector *pros = pmanager->GetProcessList();
            for (int i=0; i<pros->size(); ++i) {
                if ((*pros)[i]->GetProcessType() == fDecay) {
                    // G4cout << "Removing original decay process from " << particle->GetPDGEncoding() << "\n";
                    pmanager->RemoveProcess(i);
                    break;
                }
            }
            
            // Since we have removed the original decay process from the tau, we can now set a new one, 
            // which has an external decayer
            // G4cout << "Adding new decay process to " << particle->GetPDGEncoding() << "\n";
            pmanager->AddProcess(theDecayProcess);
            // idxPostStep and idxAtRest are global values in Geant4 that specify when the process should be used
            // For decays, idxPostStep and idxAtRest are the values that should be set (see https://github.com/strigazi/athena/blob/7057be3e1a6d385fbca8d043fd9ab02a9780fae6/Simulation/G4Extensions/RHadrons/src/RHadronsPhysicsTool.cxx#L208)
            pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
            pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
        }
    }

    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
