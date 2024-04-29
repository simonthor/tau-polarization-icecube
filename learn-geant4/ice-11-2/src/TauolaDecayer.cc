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

#include "TauolaDecayer.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TauolaDecayer::simpleBoost(Tauolapp::TauolaParticle *tau, Tauolapp::TauolaParticle *target)
{
  double p1=tau->getPx();
  double p2=tau->getPy();
  double p3=tau->getPz();
  double E =tau->getE();
  double m =tau->getMass();

  double betx=p1/m;
  double bety=p2/m;
  double betz=p3/m;

  double gam=E/m;

  double pb=betx*target->getPx()+bety*target->getPy()+betz*target->getPz();

  target->setPx(target->getPx()+betx*(target->getE()+pb/(gam+1.0)));
  target->setPy(target->getPy()+bety*(target->getE()+pb/(gam+1.0)));
  target->setPz(target->getPz()+betz*(target->getE()+pb/(gam+1.0)));

  target->setE(target->getE()*gam+pb);
}

TauolaDecayer::TauolaDecayer()
   : G4VExtDecayer("TauolaDecayer")
{

   Tauolapp::Tauola::initialize();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TauolaDecayer::~TauolaDecayer()
{

   // delete fDecayer;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DecayProducts* TauolaDecayer::ImportDecayProducts(const G4Track& track)
{

   Tauolapp::TauolaHEPEVTEvent * evt = new Tauolapp::TauolaHEPEVTEvent();
   
   G4DecayProducts* dproducts = nullptr;   
   
   G4ParticleDefinition* pd = track.GetDefinition();
   int    pdgid   = pd->GetPDGEncoding();
   
   // check if pdgid is a tau or ant-tau
   if ( pdgid != 15 && pdgid != -15 )
   {
      G4cout << " Particle of pdgid = " << pdgid 
             << " is NOT a tau or anti-tau" << G4endl;
      return dproducts;
   }

   // NOTE: Energy should be in GeV 
   // TauolaHEPEVTParticle 	( 	int  	pdgid, int status, double px, double py, double pz, double e, double m, int ms, int me, int ds, int de )	
   Tauolapp::TauolaHEPEVTParticle *tau  = new Tauolapp::TauolaHEPEVTParticle(pdgid, 1, 
                                                         track.GetMomentum().x() / CLHEP::GeV, 
                                                         track.GetMomentum().y() / CLHEP::GeV,  
                                                         track.GetMomentum().z() / CLHEP::GeV,
                                                         track.GetDynamicParticle()->GetTotalEnergy() / CLHEP::GeV,
                                                         pd->GetPDGMass() / CLHEP::GeV,
                                                         -1, -1, -1, -1);
   evt->addParticle(tau);
   // If the polarization vector is defined in the lab frame, boost is set to true.
   // NOTE it might have the same effect if this is called only once, instead of before every decay. This is however unclear from the Tauola documentation
   if (boost) {
      Tauolapp::Tauola::setBoostRoutine(simpleBoost);
   }
   // specify polarization, if any
   // TODO read this from Geant4, from a file, or something. 
   // It is possible to set the polarization of particles in Geant4 (even though it is not used).
   // That might be possible to add in IceTray. Then, one can extract this information and calculate the decays

   double polx = 0.;
   double poly = 0.;
   double polz = 0.;
   // To have a fully left-handed polarization, uncomment the lines below and comment the three lines above
   // double polx = -tau->getPx() / tau->getP();
   // double poly = -tau->getPy() / tau->getP();
   // double polz = -tau->getPz() / tau->getP();
   
   Tauolapp::Tauola::decayOne(tau, 
      true, // This will first undecay the particle and then decay it again. In this case, it should not have an effect if it is true or false, since we manually defined the particle
      polx, poly, polz);
   
   int npart_after_decay = evt->getParticleCount();
   
   // create & fill up decay products
   //
   dproducts = new G4DecayProducts(*(track.GetDynamicParticle()));
   
   // create G4DynamicParticle out of each fDecayer->event entry (except the 1st one)
   // and push into dproducts
   
   // Iterate over the HepEvtEvent to extract all particles except the first one (which is the tau itself)
   G4cout << "Tau decays to ";
   for ( int ip=0; ip<npart_after_decay; ++ip )
   {
      Tauolapp::TauolaHEPEVTParticle* p = evt->getParticle(ip);
      // only select final state decay products (direct or via subsequent decays);
      // skip all others
      //
      // NOTE: in general, final state decays products will have 
      //       positive status code between 91 and 99 
      //       (in case such information could be of interest in the future)
      //
      // TODO check that this status check actually works
      if (p->getStatus() != 1 ) continue;
            
      G4ParticleDefinition* pddec = 
         G4ParticleTable::GetParticleTable()->FindParticle(p->getPdgID());
      if ( !pddec ) {
         G4cout << " Warning: particle with PDG ID " << p->getPdgID() << "could not be found in Geant4\n";
      }
      G4cout << p->getPdgID() << ", ";

      G4ThreeVector momentum = G4ThreeVector( p->getPx() * CLHEP::GeV,
                                              p->getPy() * CLHEP::GeV,
                                              p->getPz() * CLHEP::GeV ); 
      dproducts->PushProducts( new G4DynamicParticle( pddec, momentum) ); 
   }
   G4cout << "\n";
   return dproducts;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

