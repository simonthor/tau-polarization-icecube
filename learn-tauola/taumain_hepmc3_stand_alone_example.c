/**
 * Example of use of tauola C++ on pre-constructed HepMC3 event.
 * HepMC3 events  e+, e- -> tau+ tau- are constructed.
 * Taus are subsequently decayed via tauola.
 * 
 * Note that this is a simple example to demonstrate HepMC3
 * integration. For more details on Tauola use
 * see other examples.
 *
 * @author Tomasz Przedzinski
 * @date 15 January 2020
 */

#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMC3Event.h"

// HepMC3 headers
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"

using namespace std;
using namespace Tauolapp;
using namespace HepMC3;

/* Create simple e+ e- -> tau+ tau- event */
GenEvent make_simple_tau_event() {

  // Create event
  GenEvent evt(Units::GEV,Units::CM);

  // Create four vectors for the electrons
  double ele_mass_sqr = parmas_.amell*parmas_.amell;
  FourVector momentum_e1(0,0,0,0);
  FourVector momentum_e2(0,0,0,0);

  momentum_e1.setPz(-2);  //change these
  momentum_e2.setPz(3.5); //as needed

  momentum_e1.setE(sqrt(momentum_e1.pz()*momentum_e1.pz()+ele_mass_sqr));
  momentum_e2.setE(sqrt(momentum_e2.pz()*momentum_e2.pz()+ele_mass_sqr));

  // Create four vectors for tau taus
  double tau_mass = parmas_.amtau;
  FourVector momentum_tau1(0,0,0,0);
  FourVector momentum_tau2(0,0,0,0);

  // Make particles
  GenParticlePtr e1   = make_shared<GenParticle>(momentum_e1,   -11, 3);
  GenParticlePtr e2   = make_shared<GenParticle>(momentum_e2,    11, 3);
  GenParticlePtr tau1 = make_shared<GenParticle>(momentum_tau1, -15, 1);
  GenParticlePtr tau2 = make_shared<GenParticle>(momentum_tau2,  15, 1);
  
  // Set the masses
  e1->set_generated_mass(parmas_.amell);
  e2->set_generated_mass(parmas_.amell);
  tau1->set_generated_mass(tau_mass);
  tau2->set_generated_mass(tau_mass);

  // Make the vertex
  GenVertexPtr v1 = make_shared<GenVertex>();
  v1->add_particle_in (e1);
  v1->add_particle_in (e2);
  v1->add_particle_out(tau1);
  v1->add_particle_out(tau2);

  evt.add_vertex(v1);

  // Calculate center of mass frame
  FourVector cms(0,0,(momentum_e1.pz()+momentum_e2.pz()),
			  momentum_e1.e()+momentum_e2.e());

  GenParticlePtr cms_particle = make_shared<GenParticle>(cms,0,0);

  // Make TauolaParticles for boosting
  TauolaHepMC3Particle cms_boost(cms_particle);

  TauolaHepMC3Particle first_tau(tau1);
  TauolaHepMC3Particle second_tau(tau2);

  double tau_energy = 0.5*sqrt(cms.e()*cms.e() - (cms.px()*cms.px()
			+ cms.py()*cms.py()+cms.pz()*cms.pz()));
		
  first_tau.setE(tau_energy);
  first_tau.setPx((1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass*tau_mass));
  first_tau.setPy((1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass*tau_mass));

  second_tau.setE(tau_energy);
  second_tau.setPx(-1*(1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass*tau_mass));
  second_tau.setPy(-1*(1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass*tau_mass));

  first_tau.boostFromRestFrame(&cms_boost);
  second_tau.boostFromRestFrame(&cms_boost);

  return evt;
}

int main(void){

  int NumberOfEvents = 10;

  Tauola::initialize();

  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {

    GenEvent evt = make_simple_tau_event();

    cout << "BEFORE:" << endl;
    Print::listing(evt);
    //Print::content(evt); // Prints detailed event content
    
    TauolaHepMC3Event t_event(&evt);
    t_event.decayTaus();

    cout << "AFTER:" << endl;
    Print::listing(evt);
    //Print::content(evt); // Prints detailed event content
  }
}

