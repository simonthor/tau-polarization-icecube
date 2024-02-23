/**
 *
 * @author Simon Thor
 * @date 2024-02-22
 */

// Tauola headers
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMC3Event.h"
// #include "Tauola/TauolaHepMC3Particle.h"

// HepMC3 headers
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

using namespace std;
using namespace Tauolapp;
using namespace HepMC3;

int main(int argc, char **argv){

  if( argc<3 ) {
      std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> <output_file> [polx poly polz]" << std::endl;
      exit(-1);
  }

  vector<double> pol_vec = {};

  if (argc == 6) {
    pol_vec.push_back(atof(argv[3]));
    pol_vec.push_back(atof(argv[4]));
    pol_vec.push_back(atof(argv[5]));
    cout << "Manual polarization is used" << endl;
  }
  else {
    cout << "Polarization is inferred" << endl;
  }

  int events_parsed = 0;

  Tauola::initialize();

  ReaderAscii input_file (argv[1]);
  WriterAscii output_file (argv[2]);

  // for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
  while (!input_file.failed()) {
    GenEvent evt(Units::GEV,Units::MM);
    // Read event from input file
    input_file.read_event(evt);

    // If reading failed - exit loop
    if( input_file.failed() ) break;
    
    // cout << "BEFORE:" << endl;
    // Print::listing(evt);
    //Print::content(evt); // Prints detailed event content
    
    TauolaHepMC3Event t_event(&evt);
    if (pol_vec.size() == 0) {
      t_event.decayTaus();
    
    } else {
      // Iterate over all particles in the event, find the tau particle (abs(pdg_id) == 15) and decay it using decayOne
      for (auto p : evt.particles()) {
        if (abs(p->pdg_id()) == 15) {
          TauolaHepMC3Particle *htau = new TauolaHepMC3Particle(p);
          Tauola::decayOne(htau, false, pol_vec[0], pol_vec[1], pol_vec[2]);
          break;
        }
      }
    }
    // cout << "AFTER:" << endl;
    // Print::listing(evt);
    //Print::content(evt); // Prints detailed event content
    output_file.write_event(evt);
    
    ++events_parsed;
    if (events_parsed % 1000 == 0) {
      std::cout << "Parsed " << events_parsed << " events" << std::endl;
    }
  }
}

