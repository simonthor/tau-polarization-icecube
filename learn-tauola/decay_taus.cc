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

int main(int argc, char **argv) {

  if( argc<4 ) {
      std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> <output_file> <polfile>" << std::endl;
      exit(-1);
  }

  ifstream polfile(argv[3]);
  // Read the first line of the file, which contains the column names
  string line;
  getline(polfile, line);
  // Split the line between , and determine which indices correspond to the polarization vector
  vector<int> pol_vec_i;
  stringstream ss(line);
  string token;
  int index_counter = 0;
  while (getline(ss, token, ',')) {
    if (token == "polx" || token == "poly" || token == "polz") {
      pol_vec_i.push_back(index_counter);
    }
    ++index_counter;
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
    
    // Iterate over all particles in the event, find the tau particle (abs(pdg_id) == 15) and decay it using decayOne
    for (auto p : evt.particles()) {
      if (abs(p->pdg_id()) == 15) {
        TauolaHepMC3Particle *htau = new TauolaHepMC3Particle(p);
        
        // Read a line in the csv file and get the polarization vector
        getline(polfile, line);
        stringstream ss(line);
        string token;
        vector<double> pol_vec;
        int index_counter = 0;

        while (getline(ss, token, ',')) {
          // Check that the event number in the csv file matches the event number in the HepMC3 GenEvent
          if (index_counter == 0) {
            if (stoi(token) != evt.event_number()) {
              cout << "Event number in csv file does not match event number in HepMC3 GenEvent\n"
              << "Instead of " << evt.event_number() << " found " << token << " in csv file\n";
              exit(-1);
            }
          }

          if (index_counter == pol_vec_i[0] || index_counter == pol_vec_i[1] || index_counter == pol_vec_i[2]) {
            pol_vec.push_back(stod(token));
          }
          ++index_counter;
        }
        
        Tauola::decayOne(htau, false, pol_vec[0], pol_vec[1], pol_vec[2]);
        break;
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

