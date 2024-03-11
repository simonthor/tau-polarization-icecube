/**
 *
 * @author Simon Thor
 * @date 2024-02-22
 */

// Tauola headers
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMC3Event.h"
#include "Tauola/Log.h"

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


vector<double> get_pol_vec(string& line, vector<int>& pol_vec_i, GenEvent& evt) {
  // Read a line in the csv file and get the polarization vector
  stringstream ss(line);
  string token;
  vector<double> pol_vec;
  int index_counter = 0;

  while (getline(ss, token, ',')) {
    // cout << "Token: " << token << endl;
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

  return pol_vec;
}

int main(int argc, char **argv) {
  ifstream polfile;
  bool polfile_given = false;
  vector<double> pol_vec;
  vector<int> pol_vec_i;
  string line;
  // No limit for the number of warnings
  Log::SetWarningLimit(0);

  if( argc<4 ) {
      std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> <output_file> [polx poly polz] [polfile]" << std::endl;
      exit(-1);
  }
  else if (argc == 7) {
    polfile_given = true;
    polfile = ifstream(argv[6]);
    // If the file cannot be opened, print an error message and exit
    if (!polfile.is_open()) {
      cout << "Could not open file " << argv[3] << endl;
      exit(-1);
    }
    pol_vec_i = {stoi(argv[3]), stoi(argv[4]), stoi(argv[5])};
    // Skip first line
    getline(polfile, line);
    // cout << "Polarization vector indices: " << pol_vec_i[0] << " " << pol_vec_i[1] << " " << pol_vec_i[2] << endl;
  }
  else {
    pol_vec = {stod(argv[3]), stod(argv[4]), stod(argv[5])};
  }

  int events_parsed = 0;

  Tauola::initialize();
  // TODO disable radiative correction?
  
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
        
        if (polfile_given) {
          getline(polfile, line);
          pol_vec = get_pol_vec(line, pol_vec_i, evt);
          // If the polarization vector has a norm > 1, normalize it to have a norm of 1. 
          // This is caused by floating point errors
          double norm = pol_vec[0]*pol_vec[0] + pol_vec[1]*pol_vec[1] + pol_vec[2]*pol_vec[2];
          while (norm > 1) {
            pol_vec[0] /= norm;
            pol_vec[1] /= norm;
            pol_vec[2] /= norm;
            norm = pol_vec[0]*pol_vec[0] + pol_vec[1]*pol_vec[1] + pol_vec[2]*pol_vec[2];
          }
        }
        // Decay the particle with the specified polarization
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

