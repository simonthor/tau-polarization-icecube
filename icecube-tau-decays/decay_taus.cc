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


/*
   Simple boost routine example.
   Note that this routine will not respect direction of polarization vector along boost direction,
   thus (0,0,1) does not mean helicity state. This is simply z component of spin for tau
   with momentum px,py,pz. Can be thus helicitywise anything in range (-1,1).
*/
void simpleBoost(TauolaParticle *tau, TauolaParticle *target)
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

  if(argc <= 3) {
      std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> <output_file> [polx poly polz] [polfile] [-r]\n"
      << "HepMC3_input_file : The HepMC3 input file name with tau leptons that have not been decayed.\n"
      << "                    The tau lepton status should be 1.\n"
      << "output_file : the HepMC3 output file name with tau leptons decayed.\n"
      << "polx poly polz : The polarization vector of the tau lepton.\n"
      << "                 If the polfile is not given, the polarization vector is set to these values.\n"
      << "                 If polfile is given, these values are used as the column indices for the polarization vector in the csv file.\n"
      << "polfile : The csv file with the polarization vectors for the tau leptons.\n" 
      << "          The first row should contain the column names.\n"
      << "-r : If this flag is given, radiative corrections are turned off.\n"
      << "-b : Set boost routine to simpleBoost. This is used if the input tau polarization is in the lab frame.\n";

      exit(-1);
  }

  bool radiation = true;
  bool set_boost = false;
  // If any of the arguments are -r, set radiative corrections to false and remove it from argv
  for (int i = 3; i < argc; ++i) {
    bool argfound = false;

    if (string(argv[i]) == "-r") {
      radiation = false;
      argfound = true;
      cout << "Radiative corrections are turned off" << std::endl;
      break;
    }
    else if (string(argv[i]) == "-b") {
      cout << "Boost routine is activated" << std::endl;
      set_boost = true;
      argfound = true;
    }

    if (argfound) {
      for (int j = i; j < argc-1; ++j) {
          argv[j] = argv[j+1];
        }
        --argc;
    }
  }

  // Print argv
  for (int i = 0; i < argc; ++i) {
    cout << "argv[" << i << "]: " << argv[i] << endl;
  }
  if (argc == 7) {
    polfile_given = true;
    polfile = ifstream(argv[6]);
    // If the file cannot be opened, print an error message and exit
    if (!polfile.is_open()) {
      cout << "Could not open file " << argv[6] << endl;
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

  // disable radiative correction if -r flag is given
  Tauola::setRadiation(radiation);

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
          // Extract event number. It is all text in the line before the first comma
          // int event_num = stoi(line.substr(0, line.find(',')));
          // cout << "event_num: " << event_num << " pol_vec: " << pol_vec[0] << " " << pol_vec[1] << " " << pol_vec[2] << endl;
          while (norm > 1) {
            pol_vec[0] /= norm;
            pol_vec[1] /= norm;
            pol_vec[2] /= norm;
            norm = pol_vec[0]*pol_vec[0] + pol_vec[1]*pol_vec[1] + pol_vec[2]*pol_vec[2];
          }
          // cout << "pol_vec after normalization: " << pol_vec[0] << " " << pol_vec[1] << " " << pol_vec[2] << endl;
        }
        if (set_boost) {
          Tauola::setBoostRoutine(simpleBoost);
        }
        // Decay the particle with the specified polarization
        Tauola::decayOne(htau, true, pol_vec[0], pol_vec[1], pol_vec[2]);
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

