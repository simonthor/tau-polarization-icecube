/* @author Simon Thor
 * Read a ROOT file generated by GENIE.
 * Take the incoming nu_tau, incoming oxygen, outgoing tau.
 * Store the 4-momentum and polarization of all of these particles in a csv file.
 * Ignore all events with missing tau polarization.
*/
#include <string>
#include <iostream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/CmdLnArgParser.h"

using std::string;
using namespace genie;

// void GetCommandLineArgs (int argc, char ** argv);

int    n_events = -1;
string input_filename = "../data/gntp.0.ghep.root";
string output_filename = "../data/genie_tau_pol_data.csv";
//___________________________________________________________________
int get_pol_genie() {
  
  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(input_filename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (n_events > 0) ?
        TMath::Min(n_events, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  // Open a csv file where the interesting data will be stored
  std::ofstream csv_file;
  csv_file.open(output_filename);
  
  csv_file << "event_num,pdg,E,px,py,pz,polx,poly,polz\n";

  //
  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);
    // Print the whole event:
    // LOG("myAnalysis", pNOTICE) << event;

    GHepParticle * p = 0;
    TIter event_iter(&event);

    std::vector<GHepParticle*> interesting_particles;
    interesting_particles.clear();

    // Number of taus in an event
    int n_taus = 0;


    //
    // Loop over all particles in this event
    //
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
      if (
        // If the particle is a tau
        ((p->Status() == kIStStableFinalState) && (p->Pdg() == kPdgTau) && (p->PolzIsSet()))
        // If the particle is an incoming particle
        || (p->Status() == kIStInitialState)
      ) {
        interesting_particles.push_back(p);
        if (p->Pdg() == kPdgTau) {
          n_taus++;
        }
      }

    }// end loop over particles	
    
    // If there are no taus in the event, skip it
    if (n_taus != 1) {
      continue;
    }

    for (auto p : interesting_particles) {
      TVector3 pol;
      p->GetPolarization(pol);
      // Write the 4-momentum and polarization of the particle to the csv file
      csv_file << i << "," << p->Pdg() << "," 
        << p->E() << "," << p->Px() << "," << p->Py() << "," << p->Pz() << "," 
        << pol.X() << "," << pol.Y() << "," << pol.Z() << "\n";
    }

    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  // Print total number of events
  LOG("myAnalysis", pNOTICE) << "Processed " << nev << " events";

  // close input GHEP event file
  file.Close();
  // Close output csv file
  csv_file.close();

  LOG("myAnalysis", pNOTICE)  << "Done!";

  return 0;
}
