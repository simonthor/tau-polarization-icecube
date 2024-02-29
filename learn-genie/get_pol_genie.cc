#include <string>

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

int    gOptNEvt = -1;
string gOptInpFilename = "gntp.0.ghep.root";

//___________________________________________________________________
int get_pol_genie()
{
  // GetCommandLineArgs (argc, argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  int good_tau_events = 0;
  int too_many_tau_events = 0;
  int unpolarized_tau_events = 0;
  int no_tau_events = 0;

  //
  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);

    // LOG("myAnalysis", pNOTICE) << event;

    //
    // Put your event analysis code here 
    //
    // ... ... ... ... ...
    // ... ... ... ... ...
    //
    // 

    

    //
    // Loop over all particles in this event
    //

    GHepParticle * p = 0;
    TIter event_iter(&event);

    // Number of taus in an event
    int n_taus = 0;

    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
       //
       // Put your event analysis code here 
       //
       // ... ... ... ... ...
       // ... ... ... ... ...
       //
       //

       // EXAMPLE: Print out the energy of all final state pions.
      if ((p->Status() == kIStStableFinalState) && (p->Pdg() == kPdgTau)) {
        LOG("myAnalysis", pNOTICE)  
          << "Got a : " << p->Name() << " with E = " << p->E() << " GeV"
          << " and polarization: \n" 
          << " polar angle " << p->PolzPolarAngle()
          << " azimuth angle " << p->PolzAzimuthAngle()
          << " polarization is set: " << p->PolzIsSet() << "\n";
          //<< " polarization vector: " << p->GetPolarization(); //.X() << " " << p->GetPolarization().Y() << " " << p->GetPolarization().Z() << " " << p->GetPolarization().T();
        n_taus++;
        if (p->PolzIsSet()) {
          good_tau_events++;
        }
        else {
          unpolarized_tau_events++;
        }
      }

    }// end loop over particles	
    if (n_taus == 0) {
      LOG("myAnalysis", pNOTICE) << "No taus in this event";
      no_tau_events++;
    }
    else if (n_taus > 1) {
      LOG("myAnalysis", pNOTICE) << "Number of taus in this event: " << n_taus;
      too_many_tau_events++;
    }

    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  LOG("myAnalysis", pNOTICE) << "Number of good tau events: " << good_tau_events;
  LOG("myAnalysis", pNOTICE) << "Number of unpolarized tau events: " << unpolarized_tau_events;
  LOG("myAnalysis", pNOTICE) << "Number of no tau events: " << no_tau_events;
  LOG("myAnalysis", pNOTICE) << "Number of too many tau events: " << too_many_tau_events;
  
  // close input GHEP event file
  file.Close();

  LOG("myAnalysis", pNOTICE)  << "Done!";

  return 0;
}
