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
#include "Physics/DeepInelastic/XSection/BYStrucFunc.h"
#include "Physics/Resonance/XSection/BergerSehgalRESPXSec2014.h"


using std::string;
using namespace genie;

// void GetCommandLineArgs (int argc, char ** argv);

// int start_ev = 500000;
// int end_ev = 1000000;
// string input_filename = "../data/gntp.3.ghep.root";

//___________________________________________________________________
int compute_xsec_genie(string input_filename, string output_eventfile, int start_ev, int end_ev) {
  
  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(input_filename.c_str(),"READ");

  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );

  if(!tree) return 1;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (end_ev > 0) ?
        TMath::Min(end_ev, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  // Open a csv file with event information, if the filename is specified
  std::ofstream event_file;
  // Write the header of the event file. Only do this if start_ev == 0
  if (start_ev == 0) {
    event_file.open(output_eventfile);
    event_file << "event_num,qel,res,dis,sea,hitqrk,xs,Q2s,hitnuc,atom,Mnuc,charm,resid,F1,F2,F3,F4,F5,sigmm,sigpp\n";
    event_file.close();
  }

  BYStrucFunc* structFunc = new BYStrucFunc();
  structFunc->Configure("Default");

  const DISStructureFuncModelI* fDISSFModel =
    dynamic_cast<const DISStructureFuncModelI *> (structFunc);
  
  DISStructureFunc* fDISSF = new DISStructureFunc();

  fDISSF->SetModel(fDISSFModel); // <-- attach algorithm
  
  BergerSehgalRESPXSec2014* resXSecCalculator = new BergerSehgalRESPXSec2014();
  resXSecCalculator->Configure("NoPauliBlock");
  
  //
  // Loop over all events
  //
  for(int i = start_ev; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);
    // Print the whole event:
    
    Interaction * in = event.Summary();
    const ProcessInfo & proc = in->ProcInfo();
    const XclsTag & xclsv = in->ExclTag();

    bool qelcc = proc.IsQuasiElastic() && proc.IsWeakCC();
    bool rescc = proc.IsResonant() && proc.IsWeakCC();
    bool discc = proc.IsDeepInelastic() && proc.IsWeakCC();
    
    const InitialState & init_state = in->InitState();
    const Target & tgt = init_state.Tgt();
    
    bool sea  = tgt.HitSeaQrk();
    int hitqrk = tgt.HitQrkPdg();

    Kinematics* kinematics = in->KinePtr();
    kinematics->UseSelectedKinematics();

    // Write event information to event file
    event_file.open(output_eventfile, std::ofstream::app);
    event_file << i << ",";
    // Set number of decimal places to 16
    event_file.precision(16);
    
    event_file << (qelcc ? "True" : "False") << "," << (rescc ? "True" : "False") << "," << (discc ? "True" : "False") << "," 
      << (sea ? "True" : "False") << "," << hitqrk << "," 
      << kinematics->x() << "," << kinematics->Q2() << "," 
      << tgt.HitNucPdg() << "," << tgt.A() << "," << tgt.HitNucP4().M() << "," << (xclsv.IsCharmEvent() ? "True" : "False") << "," << xclsv.Resonance();

    if (discc) {

      if (i < 10) {
        LOG("DISInt", pNOTICE) << *in << "\n";
      }

      fDISSF->Calculate(in);
      event_file << "," << fDISSF->F1() << "," << fDISSF->F2() << "," << fDISSF->F3() << "," << fDISSF->F4() << "," << fDISSF->F5() << ",,";
    }
    else if (rescc) {
      event_file << ",,,,,,";
      resXSecCalculator->XSec(in, kPSWQ2fE);
      std::vector<double> sigs = resXSecCalculator->GetSigs();
      event_file << sigs[0] << "," << sigs[1];
    }
    
    event_file << "\n";
    event_file.close();

    // clear current mc event record
    mcrec->Clear();

    // if (i % 100000 == 0) {
    //   LOG("myAnalysis", pNOTICE) << "Processed " << i << " events";
    // }
  }//end loop over events

  // Print total number of events
  LOG("myAnalysis", pNOTICE) << "Processed " << nev << " events";

  // close input GHEP event file
  file.Close();

  LOG("myAnalysis", pNOTICE)  << "Done!";

  return 0;
}
