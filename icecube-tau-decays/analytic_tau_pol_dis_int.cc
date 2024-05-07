// #include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
// #include "Physics/PartonDistributions/GRV98LO.h"

// #include "Framework/Algorithm/AlgConfigPool.h"
// #include "Framework/Conventions/GBuild.h"
// #include "Framework/Conventions/RefFrame.h"
// #include "Framework/Messenger/Messenger.h"
// #include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
#include <TMath.h>
#include <iostream>
#include <stdexcept>

#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/DeepInelastic/XSection/BYStrucFunc.h"
#include "Physics/DeepInelastic/XSection/DISStructureFunc.h"

using namespace genie;


int analytic_tau_pol_dis_int(std::string input_file, std::string output_file) {  
  // Open input file, and read the first line. 
  std::ifstream file(input_file);
  std::string line, word;
  std::getline(file, line);
  std::stringstream s(line);
  int index = 0;

  // TODO do a better solution for this using e.g. a map
  int xcol, ycol, Q2col, Wcol, nuc_pdgcol, discol, Acol, hitqrkcol, seacol, Encol, pxncol, pyncol, pzncol, Enucol, pxnucol, pynucol, pznucol;
  xcol = ycol = Q2col = Wcol = nuc_pdgcol = discol = Acol = hitqrkcol = seacol = Encol = pxncol = pyncol = pzncol = Enucol = pxnucol = pynucol = pznucol = 0;

  // Identify the index where the value matches the column name
  while (std::getline(s, word, ',')) {
    // std::cout << word << std::endl;
    if (word == "xs") {
        xcol = index;
    } else if (word == "ys") {
        ycol = index;
    } else if (word == "Q2s") {
        Q2col = index;
    } else if (word == "Ws") {
        Wcol = index;
    } else if (word == "hitnuc") {
        nuc_pdgcol = index;
    } else if (word == "dis") {
        discol = index;
    } else if (word == "atom") {
        Acol = index;
    } else if (word == "hitqrk") {
        hitqrkcol = index;
    } else if (word == "sea") {
        seacol = index;
    } else if (word == "En") {
        Encol = index;
    } else if (word == "pxn") {
        pxncol = index;
    } else if (word == "pyn") {
        pyncol = index;
    } else if (word == "pzn") {
        pzncol = index;
    } else if (word == "Enu") {
        Enucol = index;
    } else if (word == "pxnu") {
        pxnucol = index;
    } else if (word == "pynu") {
        pynucol = index;
    } else if (word == "pznu") {
        pznucol = index;
    }
    index++;
  }

  std::cout << "xcol: " << xcol << ", ycol" << ycol << ", Q2col: " << Q2col << ", Wcol: " << Wcol
    << ", nuc_pdgcol: " << nuc_pdgcol << ", discol: " << discol 
    << ", Acol: " << Acol << ", hitqrkcol: " << hitqrkcol
    << ", seacol: " << seacol
    << ", Encol: " << Encol << ", pxncol: " << pxncol << ", pyncol: " << pyncol << ", pzncol: " << pzncol
    << ", Enucol: " << Enucol << ", pxnucol: " << pxnucol << ", pynucol: " << pynucol << ", pznucol: " << pznucol
    << std::endl;
  
  if ((xcol == 0) && (Q2col == 0)) {
    std::cerr << "Could not find 'xs' and 'Q2s' in the header of the input file" << std::endl;
    return 1;
  }
  
  BYStrucFunc* structFunc = new BYStrucFunc();
  structFunc->Configure("Default");

  const DISStructureFuncModelI* fDISSFModel =
    dynamic_cast<const DISStructureFuncModelI *> (structFunc);
  
  DISStructureFunc* fDISSF = new DISStructureFunc();

  fDISSF->SetModel(fDISSFModel); // <-- attach algorithm
  
  // Open output file and write header to it
  std::ofstream ofile(output_file);
  ofile << line << ",F1,F2,F3,F4,F5" << std::endl;
  
  double x, y, Q2val, W, En, pxn, pyn, pzn, Enu, pxnu, pynu, pznu;
  int nuc_pdgc, A, hitqrk;
  bool dis, sea;
  int event_num = 0;
  // Read the rest of the lines and compute the PDF at each (x, Q2)
  while (std::getline(file, line)) {
    std::stringstream s(line);
    index = 0;
    while (std::getline(s, word, ',')) {
      if (index == xcol) {
          x = std::stod(word);
      } else if (index == ycol) {
          y = std::stod(word);
      } else if (index == Q2col) {
          Q2val = std::stod(word);
      } else if (index == Wcol) {
          W = std::stod(word);
      } else if (index == nuc_pdgcol) {
          nuc_pdgc = std::stoi(word);
      } else if (index == discol) {
          // Set dis to True or False
          dis = word == "True";
      } else if (index == Acol) {
          A = std::stoi(word);
      } else if (index == hitqrkcol) {
          hitqrk = std::stoi(word);
      } else if (index == seacol) {
          sea = word == "True";
      } else if (index == Encol) {
          En = std::stod(word);
      } else if (index == pxncol) {
          pxn = std::stod(word);
      } else if (index == pyncol) {
          pyn = std::stod(word);
      } else if (index == pzncol) {
          pzn = std::stod(word);
      } else if (index == Enucol) {
          Enu = std::stod(word);
      } else if (index == pxnucol) {
          pxnu = std::stod(word);
      } else if (index == pynucol) {
          pynu = std::stod(word);
      } else if (index == pznucol) {
          pznu = std::stod(word);
      }
      index++;
    }
    
    if (!dis) {
      ofile << line << ",,,,,\n";
      continue;
    }

    // Create ProcInfo object
    ProcessInfo proc_info(kScDeepInelastic, kIntWeakCC);
    if (!proc_info.IsDeepInelastic() || !proc_info.IsWeakCC()) {
      LOG("myAnalysis", pERROR) << "Not a DIS CC event";
    }

    // Create InitialState object
    // TODO read this from file instead
    int Z = (A == 16) ? 8 : 1;
    int nu_pdg = 16;
    Target tgt(Z, A);

    InitialState init_state(tgt, nu_pdg);
    
    TLorentzVector nu4m(pxnu, pynu, pznu, Enu);

    init_state.SetProbeP4(nu4m);

    Interaction interaction(init_state, proc_info);

    Target* tgt_ptr = interaction.InitStatePtr()->TgtPtr();
    tgt_ptr->SetHitNucPdg(nuc_pdgc);
    TLorentzVector nucleon4m(pxn, pyn, pzn, En);
    tgt_ptr->SetHitNucP4(nucleon4m);
    tgt_ptr->SetHitQrkPdg(hitqrk);
    tgt_ptr->SetHitSeaQrk(sea);

    Kinematics* kine = interaction.KinePtr();
    // Maybe I only need to set 3 of these
    kine->Setx(x, true);
    kine->Sety(y, true);
    kine->SetQ2(Q2val, true);
    kine->SetW(W, true);
    
    kine->UseSelectedKinematics();
    
    if (event_num < 10) {
      LOG("myAnalysis", pINFO) << interaction;
    }

    // Compute structure functions
    fDISSF->Calculate(&interaction);

    // Write to output file
    ofile << line << "," << fDISSF->F1() << "," << fDISSF->F2() << "," << fDISSF->F3() << "," << fDISSF->F4() << "," << fDISSF->F5() << std::endl;
    ++event_num;
  }

  ofile.close();
  file.close();
  
  return 0;
}

