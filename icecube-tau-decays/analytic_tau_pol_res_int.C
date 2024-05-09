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
#include "Physics/Resonance/XSection/BergerSehgalRESPXSec2014.h"


using namespace genie;


int analytic_tau_pol_res_int(std::string input_file, std::string output_file) {  
  // Open input file, and read the first line. 
  std::ifstream file(input_file);
  std::string line, word;
  std::getline(file, line);
  std::stringstream s(line);
  int index = 0;

  // TODO do a better solution for this using e.g. a map
  // Create a map with strings as keys, where the key is the column name, and the value is the index
  std::map<string, int> colname2index{
    {"xs", -1},
    {"ys", -1},
    {"Q2s", -1},
    {"Ws", -1},
    {"hitnuc", -1},
    {"res", -1},
    {"atom", -1},
    {"hitqrk", -1},
    {"sea", -1},
    {"En", -1},
    {"pxn", -1},
    {"pyn", -1},
    {"pzn", -1},
    {"Enu", -1},
    {"pxnu", -1},
    {"pynu", -1},
    {"pznu", -1},
    {"resid", -1}
  };
  
  // Identify the index where the value matches the column name
  while (std::getline(s, word, ',')) {
    // std::cout << word << std::endl;
    for (auto& kv : colname2index) {
      if (word == kv.first) {
        kv.second = index;
        std::cout << kv.first << ": " << kv.second << " ";
      }
    }
    index++;
  }
  std::cout << std::endl;
  
  // If any value is -1, then the column was not found
  for (auto& kv : colname2index) {
    if (kv.second == -1) {
      std::cerr << "Could not find '" << kv.first << "' in the header of the input file" << std::endl;
      return 1;
    }
  }
  
  // From https://github.com/GENIE-MC/Generator/blob/afed594388703bf4139fda97fbc3cfcdbe23e78a/src/contrib/corey/testXsec.C#L71
  // AlgConfigPool * pool = AlgConfigPool::Instance();
  // const Registry * config = pool->FindRegistry("genie::BergerSehgalRESPXSec2014/Default");

  // our xsec calculator
  BergerSehgalRESPXSec2014* xsecCalculator = new BergerSehgalRESPXSec2014();
  // Might be better to use Pauli blocking. Probably no difference because this happens after the sigpp, sigmm calculations
  // TODO try turning it on and off and compare outputs
  xsecCalculator->Configure("NoPauliBlock");

  //   const BergerSehgalRESPXSec2014* xsecCalculator = ;
  //   xsecCalculator->Configure("Default");
  
  // Open output file and write header to it
  std::ofstream ofile(output_file);
  ofile << line << ",sigmm,sigpp" << std::endl;
  
  double x, y, Q2val, W, En, pxn, pyn, pzn, Enu, pxnu, pynu, pznu;
  int nuc_pdgc, A, hitqrk;
  Resonance_t resid;
  bool res, sea;

  int event_num = 0;
  // Read the rest of the lines and compute the PDF at each (x, Q2)
  while (std::getline(file, line)) {
    std::stringstream s(line);
    index = 0;

    while (std::getline(s, word, ',')) {
      // if (event_num == 21) {
      //   LOG("myAnalysis", pINFO) << word;
      // }
      if (index == colname2index["res"]) {
        // Set res to True or False
        res = word == "True";
        if (!res) {
          break;
        }
      }

      if (index == colname2index["xs"]) {
        x = std::stod(word);
      } else if (index == colname2index["ys"]) {
          y = std::stod(word);
      } else if (index == colname2index["Q2s"]) {
          Q2val = std::stod(word);
      } else if (index == colname2index["Ws"]) {
          W = std::stod(word);
      } else if (index == colname2index["hitnuc"]) {
          // LOG("myAnalysis", pINFO) << "hitnuc" << word;
          nuc_pdgc = std::stoi(word);
      } else if (index == colname2index["atom"]) {
          // LOG("myAnalysis", pINFO) << "atom" << word;
          A = std::stoi(word);
      } else if (index == colname2index["hitqrk"]) {
          // LOG("myAnalysis", pINFO) << "hitqrk" << word;
          hitqrk = std::stoi(word);
      } else if (index == colname2index["sea"]) {
          sea = word == "True";
      } else if (index == colname2index["En"]) {
          En = std::stod(word);
      } else if (index == colname2index["pxn"]) {
          pxn = std::stod(word);
      } else if (index == colname2index["pyn"]) {
          pyn = std::stod(word);
      } else if (index == colname2index["pzn"]) {
          pzn = std::stod(word);
      } else if (index == colname2index["Enu"]) {
          Enu = std::stod(word);
      } else if (index == colname2index["pxnu"]) {
          pxnu = std::stod(word);
      } else if (index == colname2index["pynu"]) {
          pynu = std::stod(word);
      } else if (index == colname2index["pznu"]) {
          pznu = std::stod(word);
      } else if (index == colname2index["resid"]) {
          resid = static_cast<Resonance_t>(std::stoi(word));
      }
      index++;
    }
    if (!res) {
      ofile << line << ",,\n";
      continue;
    }
    // TODO combine this with analytic_tau_pol_dis_int.cc
    
    // Create ProcInfo object
    ProcessInfo proc_info(kScResonant, kIntWeakCC);
    if (!proc_info.IsResonant() || !proc_info.IsWeakCC()) {
      LOG("myAnalysis", pERROR) << "Not a RES CC event";
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
    
    // Set ExclTag
    XclsTag excl_tag;
    excl_tag.SetResonance(resid);
    interaction.SetExclTag(excl_tag);

    // Set the kinematics
    Kinematics* kine = interaction.KinePtr();
    // Maybe I only need to set 3 of these
    kine->Setx(x, true);
    kine->Sety(y, true);
    kine->SetQ2(Q2val, true);
    kine->SetW(W, true);
    
    kine->UseSelectedKinematics();
    // LOG("myAnalysis", pINFO) << event_num;
    if (event_num < 10) {
      LOG("myAnalysis", pINFO) << interaction;
    }

    // Compute cross section. 
    // NOTE The second variable changes which coordinates to calculate the differential cross section over. 
    // I don't think this matters to me, since I calculate the polarization for every particle individually.
    xsecCalculator->XSec(&interaction, kPSWQ2fE);
    std::vector<double> sigs = xsecCalculator->GetSigs();

    // Write to output file
    ofile << line << "," << sigs[0] << "," << sigs[1] << std::endl;
    ++event_num;
  }

  ofile.close();
  file.close();
  
  return 0;
}

