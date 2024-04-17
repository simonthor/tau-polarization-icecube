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
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/PartonDistributions/GRV98LO.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/PhysUtils.h"

using namespace genie;

// CKM mixing matrix constants
const double fVud2 = TMath::Power(0.97417, 2.);
const double fVus2 = TMath::Power(0.2248, 2.);
const double fVcd2 = TMath::Power(0.220, 2.);
const double fVcs2 = TMath::Power(0.995, 2.);
// include an R (~FL) factor into the calculation
const bool fIncludeR = true;
// include a nuclear factor (accounting for shadowing/anti-shadowing)
const bool fIncludeNuclMod = true;
// include corrections for calculating relation between 2xF1 and F2
const bool fUse2016Corrections = true;
// Select value for Q2 cutoff in relation between 2xF1 and F2
const double fLowQ2CutoffF1F2 = 0.8;
// BY constants
const double fA = 0.538;
const double fB = 0.305;
const double fCsU = 0.363;
const double fCsD = 0.621;
const double fCv1U = 0.291;
const double fCv2U = 0.189;
const double fCv1D = 0.202;
const double fCv2D = 0.255;
// Minimum Q2 value for PDF calculation
// Source: https://github.com/GENIE-MC/Generator/blob/1817b89cf815c0694c187f874288f1c1be1e712a/config/GRV98LO.xml#L11
const double fQ2min = 0.800;
// Charm mass, in GeV, according to Wikipedia
const double fMc = 1.27; 


double ScalingVar(double x, double Q2, const double fA, const double fB) {
// The modified Björken x that is used for BY
  double a  = TMath::Power( 2*constants::kProtonMass*x, 2 ) / Q2;
  double xw =  2*x*(Q2+fB) / (Q2*(1.+TMath::Sqrt(1+a)) +  2*fA*x);
  return xw;
}


double NuclMod(const bool fIncludeNuclMod, double x, int A)
{
// Nuclear modification to Fi

  double f = 1.;
  
  if(fIncludeNuclMod) {
    //   The x used for computing the DIS Nuclear correction factor should be the
    //   experimental x, not the rescaled x or off-shell-rest-frame version of x
    //   (i.e. selected x).  Since we do not have access to experimental x at this
    //   point in the calculation, just use selected x.
    f = utils::nuclear::DISNuclFactor(x, A);
  }

  return f;
}


double R(const bool fIncludeR, double x, double Q2) {
// Computes R ( ~ longitudinal structure function FL = R * 2xF1)
// The scaling variable can be overwritten to include corrections
//   The x used for computing the DIS Nuclear correction factor should be the
//   experimental x, not the rescaled x or off-shell-rest-frame version of x
//   (i.e. selected x).  Since we do not have access to experimental x at this
//   point in the calculation, just use selected x.
  if(fIncludeR) {
    double Rval = utils::phys::RWhitlow(x, Q2);
    return Rval;
  }
  return 0;
}


void KFactors(double Q2, double & kuv, double & kdv, double & kus, double & kds) {
// Compute the BY K factors for u(valence), d(valence), u(sea), d(sea)
// Note: it sets the values in-place

  double GD  = 1. / TMath::Power(1.+Q2/0.71, 2); // p elastic form factor
  double GD2 = TMath::Power(GD,2);

  kuv = (1.-GD2)*(Q2+fCv2U)/(Q2+fCv1U); // K - u(valence)
  kdv = (1.-GD2)*(Q2+fCv2D)/(Q2+fCv1D); // K - d(valence)
  kus = Q2/(Q2+fCsU);                   // K - u(sea)
  kds = Q2/(Q2+fCsD);                   // K - d(sea)
}


std::vector<double> CalculatePDFs(double x, double Q2val, double M, int nuc_pdgc, PDF* fPDF, PDF* fPDFc) {
  // Clean-up previous calculation
  fPDF->Reset();
  fPDFc->Reset();
  double Q2pdf = TMath::Max(Q2val, fQ2min);
  // Compute PDFs at (x,Q2)
  fPDF->Calculate(x, Q2pdf);

  bool above_charm =
      utils::kinematics::IsAboveCharmThreshold(x, Q2val, M, fMc);
  if(above_charm) {
    // compute the slow rescaling var
    double xc = utils::kinematics::SlowRescalingVar(x, Q2val, M, fMc);
    if(xc<0 || xc>1) {
      std::cerr << "Unphysical slow rescaling var: xc = " << xc;
      throw std::invalid_argument("Unphysical slow rescaling var" + std::to_string(xc));
    } else {
      // compute PDFs at (xc,Q2)
      fPDFc->Calculate(xc, Q2pdf);
    }
  }
  // Compute the K factors
  double kval_u = 1.;
  double kval_d = 1.;
  double ksea_u = 1.;
  double ksea_d = 1.;

  KFactors(Q2val, kval_u, kval_d, ksea_u, ksea_d);

  // std::cout << "K factors: " << kval_u << " " << kval_d << " " << ksea_u << " " << ksea_d << std::endl;
  // Apply the K factors
  //
  // Always scale d pdfs with d kfactors and u pdfs with u kfactors.
  // Don't swap the applied kfactors for neutrons.
  // Debdatta & Donna noted (Sep.2006) that a similar swap in the neugen
  // implementation was the cause of the difference in nu and nubar F2
  //
  fPDF->ScaleUpValence   (kval_u);
  fPDF->ScaleDownValence (kval_d);
  fPDF->ScaleUpSea       (ksea_u);
  fPDF->ScaleDownSea     (ksea_d);
  fPDF->ScaleStrange     (ksea_d);
  fPDF->ScaleCharm       (ksea_u);
  if(above_charm) {
    fPDFc->ScaleUpValence   (kval_u);
    fPDFc->ScaleDownValence (kval_d);
    fPDFc->ScaleUpSea       (ksea_u);
    fPDFc->ScaleDownSea     (ksea_d);
    fPDFc->ScaleStrange     (ksea_d);
    fPDFc->ScaleCharm       (ksea_u);
  }

  // Rules of thumb
  // ---------------------------------------
  // - For W+ exchange use: -1/3|e| quarks and -2/3|e| antiquarks
  // - For W- exchange use:  2/3|e| quarks and  1/3|e| antiquarks
  // - For each qi -> qj transition multiply with the (ij CKM element)^2
  // - Use isospin symmetry to get neutron's u,d from proton's u,d
  //    -- neutron d = proton u
  //    -- neutron u = proton d
  // - Use u = usea + uvalence. Same for d
  // - For s,c use q=qbar
  // - For t,b use q=qbar=0

  double fuv   = fPDF  -> UpValence();
  double fus   = fPDF  -> UpSea();
  double fdv   = fPDF  -> DownValence();
  double fds   = fPDF  -> DownSea();
  double fs    = fPDF  -> Strange();
  double fc    = 0.;
  
  // will be 0 if < charm threshold
  double fuv_c = fPDFc -> UpValence();   
  double fus_c = fPDFc -> UpSea();
  double fdv_c = fPDFc -> DownValence();
  double fds_c = fPDFc -> DownSea();    
  double fs_c  = fPDFc -> Strange();    
  double fc_c  = fPDFc -> Charm();      
  
  // if (above_charm) {
  //   // std::cout << "fuv: " << fuv << " fus: " << fus << " fdv: " << fdv << " fds: " << fds << " fs: " << fs << " fc: " << fc << std::endl;
  //   std::cout << "fuv_c: " << fuv_c << " fus_c: " << fus_c << " fdv_c: " << fdv_c << " fds_c: " << fds_c << " fs_c: " << fs_c << " fc_c: " << fc_c << std::endl;
  // }

  // The above are the proton parton density function. Get the PDFs for the
  // hit nucleon (p or n) by swapping u<->d if necessary
  bool isP = pdg::IsProton  (nuc_pdgc);
  bool isN = pdg::IsNeutron (nuc_pdgc);
  assert(isP || isN);

  double tmp = 0;
  if (isN) {  // swap u <-> d
    tmp = fuv;   fuv   = fdv;   fdv   = tmp;
    tmp = fus;   fus   = fds;   fds   = tmp;
    tmp = fuv_c; fuv_c = fdv_c; fdv_c = tmp;
    tmp = fus_c; fus_c = fds_c; fds_c = tmp;
  }

  return {fuv, fus, fdv, fds, fs, fc, fuv_c, fus_c, fdv_c, fds_c, fs_c, fc_c};
}


std::vector<double> structureFunctions(std::vector<double> pdf_values, double Q2val, double x, double bjx, int A, int hitqrk, bool sea) 
{
  // Compute the structure functions
  double fF1 = 0;
  double fF2 = 0;
  double fF3 = 0;
  double fF4 = 0;
  double fF5 = 0;
  double fF6 = 0;
  
  double fuv   = pdf_values[0];
  double fus   = pdf_values[1];
  double fdv   = pdf_values[2];
  double fds   = pdf_values[3];
  double fs    = pdf_values[4];
  double fc    = pdf_values[5];
  double fuv_c = pdf_values[6];
  double fus_c = pdf_values[7];
  double fdv_c = pdf_values[8];
  double fds_c = pdf_values[9];
  double fs_c  = pdf_values[10];
  double fc_c  = pdf_values[11];

  double switch_uv    = 1.;
  double switch_us    = 1.;
  double switch_ubar  = 1.;
  double switch_dv    = 1.;
  double switch_ds    = 1.;
  double switch_dbar  = 1.;
  double switch_s     = 1.;
  double switch_sbar  = 1.;
  double switch_c     = 1.;
  double switch_cbar  = 1.;

  // TODO read this from the file (file name, separate input parameter, line or something else) later
  bool is_nu       = true;
  bool is_nubar    = false;
  
  if(hitqrk != 0) {
    switch_uv    = 0.;
    switch_us    = 0.;
    switch_ubar  = 0.;
    switch_dv    = 0.;
    switch_ds    = 0.;
    switch_dbar  = 0.;
    switch_s     = 0.;
    switch_sbar  = 0.;
    switch_c     = 0.;
    switch_cbar  = 0.;

    int  qpdg = hitqrk;

    bool is_u    = pdg::IsUQuark     (qpdg);
    bool is_ubar = pdg::IsAntiUQuark (qpdg);
    bool is_d    = pdg::IsDQuark     (qpdg);
    bool is_dbar = pdg::IsAntiDQuark (qpdg);
    bool is_s    = pdg::IsSQuark     (qpdg);
    bool is_sbar = pdg::IsAntiSQuark (qpdg);
    bool is_c    = pdg::IsCQuark     (qpdg);
    bool is_cbar = pdg::IsAntiCQuark (qpdg);

    if      (!sea && is_u   ) { switch_uv   = 1; }
    else if ( sea && is_u   ) { switch_us   = 1; }
    else if ( sea && is_ubar) { switch_ubar = 1; }
    else if (!sea && is_d   ) { switch_dv   = 1; }
    else if ( sea && is_d   ) { switch_ds   = 1; }
    else if ( sea && is_dbar) { switch_dbar = 1; }
    else if ( sea && is_s   ) { switch_s    = 1; }
    else if ( sea && is_sbar) { switch_sbar = 1; }
    else if ( sea && is_c   ) { switch_c    = 1; }
    else if ( sea && is_cbar) { switch_cbar = 1; }
    else throw std::invalid_argument("Unknown quark type");

     // make sure user inputs make sense
    if(is_nu    && is_u   ) throw std::invalid_argument("is_nu && is_u");
    if(is_nu    && is_c   ) throw std::invalid_argument("is_nu && is_c");
    if(is_nu    && is_dbar) throw std::invalid_argument("is_nu && is_dbar");
    if(is_nu    && is_sbar) throw std::invalid_argument("is_nu && is_sbar");
    if(is_nubar && is_ubar) throw std::invalid_argument("is_nubar && is_ubar");
    if(is_nubar && is_cbar) throw std::invalid_argument("is_nubar && is_cbar");
    if(is_nubar && is_d   ) throw std::invalid_argument("is_nubar && is_d");
    if(is_nubar && is_s   ) throw std::invalid_argument("is_nubar && is_s");
  }
  double F2val=0, xF3val=0;

  double q=0, qbar=0;

  if (is_nu) {
    q    = ( switch_dv * fdv   + switch_ds * fds   ) * fVud2 +
            ( switch_s  * fs                        ) * fVus2 +
            ( switch_dv * fdv_c + switch_ds * fds_c ) * fVcd2 +
            ( switch_s  * fs_c                      ) * fVcs2;

    qbar = ( switch_ubar * fus  ) * fVud2 +
            ( switch_ubar * fus  ) * fVus2 +
            ( switch_cbar * fc_c ) * fVcd2 +
            ( switch_cbar * fc_c ) * fVcs2;
  }
  else
  if (is_nubar) {
    q    = ( switch_uv * fuv + switch_us * fus    ) * fVud2 +
            ( switch_uv * fuv + switch_us * fus    ) * fVus2 +
            ( switch_c  * fc_c                     ) * fVcd2 +
            ( switch_c  * fc_c                     ) * fVcs2;

    qbar = ( switch_dbar * fds_c ) * fVcd2 +
            ( switch_dbar * fds   ) * fVud2 +
            ( switch_sbar * fs    ) * fVus2 +
            ( switch_sbar * fs_c  ) * fVcs2;
  }
  else {
    std::cerr << "ERROR: both is_nu and is_nubar are false" << std::endl;
    throw std::invalid_argument("Both is_nu and is_nubar are false");
  }

  F2val  = 2*(q+qbar);
  xF3val = 2*(q-qbar);

  // Nuclear modification coefficient. Note that the regular Björken x is used here
  double f = NuclMod(fIncludeNuclMod, bjx, A);
  
  // R ~ FL
  double r = R(fIncludeR, bjx, Q2val); 


  if(fUse2016Corrections) {
    //It was confirmed by A.Bodek that the modified scaling variable
    //should just be used to compute the strucure functions F2 and xF3,
    //but that the usual Bjorken x should be used for the relations
    //between the structure functions.
    //For the same reason remove the freezing of Q2 at 0.8 for those relations,
    //although it has not been explicitly asked to A.Bodek if it should be done.

    double a = TMath::Power(bjx,2.) / TMath::Max(Q2val, fLowQ2CutoffF1F2);
    // kNucleonMass2 is a constant from GENIE
    double c = (1. + 4. * constants::kNucleonMass2 * a) / (1.+r);

    fF3 = f * xF3val/bjx;
    fF2 = f * F2val;
    fF1 = fF2 * 0.5*c/bjx;
    fF5 = fF2/bjx;           // Albright-Jarlskog relation
    fF4 = 0.;                // Nucl.Phys.B 84, 467 (1975)
  }
  else {
    double a = TMath::Power(x,2.) / TMath::Max(Q2val, fLowQ2CutoffF1F2);
    double c = (1. + 4. * constants::kNucleonMass2 * a) / (1.+r);

    fF3 = f * xF3val / x;
    fF2 = f * F2val;
    fF1 = fF2 * 0.5 * c / x;
    fF5 = fF2 / x;         // Albright-Jarlskog relation
    fF4 = 0.;              // Nucl.Phys.B 84, 467 (1975)
  }

  return {fF1, fF2, fF3, fF4, fF5};
}


int analytic_tau_pol_dis(std::string input_file, std::string output_file) {
  // PDF model
  const PDFModelI * pdf_model =
        dynamic_cast<const PDFModelI *> (new GRV98LO());
  PDF* fPDF =  new PDF();
  // PDF used when above charm production threshold
  PDF* fPDFc = new PDF();
  
  fPDF->SetModel(pdf_model);
  fPDFc->SetModel(pdf_model);
  
  // Open input file, and read the first line. 
  std::ifstream file(input_file);
  std::string line, word;
  std::getline(file, line);
  std::stringstream s(line);
  int index = 0;

  int xcol, Q2col, nuc_pdgcol, discol, Acol, Mcol, hitqrkcol, seacol;
  xcol = Q2col = nuc_pdgcol = discol = Acol = Mcol = hitqrkcol = seacol = 0;

  // Identify the index where the value is "x" and "Q2"
  while (std::getline(s, word, ',')) {
    // std::cout << word << std::endl;
    if (word == "xs") {
        xcol = index;
    } else if (word == "Q2s") {
        Q2col = index;
    } else if (word == "hitnuc") {
        nuc_pdgcol = index;
    } else if (word == "dis") {
        discol = index;
    } else if (word == "atom") {
        Acol = index;
    } else if (word == "Mnuc") {
        Mcol = index;
    } else if (word == "hitqrk") {
        hitqrkcol = index;
    } else if (word == "sea") {
        seacol = index;
    }
    index++;
  }

  std::cout << "xcol: " << xcol << ", Q2col: " << Q2col 
    << ", nuc_pdgcol: " << nuc_pdgcol << ", discol: " << discol 
    << ", Acol: " << Acol << ", Mcol: " << Mcol << ", hitqrkcol: " << hitqrkcol
    << ", seacol: " << seacol
    << std::endl;
  
  if ((xcol == 0) && (Q2col == 0)) {
    std::cerr << "Could not find 'xs' and 'Q2s' in the header of the input file" << std::endl;
    return 1;
  }
  
  // Open output file and write header to it
  std::ofstream ofile(output_file);
  ofile << line << ",F1,F2,F3,F4,F5" << std::endl;
  
  double x, Q2val, M;
  int nuc_pdgc, A, hitqrk;
  bool dis, sea;
  // Read the rest of the lines and compute the PDF at each (x, Q2)
  while (std::getline(file, line)) {
    std::stringstream s(line);
    index = 0;
    while (std::getline(s, word, ',')) {
      if (index == xcol) {
          x = std::stod(word);
      } else if (index == Q2col) {
          Q2val = std::stod(word);
      } else if (index == nuc_pdgcol) {
          nuc_pdgc = std::stoi(word);
      } else if (index == discol) {
          // Set dis to True or False
          dis = word == "True";
      } else if (index == Acol) {
          A = std::stoi(word);
      } else if (index == Mcol) {
          M = std::stod(word);
      } else if (index == hitqrkcol) {
          hitqrk = std::stoi(word);
      } else if (index == seacol) {
          sea = word == "True";
      }
      index++;
    }
    
    if (!dis) {
      ofile << line << ",,,,,\n";
      continue;
    }

    // Save regular Björken x in a variable, as a modified x will be used for many other calculations
    double bjx = x;

    // Modify x to use the BY scaling variable
    x = ScalingVar(x, Q2val, fA, fB);
    // Q2 is left as it is
    // std::cout << "Björken x: " << bjx << ", Modified x: " << x << ", Q2: " << Q2val << std::endl;
    // Compute PDFs
    std::vector<double> pdf_values = CalculatePDFs(x, Q2val, M, nuc_pdgc, fPDF, fPDFc);
    // Compute structure functions
    std::vector<double> Fs = structureFunctions(pdf_values, Q2val, x, bjx, A, hitqrk, sea);
    
    // Write to output file
    ofile << line << "," << Fs[0] << "," << Fs[1] << "," << Fs[2] << "," << Fs[3] << "," << Fs[4] << std::endl;
  }

  delete fPDF;
  delete fPDFc;

  ofile.close();
  file.close();
  
  return 0;
}

