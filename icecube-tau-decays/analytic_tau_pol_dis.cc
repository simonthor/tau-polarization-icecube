// #include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
// #include "Physics/PartonDistributions/GRV98LO.h"

// #include "Framework/Algorithm/AlgConfigPool.h"
// #include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
// #include "Framework/Conventions/RefFrame.h"
// #include "Framework/Messenger/Messenger.h"
// #include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
#include <TMath.h>
#include <iostream>
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/PartonDistributions/GRV98LO.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
// #include "Physics/NuclearState/NuclearUtils.h"
// #include "Framework/Utils/PhysUtils.h"

using namespace genie;


int analytic_tau_pol_dis() {
    const PDFModelI * pdf_model =
         dynamic_cast<const PDFModelI *> (new GRV98LO());
    PDF* fPDF =  new PDF();
    // PDF used when above charm production threshold
    PDF* fPDFc = new PDF();
    
    fPDF->SetModel(pdf_model);
    fPDFc->SetModel(pdf_model);
    
    fPDF->Reset();
    fPDFc->Reset();

    // Source: https://github.com/GENIE-MC/Generator/blob/1817b89cf815c0694c187f874288f1c1be1e712a/config/GRV98LO.xml#L11
    const double fQ2min = 0.800;
    // Charm mass, in GeV, according to Wikipedia
    const double fMc = 1.27; 

    // Variables to be read from a file
    double x = 0.1;
    double Q2 = 0.7; // Q2 in GeV^2
    double M = 0.938; // Hit nucleon mass, in GeV
    double nuc_pdgc = 2212; // Hit nucleon PDG code. 2212 is a proton, 2112 is a neutron

    double Q2pdf =  TMath::Max(Q2, fQ2min);

    // Compute PDFs at (x,Q2)
    fPDF->Calculate(x, Q2pdf);

    bool above_charm = utils::kinematics::IsAboveCharmThreshold(x, Q2, M, fMc);
    if(above_charm) {
        // compute the slow rescaling var
        double xc = utils::kinematics::SlowRescalingVar(x, Q2, M, fMc);
        if(xc<0 || xc>1) {
            std::cout << "Unphysical slow rescaling var: xc = " << xc;
            // return 1;
        } else {
            // compute PDFs at (xc,Q2)
            fPDFc->Calculate(xc, Q2pdf);
        }
    }

    double fuv =  fPDF  -> UpValence();
    double fus   = fPDF  -> UpSea();
    double fdv   = fPDF  -> DownValence();
    double fds   = fPDF  -> DownSea();
    double fs    = fPDF  -> Strange();
    double fc    = 0.;
    std::cout << "fuv: " << fuv << " fus: " << fus << " fdv: " << fdv << " fds: " << fds << " fs: " << fs << " fc: " << fc << std::endl;
    
    // will be 0 if < charm threshold
    double fuv_c = fPDFc -> UpValence();   
    double fus_c = fPDFc -> UpSea();
    double fdv_c = fPDFc -> DownValence();
    double fds_c = fPDFc -> DownSea();    
    double fs_c  = fPDFc -> Strange();    
    double fc_c  = fPDFc -> Charm();      
    std::cout << "fuv_c: " << fuv_c << " fus_c: " << fus_c << " fdv_c: " << fdv_c << " fds_c: " << fds_c << " fs_c: " << fs_c << " fc_c: " << fc_c << std::endl;

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

    return 0;
}


void QPMDISStrucFuncBase::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;

  // Get process info & perform various checks
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();

  int  nuc_pdgc    = tgt.HitNucPdg();
  int  probe_pdgc  = init_state.ProbePdg();
  bool is_p        = pdg::IsProton       ( nuc_pdgc    );
  bool is_n        = pdg::IsNeutron      ( nuc_pdgc    );
  bool is_nu       = pdg::IsNeutrino     ( probe_pdgc  );
  bool is_nubar    = pdg::IsAntiNeutrino ( probe_pdgc  );
  bool is_lepton   = pdg::IsLepton       ( probe_pdgc  );
  bool is_dm       = pdg::IsDarkMatter   ( probe_pdgc  );
  bool is_CC       = proc_info.IsWeakCC();
  bool is_NC       = proc_info.IsWeakNC();
  bool is_EM       = proc_info.IsEM();
  bool is_dmi      = proc_info.IsDarkMatter();

  if ( !is_lepton && !is_dm ) return;
  if ( !is_p && !is_n       ) return;
  if ( tgt.N() == 0 && is_n ) return;
  if ( tgt.Z() == 0 && is_p ) return;

  // Flags switching on/off quark contributions so that this algorithm can be
  // used for both l + N -> l' + X, and l + q -> l' + q' level calculations
  
  // TODO maybe use this for my own calculations, if quark level interactions are happening. For now, ignore.
//   double switch_uv    = 1.;
//   double switch_us    = 1.;
//   double switch_ubar  = 1.;
//   double switch_dv    = 1.;
//   double switch_ds    = 1.;
//   double switch_dbar  = 1.;
//   double switch_s     = 1.;
//   double switch_sbar  = 1.;
//   double switch_c     = 1.;
//   double switch_cbar  = 1.;

//   if(tgt.HitQrkIsSet()) {

//      switch_uv    = 0.;
//      switch_us    = 0.;
//      switch_ubar  = 0.;
//      switch_dv    = 0.;
//      switch_ds    = 0.;
//      switch_dbar  = 0.;
//      switch_s     = 0.;
//      switch_sbar  = 0.;
//      switch_c     = 0.;
//      switch_cbar  = 0.;

//      int  qpdg = tgt.HitQrkPdg();
//      bool sea  = tgt.HitSeaQrk();

//      bool is_u    = pdg::IsUQuark     (qpdg);
//      bool is_ubar = pdg::IsAntiUQuark (qpdg);
//      bool is_d    = pdg::IsDQuark     (qpdg);
//      bool is_dbar = pdg::IsAntiDQuark (qpdg);
//      bool is_s    = pdg::IsSQuark     (qpdg);
//      bool is_sbar = pdg::IsAntiSQuark (qpdg);
//      bool is_c    = pdg::IsCQuark     (qpdg);
//      bool is_cbar = pdg::IsAntiCQuark (qpdg);

//      if      (!sea && is_u   ) { switch_uv   = 1; }
//      else if ( sea && is_u   ) { switch_us   = 1; }
//      else if ( sea && is_ubar) { switch_ubar = 1; }
//      else if (!sea && is_d   ) { switch_dv   = 1; }
//      else if ( sea && is_d   ) { switch_ds   = 1; }
//      else if ( sea && is_dbar) { switch_dbar = 1; }
//      else if ( sea && is_s   ) { switch_s    = 1; }
//      else if ( sea && is_sbar) { switch_sbar = 1; }
//      else if ( sea && is_c   ) { switch_c    = 1; }
//      else if ( sea && is_cbar) { switch_cbar = 1; }
//      else return;

//      // make sure user inputs make sense
//     if(is_nu    && is_CC && is_u   ) return;
//     if(is_nu    && is_CC && is_c   ) return;
//     if(is_nu    && is_CC && is_dbar) return;
//     if(is_nu    && is_CC && is_sbar) return;
//     if(is_nubar && is_CC && is_ubar) return;
//     if(is_nubar && is_CC && is_cbar) return;
//     if(is_nubar && is_CC && is_d   ) return;
//     if(is_nubar && is_CC && is_s   ) return;
//   }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Applying all PDF K-factors abd scaling variable corrections

  this -> CalcPDFs (interaction);

  //
  // Compute structure functions for the EM, NC and CC cases
  //

  double F2val=0, xF3val=0;
  
  // NOTE for me, only CC is relevant
  // ***  CHARGED CURRENT

  if(is_CC) {
    double q=0, qbar=0;

    if (is_nu) {
      q    = ( switch_dv * fdv   + switch_ds * fds   ) * fVud2 + // CKM mixing matrix value ^2
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
      return;
    }

    F2val  = 2*(q+qbar);
    xF3val = 2*(q-qbar);
  }

  double Q2val = this->Q2        (interaction);
  double x     = this->ScalingVar(interaction);

  double f = 1.;

  if (fIncludeNuclMod) {
    int A = 18; // Oxygen. 
    // TODO change this to be a parameter read from a file
    f = utils::nuclear::DISNuclFactor(x,A);
  }

  // longitudinal structure function FL = R * 2xF1
  double r = 0;
  if (fIncludeR) {
     double r = utils::phys::RWhitlow(x, Q2val);
  }


  //It was confirmed by A.Bodek that the modified scaling variable
  //should just be used to compute the strucure functions F2 and xF3,
  //but that the usual Bjorken x should be used for the relations
  //between the structure functions.
  //For the same reason remove the freezing of Q2 at 0.8 for those relations,
  //although it has not been explicitly asked to A.Bodek if it should be done.
  
//   const Kinematics & kinematics = interaction->Kine();
//   double bjx = kinematics.x();
  
  double bjx = x;
  double a = TMath::Power(bjx,2.) / TMath::Max(Q2, 0.8);
  // I am guessing that kNucleonMass2 comes from Constants.h
  double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);
  
  fF3 = f * xF3val/bjx;
  fF2 = f * F2val;
  fF1 = fF2 * 0.5*c/bjx;
  fF5 = fF2/bjx;           // Albright-Jarlskog relation
  fF4 = 0.;                // Nucl.Phys.B 84, 467 (1975)
}


// NOTE check how large these corrections are for my purposes
double QPMDISStrucFuncBase::NuclMod(const Interaction * interaction) const
{
// Nuclear modification to Fi
// The scaling variable can be overwritten to include corrections

  if( interaction->TestBit(kIAssumeFreeNucleon)   ) return 1.0;
  if( interaction->TestBit(kINoNuclearCorrection) ) return 1.0;

  double f = 1.;
  if(fIncludeNuclMod) {
     const Target & tgt  = interaction->InitState().Tgt();

//   The x used for computing the DIS Nuclear correction factor should be the
//   experimental x, not the rescaled x or off-shell-rest-frame version of x
//   (i.e. selected x).  Since we do not have access to experimental x at this
//   point in the calculation, just use selected x.
     const Kinematics & kine  = interaction->Kine();
     double x  = kine.x();
     int    A = tgt.A();
     f = utils::nuclear::DISNuclFactor(x,A);
  }

  return f;
}


double BYStrucFunc::ScalingVar(const Interaction * interaction) const
{
// Overrides QPMDISStrucFuncBase::ScalingVar() to compute the BY scaling var

  const Kinematics & kine  = interaction->Kine();
  double x  = kine.x();
  double myQ2 = this->Q2(interaction);
  //myQ2 = TMath::Max(Q2,fQ2min);
  LOG("BodekYang", pDEBUG) << "Q2 at scaling var calculation = " << myQ2;

  double a  = TMath::Power( 2*kProtonMass*x, 2 ) / myQ2;
  double xw =  2*x*(myQ2+fB) / (myQ2*(1.+TMath::Sqrt(1+a)) +  2*fA*x);
  return xw;
}
//____________________________________________________________________________
void BYStrucFunc::KFactors(const Interaction * interaction,
	         double & kuv, double & kdv, double & kus, double & kds) const
{
// Overrides QPMDISStrucFuncBase::KFactors() to compute the BY K factors for
// u(valence), d(valence), u(sea), d(sea);

  double myQ2  = this->Q2(interaction);
  double GD  = 1. / TMath::Power(1.+myQ2/0.71, 2); // p elastic form factor
  double GD2 = TMath::Power(GD,2);

  kuv = (1.-GD2)*(myQ2+fCv2U)/(myQ2+fCv1U); // K - u(valence)
  kdv = (1.-GD2)*(myQ2+fCv2D)/(myQ2+fCv1D); // K - d(valence)
  kus = myQ2/(myQ2+fCsU);                   // K - u(sea)
  kds = myQ2/(myQ2+fCsD);                   // K - d(sea)
}