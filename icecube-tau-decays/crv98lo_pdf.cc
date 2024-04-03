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


int crv98lo_pdf() {
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
    double Q2 = 1.; // Q2 in GeV^2
    double M = 0.938; // Hit nucleon mass, in GeV
    double nuc_pdgc = 2112; // Hit nucleon PDG code. 2212 is a proton, 2112 is a neutron

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
    // will be 0 if < charm threshold
    double fuv_c = fPDFc -> UpValence();
    double fus_c = fPDFc -> UpSea();
    double fdv_c = fPDFc -> DownValence();
    double fds_c = fPDFc -> DownSea();    
    double fs_c  = fPDFc -> Strange();    
    double fc_c  = fPDFc -> Charm();

    // The above are the proton parton density function. Get the PDFs for the
    // hit nucleon (p or n) by swapping u<->d if necessary
    bool isP = pdg::IsProton  (nuc_pdgc);
    bool isN = pdg::IsNeutron (nuc_pdgc);
    assert(isP || isN);

    double tmp = 0;
    if (isN) {  // swap u <-> d, c <-> s
        tmp = fuv;   fuv   = fdv;   fdv   = tmp;
        tmp = fus;   fus   = fds;   fds   = tmp;
        tmp = fuv_c; fuv_c = fdv_c; fdv_c = tmp;
        tmp = fus_c; fus_c = fds_c; fds_c = tmp;
    }

    std::cout << "fuv: " << fuv << " fus: " << fus << " fdv: " << fdv << " fds: " << fds << " fs: " << fs << " fc: " << fc << std::endl;
    std::cout << "u = fuv + fus = " << fuv + fus << " d = fdv + fds = " << fdv + fds << " s = " << fs << " c = " << fc << std::endl;
    
    std::cout << "fuv_c: " << fuv_c << " fus_c: " << fus_c << " fdv_c: " << fdv_c << " fds_c: " << fds_c << " fs_c: " << fs_c << " fc_c: " << fc_c << std::endl;
    std::cout << "u_c = fuv_c + fus_c = " << fuv_c + fus_c << " d = fdv_c + fds_c = " << fdv_c + fds_c << " s_c = " << fs_c << " c_c = " << fc_c << std::endl;
    
    return 0;
}
