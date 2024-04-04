#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>

#include "Physics/PartonDistributions/PDFModelI.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/PartonDistributions/GRV98LO.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
// #include "Physics/NuclearState/NuclearUtils.h"
// #include "Framework/Utils/PhysUtils.h"

using namespace genie;


int crv98lo_pdf_eval(std::string output_file, double Q2) {
    const PDFModelI * pdf_model =
         dynamic_cast<const PDFModelI *> (new GRV98LO());
    PDF* fPDF =  new PDF();
    
    fPDF->SetModel(pdf_model);

    // Source: https://github.com/GENIE-MC/Generator/blob/1817b89cf815c0694c187f874288f1c1be1e712a/config/GRV98LO.xml#L11
    const double fQ2min = 0.800;
    // Charm mass, in GeV, according to Wikipedia
    const double fMc = 1.27;
    
    // Open output file and write header to it
    std::ofstream ofile(output_file);
    ofile << "x,Q2,fuv,fus,fdv,fds,fs,fc" << std::endl;

    // Iterate over an x from 0 to 1 with steps of 0.01
    for (double x = 0; x <= 1; x += 0.01) {
        // Reset PDF
        fPDF->Reset();
        // Calculate PDF at (x, Q2)
        double Q2pdf =  TMath::Max(Q2, fQ2min);
        fPDF->Calculate(x, Q2pdf);
        
        double fuv =  fPDF  -> UpValence();
        double fus   = fPDF  -> UpSea();
        double fdv   = fPDF  -> DownValence();
        double fds   = fPDF  -> DownSea();
        double fs    = fPDF  -> Strange();
        double fc    = fPDF -> Charm();
        // Write to output file
        ofile << x << "," << Q2 << "," << fuv << "," << fus << "," << fdv << "," << fds << "," << fs << "," << fc << std::endl;
    }

    ofile.close();
    delete fPDF;
    delete pdf_model;
    
    return 0;
}
