#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>

#include "Physics/PartonDistributions/PDFModelI.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/DeepInelastic/XSection/BYPDF.h"
#include "Physics/PartonDistributions/GRV98LO.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
// #include "Physics/NuclearState/NuclearUtils.h"
// #include "Framework/Utils/PhysUtils.h"

using namespace genie;


int bypdf(std::string output_file, double Q2) {
    BYPDF * bypdf = new BYPDF();
    bypdf->Configure("Default");

    const PDFModelI * pdf_model =
         dynamic_cast<const PDFModelI *> (bypdf);
    // std::cout << pdf_model->GetConfig() << "\n";
    PDF* fPDF =  new PDF();
    // std::cout << "Created PDF\n";

    fPDF->SetModel(pdf_model);

    // Charm mass, in GeV, according to Wikipedia
    const double fMc = 1.27;
    
    // Open output file and write header to it
    std::ofstream ofile(output_file);
    ofile << "x,Q2,fuv,fus,fdv,fds,fs,fc" << std::endl;
    
    // std::cout << "Opened file\n";

    // Iterate over an x from 0 to 1 with steps of 0.01
    for (double x = 0; x <= 1; x += 0.01) {
        // Reset PDF
        fPDF->Reset();
        // Calculate PDF at (x, Q2)
        fPDF->Calculate(x, Q2);
            
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
