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


int grv98lo_pdf(std::string input_file, std::string output_file, bool charm_correction = false) {
    BYPDF * bypdf = new BYPDF();
    bypdf->Configure("Default");

    const PDFModelI * pdf_model =
         dynamic_cast<const PDFModelI *> (bypdf);
    
    PDF* fPDF =  new PDF();
    fPDF->SetModel(pdf_model);
    
    // Open input file, and read the first line. 
    std::ifstream file(input_file);
    std::string line, word;
    std::getline(file, line);
    std::stringstream s(line);
    int index = 0;
    
    int xcol, Q2col, nuc_pdgcol, discol;
    std::string xcol_name = charm_correction ? "xs" : "xi";

    // Identify the index where the value is "x" and "Q2"
    while (std::getline(s, word, ',')) {
        // std::cout << word << std::endl;
        if (word == xcol_name) {
            xcol = index;
        } else if (word == "Q2s") {
            Q2col = index;
        } else if (word == "hitnuc") {
            nuc_pdgcol = index;
        } else if (word == "dis") {
            discol = index;
        }
        index++;
    }
    std::cout << "xcol: " << xcol << ", Q2col: " << Q2col << ", nuc_pdgcol: " << nuc_pdgcol << ", discol: " << discol << std::endl;

    // Open output file and write header to it
    std::ofstream ofile(output_file);
    ofile << line << ",fuv,fus,fdv,fds,fs,fc" << std::endl;

    double x, Q2;
    int nuc_pdgc;
    bool dis;
    // Read the rest of the lines and compute the PDF at each (x, Q2)
    while (std::getline(file, line)) {
        std::stringstream s(line);
        index = 0;
        while (std::getline(s, word, ',')) {
            if (index == xcol) {
                x = std::stod(word);
            } else if (index == Q2col) {
                Q2 = std::stod(word);
            } else if (index == nuc_pdgcol) {
                nuc_pdgc = std::stoi(word);
            } else if (index == discol) {
                // Set dis to True or False
                dis = word == "True";
            }
            index++;
        }

        if (!dis) {
            ofile << line << ",,,,,,\n";
            continue;
        }

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
        // If collision is with neutron, swap u<->d
        bool isP = pdg::IsProton  (nuc_pdgc);
        bool isN = pdg::IsNeutron (nuc_pdgc);
        assert(isP || isN);
        double tmp = 0;
        if (isN) {
            tmp = fuv;   fuv   = fdv;   fdv   = tmp;
            tmp = fus;   fus   = fds;   fds   = tmp;
        }
        // Write to output file
        ofile << line << "," << fuv << "," << fus << "," << fdv << "," << fds << "," << fs << "," << fc << std::endl;
    }

    ofile.close();
    file.close();
    delete fPDF;
    delete pdf_model;
    
    return 0;
}
