#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include "Pythia8/Pythia.h"

struct SimpleParticle {
    int event_num;
    int pdg;
    double E;
    double px;
    double py;
    double pz;
};

std::vector<SimpleParticle> readCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<SimpleParticle> particles;
    std::string line;

    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;

        SimpleParticle particle;
    
        std::getline(iss, token, ',');
        particle.event_num = std::stoi(token);

        std::getline(iss, token, ',');
        particle.pdg = std::stoi(token);

        std::getline(iss, token, ',');
        particle.E = std::stod(token);

        std::getline(iss, token, ',');
        particle.px = std::stod(token);

        std::getline(iss, token, ',');
        particle.py = std::stod(token);

        std::getline(iss, token);
        particle.pz = std::stod(token);

        particles.push_back(particle);
    }

    return particles;
}


int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " input_file output_file polarization" << std::endl;
        return 1;
    }

    const std::string input_filename = argv[1];
    const std::string output_filename = argv[2];

    const double polarization = std::stod(argv[3]);

    // Read particles from CSV
    std::vector<SimpleParticle> particles = readCSV(input_filename);

    // Initialize Pythia8
    Pythia8::Pythia pythia;

    // Set tau decays mode and polarization
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("ProcessLevel:resonanceDecays=on");
    pythia.readString("TauDecays:mode = 3");
    pythia.readString("TauDecays:tauPolarization = " + std::to_string(polarization));

    // Disable some decays
    pythia.readString("111:onMode = off");
    pythia.readString("310:onMode = off");
    pythia.readString("221:onMode = off");

    pythia.init();

    // Open output file in append mode
    std::ofstream outfile(output_filename);

    outfile << "event_num,pdg,E,px,py,pz\n";

    for (const auto& particle : particles) {
        // Write dummy nucleus and neutrino
        outfile << particle.event_num << "," << 1000080160 << "," 
                << 1 << "," << 0 << "," << 0 << "," << 0 << "\n";
        outfile << particle.event_num << "," << 16 << "," 
                << 1 << "," << 0 << "," << 0 << "," << 1 << "\n";
        // Write tau lepton
        outfile << particle.event_num << ","
                << particle.pdg << ","
                << std::setprecision(16) << particle.E << ","
                << std::setprecision(16) << particle.px << ","
                << std::setprecision(16) << particle.py << ","
                << std::setprecision(16) << particle.pz << "\n";

        // Reset Pythia8 event
        pythia.event.reset();

        // Calculate tau mass
        double tau_mass = std::sqrt(particle.E * particle.E - (particle.px * particle.px + particle.py * particle.py + particle.pz * particle.pz));

        // Add tau lepton to the event
        pythia.event.append(particle.pdg, 1, 0, 0, particle.px, particle.py, particle.pz, particle.E, tau_mass);

        // Generate the decay
        pythia.next();

        // Access the decay products
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Pythia8::Particle& decay_product = pythia.event[i];
            if (decay_product.isFinal()) {
                outfile << particle.event_num << ","
                        << decay_product.id() << ","
                        << std::setprecision(16) << decay_product.e() << ","
                        << std::setprecision(16) << decay_product.px() << ","
                        << std::setprecision(16) << decay_product.py() << ","
                        << std::setprecision(16) << decay_product.pz() << "\n";
            }
        }
    }

    pythia.stat();

    return 0;
}
