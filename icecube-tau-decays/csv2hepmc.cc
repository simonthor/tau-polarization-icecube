// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @example HepMC3_fileIO_example.cc
 *  @brief Test of file I/O
 *
 *  Parses HepMC3 file and saves it as a new HepMC3 file.
 *  The resulting file should be an exact copy of the input file
 *
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include <iostream>
using namespace HepMC3;

/** Main program */
int main(int argc, char **argv) {

    if( argc<3 ) {
        std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> <output_file>" << std::endl;
        exit(-1);
    }

    WriterAscii output_file(argv[2]);
    // Open a csv file for reading and iterate over the rows until a blank line is encountered
    ifstream input_file(argv[1]);
    vector<GenParticlePtr> particles;
    vector<GenVertexPtr> vertices;
    GenEvent evt(Units::GEV,Units::MM);
    while (input_file.good()) {
        string line;
        getline(input_file, line);

        // If the line is not empty, parse it
        // Create a vector with the values of the line
        vector<string> values;
        stringstream ss(line);
        string value;
        while (getline(ss, value, ',')) {
            values.push_back(value);
        }
        GenParticlePtr p1 = std::make_shared<GenParticle>( FourVector( 0.0,    0.0,   7000.0,  7000.0  ),2212,  3 );
        
        // The values are energy, px, py, pz, pdgid, status

        if (line == "") {
            particles.clear();
            vertices.clear();
            GenEvent evt(Units::GEV,Units::MM);
        }
    }

    int events_parsed = 0;

    // while(!input_file.failed()) {

        //                                                               px      py        pz       e     pdgid status
        GenParticlePtr p2 = std::make_shared<GenParticle>( FourVector( 0.750, -1.569,   32.191,  32.238),   1,  3 );
        GenParticlePtr p3 = std::make_shared<GenParticle>( FourVector( 0.0,    0.0,  -7000.0,  7000.0  ),2212,  3 );
        GenParticlePtr p4 = std::make_shared<GenParticle>( FourVector(-3.047,-19.0,    -54.629,  57.920),  -2,  3 );

        GenVertexPtr v1 = std::make_shared<GenVertex>();
        v1->add_particle_in (p1);
        v1->add_particle_out(p2);
        evt.add_vertex(v1);

        // Set vertex status if needed
        v1->set_status(4);

        GenVertexPtr v2 = std::make_shared<GenVertex>();
        v2->add_particle_in (p3);
        v2->add_particle_out(p4);
        evt.add_vertex(v2);

        GenVertexPtr v3 = std::make_shared<GenVertex>();
        v3->add_particle_in(p2);
        v3->add_particle_in(p4);
        evt.add_vertex(v3);

        GenParticlePtr p5 = std::make_shared<GenParticle>( FourVector(-3.813,  0.113, -1.833, 4.233),  22, 1 );
        GenParticlePtr p6 = std::make_shared<GenParticle>( FourVector( 1.517,-20.68, -20.605,85.925), -24, 3 );

        v3->add_particle_out(p5);
        v3->add_particle_out(p6);

        GenVertexPtr v4 =std:: make_shared<GenVertex>();
        v4->add_particle_in (p6);
        evt.add_vertex(v4);

        GenParticlePtr p7 = std::make_shared<GenParticle>( FourVector(-2.445, 28.816,  6.082,29.552),  1, 1 );
        GenParticlePtr p8 = std::make_shared<GenParticle>( FourVector( 3.962,-49.498,-26.687,56.373), -2, 1 );

        v4->add_particle_out(p7);
        v4->add_particle_out(p8);

        // If reading failed - exit loop
        // if( input_file.failed() ) break;

        // Save event to output file
        output_file.write_event(evt);

        ++events_parsed;
        if( events_parsed%100 == 0 ) std::cout<<"Events parsed: "<<events_parsed<<std::endl;
    // }

    // input_file.close();
    output_file.close();

    return 0;
}
