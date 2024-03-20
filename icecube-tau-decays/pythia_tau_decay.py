import pythia8
import numpy as np
import pandas as pd
import argparse


def main():
    # Get file name to read from and to write to from command line
    parser = argparse.ArgumentParser(description="Simulate tau decays using Pythia8")
    parser.add_argument("input_file", help="CSV file containing particle information")
    parser.add_argument("output_file", help="CSV file to write decay information to")
    parser.add_argument("-p", "--polarization",
                        dest="p",
                        type=float,
                        help="Polarization in z direction. It might be that -1, 0, or 1 are the only possible values, or maybe floats are possible", 
                        default=0)
    
    args = parser.parse_args()

    input_filename = args.input_file
    output_filename = args.output_file
    polarization = args.p

    tau_pdg_id = 15  # PDG ID for tau lepton

    # Load particles simulated by IceCube
    particle_info = pd.read_csv(input_filename)
    
    # Initialize Pythia
    pythia = pythia8.Pythia()
    # pythia.readString("ProcessLevel:resonanceDecays=on")
    # this is the trick to make Pythia8 work as "decayer"
    pythia.readString("ProcessLevel:all = off")
    pythia.readString("ProcessLevel:resonanceDecays=on")

    # shut off Pythia8 (default) verbosity
    # fDecayer->readString("Init:showAllSettings=false");
    # fDecayer->readString("Init:showChangedSettings=false");
    # fDecayer->readString("Init:showAllParticleData=false");
    # fDecayer->readString("Init:showChangedParticleData=false");
    
    # specify how many Py8 events to print out, at either level
    # in this particular case print out a maximum of 10 events
    # fDecayer->readString("Next:numberShowProcess = 0" );
    # fDecayer->readString("Next:numberShowEvent = 10");
    
    # Disable pion decay
    pythia.readString("111:onMode = off")
    pythia.readString("310:onMode = off")

    pythia.init()
    with open(output_filename, "w") as f:
        f.write("event_num,pdg,E,px,py,pz\n")

    for event_num, particles in particle_info.groupby("event_num"):
        tau = particles[particles["pdg"] == tau_pdg_id].iloc[0, :]

        for particle in particles.iloc[:3, :].itertuples():
            with open(output_filename, "a") as f:
                f.write(f"{particle.event_num:.0f},{particle.pdg:.0f},{particle.E},{particle.px},{particle.py},{particle.pz}\n")

        # Initialize the event
        pythia.event.reset()

        # Calculate tau mass from 4-momentum
        tau_mass = np.sqrt(tau.E**2 - (tau.px**2 + tau.py**2 + tau.pz**2))
        # Add tau lepton to the event
        pythia.event.append(
            tau_pdg_id, 1, 0, 0, 
            #  px      py      pz      E       m
            tau.px, tau.py, tau.pz, tau.E, tau_mass
        )
        
        # Set polarization
        pythia.event.back().pol(polarization)

        # Generate the decay
        pythia.next()

        # Access the decay products
        event_particles = pythia.event
        
        # print(f"Event number: {tau.event_num}")
        # Print information about the decay products
        for i, particle in enumerate(event_particles):
            if not particle.isFinal():
                continue

            with open(output_filename, "a") as f:
                f.write(f"{tau.event_num:.0f},{particle.id():.0f},{particle.e()},{particle.px()},{particle.py()},{particle.pz()}\n")
        
    # End Pythia
    pythia.stat()

    # Convert the decay information to a DataFrame
    # pythia_decays_df = pd.DataFrame(pythia_decays, columns=["event_num", "pdg", "E", "px", "py", "pz"])
    # pythia_decays_df.to_csv(output_filename, index=False)


if __name__ == "__main__":
    main()



