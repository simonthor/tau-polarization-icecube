# %% [markdown]
# ## Generate HepMC file
# I only store the tau lepton and the particles that creates it.

# %%
import csv
import pyhepmc
import pandas as pd
import numpy as np
from argparse import ArgumentParser


def csv2genevent(event: pd.DataFrame, nucleus_row=0, nutau_row=1, tau_row=2) -> pyhepmc.GenEvent:
    """This is a simplified version of the same function in plot_hepmc.
    In that file, all decay products are included in the decay chain (including the electron that is not in the csv file)
    Here, only the nucleus, neutrino and tau are included. This is because the tau will be decayed by Tauola."""
    evt = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
    
    # Extract particles        
    nucleus = pyhepmc.GenParticle(pyhepmc.FourVector(*event.iloc[nucleus_row, 3:6], event.iloc[nucleus_row, 2]), event.iloc[nucleus_row, 1], 3)
    nutau = pyhepmc.GenParticle(pyhepmc.FourVector(*event.iloc[nutau_row, 3:6], event.iloc[nutau_row, 2]), event.iloc[nutau_row, 1], 3)
    # Set status to 1, since we want Tauola to decay it
    tau = pyhepmc.GenParticle(pyhepmc.FourVector(*event.iloc[tau_row, 3:6], event.iloc[tau_row, 2]), event.iloc[tau_row, 1], 1)

    nucleus_out4m = nucleus.momentum + nutau.momentum - tau.momentum
    nucleus_out = pyhepmc.GenParticle(nucleus_out4m, 0, 1)
    # Create vertices
    interaction_vertex = pyhepmc.GenVertex()
    interaction_vertex.add_particle_in(nucleus)
    interaction_vertex.add_particle_in(nutau)
    interaction_vertex.add_particle_out(tau)
    interaction_vertex.add_particle_out(nucleus_out)
    evt.add_vertex(interaction_vertex)
    return evt


def main():
    # %% [markdown]
    # ## GENIE events
    # These are GENIE events that I simulated on my own

    # %%
    # Read energy from a config file
    argparser = ArgumentParser()
    argparser.add_argument("-i", type=str, required=True)  # Input file name
    argparser.add_argument("-od", type=str, required=True)  # Output file name for dat file
    argparser.add_argument("-oc", type=str, required=False, default=None)  # Output file name for csv file
    argparser.add_argument("--write-nevents", action="store_true", default=False)  # Write the number of events to the settings file (default: False)
    
    args = argparser.parse_args()

    input_filename = args.i
    output_dat = args.od
    output_csv = args.oc
    write_events = args.write_nevents

    # %%
    genie_events = pd.read_csv(input_filename)

    # %%
    # If the number of particles is not 3, then it means that the file contains decay products of the tau. 
    # We then only select the first 3 particles (incoming neutrino, nucleus, tau) in each event
    if not (genie_events.groupby("event_num").count() == 3).all().all():
        genie_events = genie_events.groupby("event_num").nth[:3]

    #%%
    if write_events:
        # Write the number of events to the yaml file at the end
        tau_n_events = genie_events["event_num"].nunique()
        
        with open("settings.yaml", "r") as f:
            file_contents = f.read()
        
        if "tau_n_events" in file_contents:
            old_tau_n_events = int(file_contents[file_contents.find("tau_n_events: ") + len("tau_n_events: "):])
            print(f"WARNING: The number of events is already written in the settings file: {old_tau_n_events}. Overwriting it with {tau_n_events}")
            file_contents_minus_row = file_contents[:file_contents.find("tau_n_events: ")]
            file_contents = file_contents_minus_row + f"tau_n_events: {tau_n_events}\n"
            with open("settings.yaml", "w") as f:
                f.write(file_contents)
        else:
            with open("settings.yaml", "a") as f:
                f.write(f"tau_n_events: {tau_n_events}\n")

    # %% [markdown]
    # ### Convert GENIE events to HepMC

    # %%
    with pyhepmc.open(output_dat, "w") as f:
        for event_num, event in genie_events.groupby("event_num"):
            evt = csv2genevent(event)
            evt.event_number = event_num
            f.write(evt)

    # %%
    # If no output csv file is given, we don't need to do anything more
    if output_csv is None:
        return
    
    # Open the file with decays and tau leptons, only select the lines with tau leptons and copy those to a new file
    with open(input_filename, "r") as f:
        with open(output_csv, "w") as f2:
            reader = csv.reader(f)
            writer = csv.writer(f2)
            writer.writerow(next(reader))
            for row in reader:
                if row[1] == "15":
                    writer.writerow(row)


if __name__ == "__main__":
    main()
