# %% [markdown]
# ## Generate HepMC file
# I only store the tau lepton and the particles that creates it.

# %%
import csv
import pyhepmc
import pandas as pd
import numpy as np
import sys
sys.path.append("../")
from utils.hepmc_io import csv2genevent
from argparse import ArgumentParser


def main():
    # %% [markdown]
    # ## GENIE events
    # These are GENIE events that I simulated on my own

    # %%
    # Read energy from a config file
    argparser = ArgumentParser()
    argparser.add_argument("energy", type=str)
    args = argparser.parse_args()

    genie_nutau_energy = args.energy

    # %%
    genie_events = pd.read_csv(f"../data/genie_tau_pol_data_e{genie_nutau_energy}.csv")

    # %%
    assert (genie_events.groupby("event_num").count() == 3).all().all(), "ERROR: not all events have 3 particles. Quitting"

    #%%
    # Write the number of events to the yaml file at the end
    tau_n_events = genie_events["event_num"].unique().size
    
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
    with pyhepmc.open(f"../data/tauola_input_genie_e{genie_nutau_energy}.dat", "w") as f:
        for event_num, event in genie_events.groupby("event_num"):
            evt = csv2genevent(event)
            evt.event_number = event_num
            f.write(evt)

    # %% [markdown]
    # Save polarization in a separate file

    # %%
    # Open the file with decays and tau leptons, only select the lines with tau leptons and copy those to a new file
    with open(f"../data/genie_tau_pol_data_e{genie_nutau_energy}.csv", "r") as f:
        with open(f"../data/genie_pol_e{genie_nutau_energy}.csv", "w") as f2:
            reader = csv.reader(f)
            writer = csv.writer(f2)
            writer.writerow(next(reader))
            for row in reader:
                if row[1] == "15":
                    writer.writerow(row)


if __name__ == "__main__":
    main()
