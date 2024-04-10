# %% [markdown]
# ## Generate HepMC file
# I only store the tau lepton and the particles that creates it.

# %%
import csv
import pyhepmc
import pandas as pd
import numpy as np
from argparse import ArgumentParser


def main():
    # %% [markdown]
    # ## GENIE events
    # These are GENIE events that I simulated on my own

    # %%
    # Read energy from a config file
    argparser = ArgumentParser()
    argparser.add_argument("-i", type=str, required=True)  # Input file name
    argparser.add_argument("-od", type=str, required=True)  # Output file name for dat file
    
    args = argparser.parse_args()

    input_filename = args.i
    output_dat = args.od
    
    # %%
    events = pd.read_csv(input_filename).query("pdg == 15")

    # %% [markdown]
    # ### Convert events to HepMC

    # %%
    with pyhepmc.open(output_dat, "w") as f:
        for event_num, event in events.groupby("event_num"):
            evt = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
            
            # Set status to 1, since we want Tauola to decay it
            tau = pyhepmc.GenParticle(pyhepmc.FourVector(*event.iloc[0, 3:6], event.iloc[0, 2]), event.iloc[0, 1], 1)

            # Create vertices
            evt.add_particle(tau)

            evt.event_number = event_num

            f.write(evt)


if __name__ == "__main__":
    main()
