import pyhepmc
import numpy as np
import pandas as pd


def load_hepmc(filename: str) -> pd.DataFrame:
    particles = []
    with pyhepmc.open(filename, "r") as f:
        # Iterate over all events
        for i, evt in enumerate(f):
            incoming_particles = []
            outgoing_particles = []
            # Find the tau and identify its daughter tau neutrino
            for vertex in evt.vertices:
                event_number = evt.event_number
                # Check if there is a tau neutrino as incoming particle. Then store all incoming particles and the outgoing tau
                if any(np.abs(p.pid) == 16 for p in vertex.particles_in):
                    for mother in sorted(vertex.particles_in, key=lambda x: x.pid, reverse=True):
                        outgoing_particles.append([event_number, mother.pid, mother.momentum.e, mother.momentum.px, mother.momentum.py, mother.momentum.pz])
                    for dauther in vertex.particles_out:
                        if np.abs(dauther.pid) == 15:
                            outgoing_particles.append([event_number, dauther.pid, dauther.momentum.e, dauther.momentum.px, dauther.momentum.py, dauther.momentum.pz])
                    continue

                # Find the final state decay particles
                for daughter in vertex.particles_out:
                    # If it is not a final state particle, skip it
                    if daughter.status != 1:
                        continue

                    outgoing_particles.append(
                        [
                            event_number,
                            daughter.pid,
                            daughter.momentum.e,
                            daughter.momentum.px,
                            daughter.momentum.py,
                            daughter.momentum.pz
                        ]
                    )
            particles.extend(incoming_particles)
            particles.extend(outgoing_particles)

            if i % 10_000 == 0:
                print(i)

    return pd.DataFrame(particles, columns=['event_num', 'pdg', 'E', 'px', 'py', 'pz'])


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