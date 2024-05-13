from collections import Counter
import pandas as pd


def branching_ratios(decay_products: pd.DataFrame) -> dict[tuple[int], float]:
    """Compute branching ratios for the decay products of the tau lepton.
    NOTE: the input dataframe should only contain decay products, no taus, incoming neutrinos etc."""
    n_taus = decay_products.query("pdg == 16").shape[0]
    n_anti_taus = decay_products.query("pdg == -16").shape[0]
    c = Counter(tuple(sorted(a.tolist())) for i, a in decay_products.groupby("event_num")["pdg"])
    
    br = {}
    
    for pdgs, n in c.items():
        if 16 in pdgs:
            br[pdgs] = n / n_taus
        elif -16 in pdgs:
            br[pdgs] = n / n_anti_taus
        else:
            raise ValueError("No tau in the event")
    
    # Sort based on the branching ratio. Highest first
    br = dict(sorted(br.items(), key=lambda item: item[1], reverse=True))
    return br