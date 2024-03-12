from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


def compare_histos(
        nutau, nutau_nopol, nutau_g4, bins, 
        labels=("polarized", "unpolarized (Tauola)", "unpolarized (IceCube)"), 
        density=None, ax=None, **kwargs):
    """Plot 3 histograms on the same axis."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4), layout="constrained")
    else:
        fig = ax.get_figure()

    for energies, particle_type in zip((nutau, nutau_nopol, nutau_g4), labels):
        ax.hist(energies, bins=bins, label=particle_type, density=density, histtype="step", lw=2)
        # TODO add error bars

    ax.set(**kwargs)
    ax.grid(True, alpha=0.5)
    ax.legend(fontsize="large")

    return fig, ax


def filter_events(decay_products: pd.DataFrame, col: str, filter_func: callable, **kwargs) -> pd.DataFrame:
    """Filter events based on the function filter_func.
    The function must take values, index as the two arguments and return a boolean value.
    The function will be applied on an array of the values in the DataFrame column col and once per group, 
    where the grouping is performed over event_num."""
    mask = decay_products.groupby("event_num")[col].agg(filter_func, **kwargs)
    allowed_events = mask[mask > 0].index
    return decay_products[decay_products["event_num"].isin(allowed_events)]


def plot_histograms(
    tauola: dict[int, pd.DataFrame], tauola_nopol: dict[int, pd.DataFrame], icecube: dict[int, pd.DataFrame], /, *,
    bins: dict[int, np.ndarray], filter_func: callable, plot_func: callable, **kwargs):
    """Plot several subplots, each one containing three histograms. 
    Different subplots correspond to different keys in the dicts (typically incoming neutrino energy)."""
    
    fig, axs = plt.subplots(ncols=len(tauola), figsize=(4*len(tauola), 4), layout="constrained")

    for e, ax in zip(tauola, axs):
        b = bins[e]
        decay_products_e = tauola[e]
        decay_products_nopol_e = tauola_nopol[e]
        decay_products_ic_e = icecube[e]

        # Only select a certain decay mode
        selected_events = filter_events(decay_products_e, "pdg", filter_func, engine="numba")
        selected_events_nopol = filter_events(decay_products_nopol_e, "pdg", filter_func, engine="numba")
        selected_events_ic = filter_events(decay_products_ic_e, "pdg", filter_func, engine="numba")

        # Plot the momentum fraction as a histogram from 0 to 1
        compare_histos(
            plot_func(selected_events), plot_func(selected_events_nopol), plot_func(selected_events_ic),
            ax=ax, bins=b, **kwargs,
        )
    return fig, axs