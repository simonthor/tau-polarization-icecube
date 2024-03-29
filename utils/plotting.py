from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


def compare_histos(
        datasets: dict[str, np.ndarray], bins: dict[str, np.ndarray],  
        density=None, ax=None, errorbar=False, **kwargs):
    """Plot 3 histograms on the same axis."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4), layout="constrained")
    else:
        fig = ax.get_figure()

    for label, v in datasets.items():
        values, _, polygons = ax.hist(v, bins=bins, label=label, density=density, histtype="step", lw=2)
        # The error bars should be the sqrt of the number of events in each bin. 
        # If density is used, the error bars should be 
        # sqrt(counts) / bin_width / sum(counts) = sqrt(values) / sqrt(sum(counts) * bin_width)
        if errorbar:
            color = polygons[0].get_edgecolor()
            ax.errorbar(
                (bins[1:] + bins[:-1]) / 2, values, 
                yerr=np.sqrt(values) / np.sqrt(np.sum(v.size) * np.diff(bins)) if density else np.sqrt(values),
                fmt="none", capsize=2, capthick=2, elinewidth=2, ecolor=color,
            )

    ax.set(**kwargs)
    ax.grid(True, alpha=0.5)
    # ax.legend()

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
    datasets: dict[str, dict[int, pd.DataFrame]], /, *,
    bins: dict[int, np.ndarray], filter_func: callable, plot_func: callable, **kwargs):
    """Plot several subplots, each one containing three histograms. 
    Different subplots correspond to different keys in the dicts (typically incoming neutrino energy)."""
    
    d1 = list(datasets.values())[0]
    fig, axs = plt.subplots(ncols=len(d1), figsize=(4*len(d1), 4), layout="constrained")

    for i, (e, ax) in enumerate(zip(d1, axs)):
        b = bins[e]
        events_to_plot = {}
        for label, df in datasets.items():
            # Only select a certain decay mode
            selected_events = filter_events(df[e], "pdg", filter_func, engine="numba")
            events_to_plot[label] = plot_func(selected_events)

        # Plot the momentum fraction as a histogram from 0 to 1
        compare_histos(events_to_plot, ax=ax, bins=b, **kwargs)
        ax.set_title(f"$E_\\nu = {e}$ GeV")
        if i==0:
            ax.legend()
        
    return fig, axs