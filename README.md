# Master thesis - tau polarization effects on tau neutrino detection
This repository contains the code for my master thesis on investigating the effects on tau polarization on IceCube measurements of tau neutrinos.

## Repository structure
### recreate-tau-appearance
This folder contains the code to recreate the [IceCube tau neutrino appearance analysis](https://arxiv.org/abs/1901.05366).

The .csv files are from the publicly released data and the [hyperplanes_example.ipynb](recreate-tau-appearance/hyperplanes_example.ipynb) notebook was included as well in the public data release. This notebook shows how to use the hyperplanes to calculate the tau neutrino appearance.

The [recreate.ipynb](recreate-tau-appearance/recreate.ipynb) notebook contains the code to recreate several plots from the analysis, mainly the tau neutrino normalization $\chi^2$ plot. The notebook is an extension of the [hyperplanes_example.ipynb](recreate-tau-appearance/hyperplanes_example.ipynb) notebook, with additional code added. There should however not be any difference regarding the detector efficiencies and application of the nuisance parameters.

### learn-tauola
Code for learning, installing, and testing Tauola.

### learn-geant4
Code for (re)learning and testing Geant4.
