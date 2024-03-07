# Compare tau decays between IceCube and Tauola
This folder contains most of the code necessary for generating and analyzing events for tau decays from IceCube and Tauola, where the tau lepton has been created from a tau neutrino interacting with ice. The code is written in C++ and Python.

## Workflow
The programs used in this directory should be executed in the following order:

![workflow](workflow.png)

The arrows show the data flow and the file names of the input/output files.

All of these steps do not have to be executed manually. Instead, one can run `generate_events.sh`, which runs all of these programs sequentially. The only thing that needs to be changed is the settings.yaml, which contains all settings for the programs. The parameters are described below:
<!-- Create a table with 3 columns, one with the parameter name, one with the type, and one for the description -->
| Parameter | Type | Description |
| --- | --- | --- |
| energy | int/float | Energy of the neutrino in GeV |
| genie_n_events | int | Number of events to generate with GENIE |
| run | int | Run number of the GENIE events. Included in the GENIE output file name. |
| pdg | int | PDG code of the neutrino. Either 16 or -16 |
| start_step | int (1-5) | Step to start from. All steps before this will be skipped. |
| tau_n_events | int | Number of events to generate with Tauola. This will be calculated automatically by a program that counts the number of GENIE events that produced a tau lepton, and should therefore usually not be changed by the user. |

The start step values correspond to these parts of the workflow:
Indicates which process to start from. 
1. from beginning.
2. skip GENIE.
3. skip conversion from root to csv.
4. skip conversion from csv to HepMC file and tau csv file generation.
5. skip Tauola simulations.