#!/bin/bash
# Bash strict mode:
set -euo pipefail
IFS=$'\n\t'
# End of strict mode

# Read settings.yaml to get all the global parameters
# Open the file and iterate over the lines
# Create a list with all energies stored in it. 
# The energies are listed one after each other in the settings file, with the first line being "energy:"
# and the subsequent lines starting with several spaces, followed by "- " and then the actual energy value
echo "Reading tauola_settings.yaml..."
energy_list=()
inside_energy_list=0
start_step=0
while IFS= read -r line; do
    if [[ $line == "energy:" ]]; then
        inside_energy_list=1
    elif [[ $inside_energy_list -eq 1 ]]; then
        if [[ $line == "  - "* ]]; then
            # Only select the part of the line after the "  - " and add it to the list
            energy_list+=(${line:4})
        else
            inside_energy_list=0
        fi
    elif [[ $line == start_step* ]]; then
        start_step=$(echo $line | cut -d' ' -f2)
    fi
done < tauola_settings.yaml

echo "Energy list: ${energy_list[@]}"
echo "inside_energy_list: $inside_energy_list"
echo "Start step: $start_step"

# iterate over all energies in energy_list
for energy in "${energy_list[@]}"; do
    # Define input csv file name
    input_csv_file=../data/test_genie_NuTau_$energy.0_GeV_particles.csv
    # Define output csv file name
    output_csv_file=../data/NuTau_$energy.0_GeV_tau.csv
    # Define Tauola input dat file name
    input_dat_file=../data/NuTau_$energy.0_GeV_tauola_input.dat
    # Define Tauola output dat file name
    output_dat_file=../data/NuTau_$energy.0_GeV_tauola_output.dat
    # Define Tauola output dat file name without polarization
    output_dat_file_nopol=../data/NuTau_$energy.0_GeV_tauola_output_nopol.dat

    if [ $start_step -lt 2 ]; then
        echo "Converting GENIE csv file to dat file..."
        # Run python script to convert the file into a dat file, and generate a csv file with only the tau leptons 
        python convert_genie.py -i $input_csv_file -oc $output_csv_file -od $input_dat_file
    fi

    echo "Running Tauola tau decay simulation with polarization..."
    # Run the Tauola tau decay simulation, with polarization
    ./decay.o $input_dat_file $output_dat_file 6 7 8 $output_csv_file &> icecube_tauola_run_e$energy.log

    echo "Running Tauola tau decay simulation without polarization..."
    # Run the Tauola tau decay simulation, without polarization
    ./decay.o $input_dat_file $output_dat_file_nopol 0 0 0 &> icecube_tauola_run_e${energy}_nopol.log
done
