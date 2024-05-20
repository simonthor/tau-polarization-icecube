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
decay_flags=""
while IFS= read -r line; do
    echo $line
    if [[ $line == "energy:" ]]; then
        inside_energy_list=1
    elif [[ $inside_energy_list -eq 1 ]]; then
        if [[ $line == "  - "* ]]; then
            # Only select the part of the line after the "  - " and add it to the list
            energy_list+=(${line:4})
            echo "Added to energy list"
        else
            inside_energy_list=0
            echo "end energy list"
        fi
    elif [[ $line == start_step* ]]; then
        start_step=$(echo $line | cut -d' ' -f2)
        echo "Set start step"
    elif [[ $line == "rad: off" ]]; then
        decay_flags="-r $decay_flags"
    elif [[ $line == "boost: on" ]]; then
        decay_flags="-b $decay_flags"
    fi
done < tauola_settings.yaml

# Remove the last character from decay_flags, if it is not empty
if [[ ${#decay_flags} -gt 0 ]]; then
    decay_flags=${decay_flags::-1}
fi

echo "Energy list: ${energy_list[@]}"
echo "inside_energy_list: $inside_energy_list"
echo "Start step: $start_step"
echo "Flags passed to decay_taus.cc: $decay_flags"

# Make file end the same as decay_flags but remove all spaces
file_end=$(echo $decay_flags | tr -d ' ')

# iterate over all energies in energy_list
for energy in "${energy_list[@]}"; do
    # Define input csv file name
    input_csv_file=../data/test_genie_NuTau_$energy.0_GeV_particles.csv
    input_csv_event_file=../data/test_genie_NuTau_$energy.0_GeV_event_info.csv
    preprocessed_csv_event_file=../data/test_genie_NuTau_${energy}.0_GeV_event_info_dis.csv
    taupol_input_file=../data/test_genie_NuTau_${energy}.0_GeV_event_info_sig.csv
    # input_csv_file=../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_particles_extended.csv
    # Define output csv file name
    output_csv_file=../data/NuTau_$energy.0_GeV_tau.csv
    # output_csv_file=../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_tau.csv
    # Define Tauola input dat file name
    input_dat_file=../data/NuTau_$energy.0_GeV_tauola_input.dat
    # input_dat_file=../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_tauola_input.dat
    # Define Tauola output dat file name
    output_dat_file=../data/NuTau_$energy.0_GeV_tauola_output$file_end.dat
    # output_dat_file=../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_tauola_output$file_end.dat
    # Define Tauola output dat file name without polarization
    output_dat_file_nopol=../data/NuTau_$energy.0_GeV_tauola_output_nopol$file_end.dat
    # output_dat_file_nopol=../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_tauola_output_nopol$file_end.dat
    # Define Tauola output dat file name with full polarization
    output_dat_file_lpol=../data/NuTau_$energy.0_GeV_tauola_output_lpol$file_end.dat
    # output_dat_file_lpol=../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_tauola_output_lpol$file_end.dat
    # output_dat_file_rpol=../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_tauola_output_rpol$file_end.dat

    if [ $start_step -lt 2 ]; then
        echo "Converting GENIE csv file to dat file..."
        # Run python script to convert the file into a dat file, and generate a csv file with only the tau leptons 
        python convert_genie.py -i $input_csv_file -od $input_dat_file
    fi
    continue
    
    if [ $start_step -lt 3 ]; then
        echo "Preparing the file for calculating the polarization..."
        python preprocess_event_info.py -ip $input_csv_file -ie $input_csv_event_file -o $preprocessed_csv_event_file
        # Calculate the GRV98LO PDF values for each DIS event from x and Q^2
        #  This adds 6 new columns: fuv, fus, fdv, fds, fs, fc
        #  Each column corresponds to the PDF value for the event for the valence up quark, sea up quark, valence down quark, sea down quark, strange quark, and charm quark respectively.
        #  These are used for the polarization calculations
        #  When the last input argument is true, the charm-corrected x value is used, as described in Hagiwara et al. (2003). This is the default.
        genie -l -b -q "grv98lo_pdf.C(\"${preprocessed_csv_event_file}\",\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_pdf.csv\",true)"
        # Calculate the Sigma values of for RES events using the Berger-Sehgal model
        #  This adds two new columns: sigmm, sigpp
        #  The longitudinal polarization is then given by (sigmm - sigpp) / (sigmm + sigpp)
        #  The goal is to also add a way of calculating the transverse polarization, but this seems a bit more complicated
        #  as it is not relevant for the GENIE cross section calculations, and therefore not implemented in the library.
        genie -l -b -q "analytic_tau_pol_res_int.C(\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_pdf.csv\",\"${taupol_input_file}\")"
    fi

    if [ $start_step -lt 4 ]; then
        echo "Calculating tau polarization for all events..."
        python tau_polarization.py -ip $input_csv_file -ie $taupol_input_file -o $output_csv_file
    fi

    if [ $start_step -lt 5 ]; then
        echo "Running Tauola tau decay simulation with realistic polarization..."
        # Run the Tauola tau decay simulation, with polarization. The columns are currently set to 1, 2, 3, as given as output from tau_polarization.py
        ./decay.o $input_dat_file $output_dat_file 1 2 3 $output_csv_file $decay_flags &> ../logfiles/icecube_tauola_run_e$energy$file_end.log
    fi

    if [ $start_step -lt 6 ]; then
        echo "Running Tauola tau decay simulation without polarization..."
        # Run the Tauola tau decay simulation, without polarization
        ./decay.o $input_dat_file $output_dat_file_nopol 0 0 0 $decay_flags &> ../logfiles/icecube_tauola_run_e${energy}_nopol.log
    fi

    if [ $start_step -lt 7 ]; then
        echo "Running Tauola tau decay simulation with fully left-handed polarization..."
        # Run the Tauola tau decay simulation, without polarization
        ./decay.o $input_dat_file $output_dat_file_lpol 0 0 -1 $decay_flags &> ../logfiles/icecube_tauola_run_e${energy}_lpol.log
    fi
    
    # if [ $start_step -lt 6 ]; then
    #     echo "Running Tauola tau decay simulation with fully right-handed polarization..."
    #     # Run the Tauola tau decay simulation, without polarization
    #     ./decay.o $input_dat_file $output_dat_file_rpol 0 0 1 $decay_flags &> ../logfiles/icecube_tauola_run_e${energy}_rpol.log
    # fi

    # if [ $start_step -lt 6 ]; then
    #     echo "Running Pythia tau decay simulation without polarization..."
    #     # Run the Tauola tau decay simulation, without polarization
    #     ./pythia_decay.o $output_csv_file ../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_pythia_output_nopol.csv 0 &> ../logfiles/pythia_tau_decays_e${energy}_nopol.log
    # fi

    # if [ $start_step -lt 7 ]; then
    #     echo "Running Pythia tau decay simulation with left-handed polarization..."
    #     # Run the Tauola tau decay simulation, without polarization
    #     ./pythia_decay.o $output_csv_file ../data/test_bare_lepton_toy_Tau_000005_${energy}.0_GeV_pythia_output_lpol.csv -1 &> ../logfiles/pythia_tau_decays_e${energy}_lpol.log
    # fi

done
