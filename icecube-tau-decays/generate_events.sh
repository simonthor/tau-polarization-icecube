#!/bin/bash
# Bash strict mode:
set -euo pipefail
IFS=$'\n\t'
# End of strict mode

g4_tau_info=0

# Read settings.yaml to get all the global parameters
# Open the file and iterate over the lines
# For each line, if it starts with energy, set the energy variable (everything after ": "), etc. The variables to set are energy, genie_n_events, run, pdg, start_step
echo "Reading settings.yaml..."
for line in $(cat settings.yaml); do
  if [[ $line == energy* ]]; then
    energy=$(echo $line | cut -d' ' -f2)
  elif [[ $line == genie_n_events* ]]; then
    genie_n_events=$(echo $line | cut -d' ' -f2)
  elif [[ $line == run* ]]; then
    run=$(echo $line | cut -d' ' -f2)
  elif [[ $line == pdg* ]]; then
    pdg=$(echo $line | cut -d' ' -f2)
  elif [[ $line == start_step* ]]; then
    start_step=$(echo $line | cut -d' ' -f2)
  elif [[ $line == "g4_tau_info: on" ]]; then
    g4_tau_info=1
  fi
done


# If start_step < 2, run the GENIE event generation
if [ $start_step -lt 2 ]; then  
    echo "Running GENIE event generation..."
    # Run gevgen with the energy as the incoming neutrino energy, the number of events, the run number, and the pdg code as the incoming neutrino PDG
    gevgen -n $genie_n_events -p $pdg -e $energy -r $run --cross-sections ../data/gxspl-NUsmall.xml --seed 2024 -t "1000080160[0.9],1000010010[0.1]" --tune G21_11a_00_000 > "genie_run_e${energy}.log"
    # Move the output file, named gntp.${run}.ghep.root, to the ../data directory
    mv gntp.${run}.ghep.root ../data
fi

# If start_step < 3, run the GENIE event generation
if [ $start_step -lt 3 ]; then
    echo "Converting ghep.root file to csv..."
    # Run get_pol_genie.C to generate a csv file.
    # If the number of events is less than 500 000, run the script once
    # If the number of events is greater than 500 000, run the script in a loop, processing 500 000 events at a time
    if [ $genie_n_events -lt 500000 ]; then
    genie -l -b -q 'get_pol_genie.C("../data/gntp.'$run'.ghep.root","../data/genie_tau_pol_data_e'$energy'.csv",0,'$genie_n_events')'
    else
        for i in {0..$genie_n_events..500000}; do
            genie -l -b -q 'get_pol_genie.C("../data/gntp.'$run'.ghep.root","../data/genie_tau_pol_data_e'$energy'_'$i'.csv",'$i',$(($i+500000)))'
        done
        # Merge all csv files, dropping the first line of all but the first file
        cat ../data/genie_tau_pol_data_e${energy}_0.csv > ../data/genie_tau_pol_data_e${energy}.csv
        for i in {500000..$genie_n_events..500000}; do
            tail -n +2 ../data/genie_tau_pol_data_e${energy}_${i}.csv >> ../data/genie_tau_pol_data_e${energy}.csv
        done
    fi
fi

if [ $start_step -lt 4 ]; then
    echo "Converting GENIE csv file to dat file..."
    # Run python script to convert the file into a dat file, and generate a csv file with only the tau leptons 
    python convert_genie.py -i ../data/genie_tau_pol_data_e${energy}.csv -od ../data/tauola_output_genie_e$energy.dat -oc ../data/genie_pol_e$energy.csv --write-nevents
fi

if [ $start_step -lt 5 ]; then
    echo "Running Tauola tau decay simulation with polarization..."
    # Run the Tauola tau decay simulation, with polarization
    ./decay.o ../data/tauola_input_genie_e$energy.dat ../data/tauola_output_genie_e$energy.dat 6 7 8 ../data/genie_pol_e$energy.csv &> ../logfiles/tauola_run_e$energy.log

    echo "Running Tauola tau decay simulation without polarization..."
    # Run the Tauola tau decay simulation, without polarization
    ./decay.o ../data/tauola_input_genie_e$energy.dat ../data/tauola_output_genie_nopol_e$energy.dat 0 0 0 &> ../logfiles/tauola_run_e${energy}_nopol.log
fi

# Extract the number of tau events from the settings.yaml file
# Open the file and iterate over the lines
# For each line, if it starts with tau_n_events, set the tau_events variable (everything after ": ")
for line in $(cat settings.yaml); do
  if [[ $line == tau_n_events* ]]; then
    tau_events=$(echo $line | cut -d' ' -f2)
  fi
done

# Activate the g4 environment to run the Geant4 simulation
# conda activate g4
source ~/Downloads/Program/geant4-v11.2.0-install/bin/geant4.sh

# Reaplce the last line of run_decays.mac to /run/beamOn $tau_events in run_decays.mac
sed -i '$s/.*/\/run\/beamOn '$tau_events'/' ../learn-geant4/ice-11-2/build/run_decays.mac

echo "Running Geant4 tau decay simulation..."

echo "event_num,pdg,E,px,py,pz" > ../data/geant4_output_e$energy.csv

# Reset tau info file
if [ $g4_tau_info -eq 1 ]; then
  echo "event_num,pdg,E,px,py,pz,x,y,z,status" > ../data/geant4_tau_output_e$energy.csv
fi

# Run Geant4 tau decay simulation
../learn-geant4/ice-11-2/build/exampleB1 ../learn-geant4/ice-11-2/build/run_decays.mac &> ../logfiles/geant4_run_e$energy.log

# Finish
echo "All done!"
