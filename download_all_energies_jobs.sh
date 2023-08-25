#!/bin/bash

# List of beam energies
beam_energies=('7.7' '19.6' '27') 
# '39' '54.4' '62.4' '130' '200')  # Add more beam energies as needed

# Suffix
suffix="MidrapMuoff"

# Loop over each beam energy
for energy in "${beam_energies[@]}"
do
    echo "Running for beam energy: $energy"
    sh download_all_centralities_jobs.sh "$energy" "$suffix"
done


