#!/bin/bash

# Check if the collision_energy is provided as a command line argument
if [ -z "$1" ]; then
    echo "Usage: $0 <collision_energy> <suffix>"
    exit 1
fi

# Get the collision_energy from the first command line argument
collision_energy="$1"
suffix="$2"

# List of folders
folders=("AuAu${collision_energy}_00_10" "AuAu${collision_energy}_10_20" "AuAu${collision_energy}_20_30" "AuAu${collision_energy}_30_40" "AuAu${collision_energy}_40_50" "AuAu${collision_energy}_50_60" "AuAu${collision_energy}_60_70" "AuAu${collision_energy}_70_80")

# Loop over each folder
for folder in "${folders[@]}"
do
    scp lpdu318@beluga.computecanada.ca:/lustre03/project/6002853/lpdu318/dilepton/${folder}_${suffix}/results/dilepton_results_${folder}_${suffix}.h5 /Users/Lipei/Downloads/dilepton_runs/BES_Beluga_run4
done

