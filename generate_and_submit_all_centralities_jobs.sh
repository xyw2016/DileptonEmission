#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <beam_energy> <suffix> <parameter_file>"
    exit 1
fi

beam_energy="$1"
suffix="$2"
parameter_file="$3"

# Get the current directory
current_dir=$(pwd)

# Loop over centrality bins from 00 to 70 with step 10
for ((centrality=0; centrality<=70; centrality+=10))
do
    # Format the centrality bin with leading zeros if needed
    formatted_centrality=$(printf "%02d" "$centrality")

    # Execute the generate_jobs.py command
    command="python $current_dir/generate_jobs.py AuAu ${beam_energy} ${formatted_centrality} $((centrality+10)) ${suffix} ${parameter_file}"
    echo "Running command: $command"
    eval "$command"

    # Submit the job using sbatch
    folder_name="AuAu${beam_energy}_${formatted_centrality}_$((centrality+10))_${suffix}"
    cd "$folder_name"
    sbatch submit_jobs.pbs

    # Return to the original directory where generate_jobs.py is located
    cd "$current_dir"
done
