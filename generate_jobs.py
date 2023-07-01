import os
import sys
from os import path, makedirs
import shutil
import subprocess


def generate_submit_jobs(script_dir, directory, walltime, n_threads):

    working_folder = os.path.join(script_dir, directory)
    submit_jobs_path = os.path.join(script_dir, directory, "submit_jobs.pbs")
    mem = 4*n_threads
    # Generate the submit_jobs.pbs script content
    script_content = """#!/usr/bin/env bash
#SBATCH -J dilepton_event_0
#SBATCH -N 1
#SBATCH -n {2:d}
#SBATCH --mem={3:d}G
#SBATCH -e test.err
#SBATCH -o test.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lipei.du@mail.mcgill.ca
#SBATCH --account=def-gale
#SBATCH -t {1:s}
#SBARCH -D {0:s}

module load nixpkgs/16.09  intel/2018.3  impi/2018.3.222
module load gsl
module load hdf5-mpi/1.8.18
module load StdEnv/2020 intel/2020.1.217 hdf5/1.12.1
module load gcc

export OMP_NUM_THREADS={2:d}
./dilepton_emission.e
    """.format(working_folder, walltime, n_threads, mem)

    # Write the script content to a file
    with open('submit_jobs.pbs', 'w') as script_file:
        script_file.write(script_content)

def create_directory(folder_name, beam_energy, centrality_low, centrality_high):
    # Create directory name
    dir_name = f"{folder_name}{beam_energy}_{centrality_low}_{centrality_high}"

    # Check if the directory already exists
    if os.path.exists(dir_name):
        print(f"Warning: Directory '{dir_name}' already exists.")

        # Prompt user for action
        user_input = input("Do you want to delete the existing directory? (y/n): ")

        if user_input.lower() == 'y':
            # Delete the existing directory
            shutil.rmtree(dir_name)
        else:
            print("Exiting the program.")
            return None

    # Create the directory
    os.makedirs(dir_name, exist_ok=True)

    # Change to the new directory
    os.chdir(dir_name)

    # Create the results directory
    os.makedirs('results', exist_ok=True)
    #os.makedirs('model_parameters', exist_ok=True)

    # Return the created directory name
    return dir_name


if __name__ == '__main__':
    # Check if all required arguments are provided
    if len(sys.argv) != 6:
        print("Usage: python script.py folder_name beam_energy centrality_low centrality_high")
        sys.exit(1)

    # Retrieve the command line arguments
    folder_name = sys.argv[1]
    beam_energy = sys.argv[2]
    centrality_low = sys.argv[3]
    centrality_high = sys.argv[4]
    par_dict = sys.argv[5] # parameters_dict_user.py

    par_diretory = path.dirname(path.abspath(par_dict))
    sys.path.insert(0, par_diretory)
    parameter_dict = __import__(par_dict.split('.py')[0].split("/")[-1])

    # Retrieve the absolute directory path of the script, dilepton/
    script_dir = os.path.abspath(os.path.dirname(__file__))

    # Retrieve the parent directory path, lpdu/
    parent_directory = os.path.dirname(script_dir)

    # Call the function to create the directory structure and get the directory name, e.g. AuAu_19.6_0_10
    dir_name = create_directory(folder_name, beam_energy, centrality_low, centrality_high)

    event_dir = os.path.join(script_dir, dir_name)

    # copy parameter python files
    par_dict_path = os.path.join(script_dir, par_dict)
    dir_path = os.path.join(script_dir, dir_name)
    shutil.copy2(par_dict_path, dir_path)

    subprocess.call("(cd {}/config; ".format(script_dir)  + "python3 parameters_dict_master.py "
        + "-path {} -par {};)".format(dir_path, path.abspath(par_dict)), shell=True)

    # Create symbolic links
    os.symlink('{}/dilepton_emission.e'.format(script_dir), '{}/dilepton_emission.e'.format(event_dir))
    #os.symlink('{}/parameters.dat'.format(script_dir), '{}/parameters.dat'.format(event_dir))
    os.symlink('{}/ph_rates'.format(script_dir), '{}/ph_rates'.format(event_dir))

    os.symlink("{0:s}/iEBE-sampler/{1:s}/event_0/EVENT_RESULTS_MCGlb{1:s}_0/hydro_results_MCGlb{1:s}_0/evolution_all_xyeta.dat".format(parent_directory,dir_name), 
        "{}/results/evolution_all_xyeta.dat".format(event_dir))
    os.symlink("{0:s}/iEBE-sampler/{1:s}/event_0/EVENT_RESULTS_MCGlb{1:s}_0/hydro_results_MCGlb{1:s}_0/music_input".format(parent_directory,dir_name), 
        "{}/results/music_input".format(event_dir))

    if "walltime" in parameter_dict.control_dict.keys():
        walltime = parameter_dict.control_dict["walltime"]

    n_threads = parameter_dict.control_dict["n_threads"]

    # Call the function to generate submit_jobs.pbs inside the created directory
    generate_submit_jobs(script_dir, dir_name, walltime, n_threads)



