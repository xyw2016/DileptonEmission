
from multiprocessing import Pool
from subprocess import call
from os import path, mkdir, remove, makedirs
from glob import glob
import sys
import time
import shutil
import re
import h5py
import numpy as np
import fnmatch

def check_an_event_is_good(event_folder, required_files_list):
    """This function checks the given event contains all required files"""
    event_file_list = glob(path.join(event_folder, "*"))
    for pattern in required_files_list:
        matching_files = glob(path.join(event_folder, pattern))
        if not matching_files:
            print("Event {} is bad, missing {} ...".format(
                event_folder, pattern),
                  flush=True)
            return False
    return True

def zip_spvn_results_into_hdf5(final_results_folder, event_id):
    """This function combines all the spvn results into hdf5"""
    results_name = "dilepton_results_{}".format(event_id)

    required_files_list = [
        'Total_dilepton_*.dat',
        'QGP_*.dat'
    ]

    status = check_an_event_is_good(final_results_folder, required_files_list)

    if status:
        curr_time = time.asctime()
        print("[{}] {} is good, converting results to hdf5".format(
            curr_time, final_results_folder),
              flush=True)

        hf = h5py.File("{0}.h5".format(results_name), "w")
        gtemp = hf.create_group("{0}".format(results_name))
        file_list = glob(path.join(final_results_folder, "*"))

        for file_path in file_list:
            file_name = file_path.split("/")[-1]
            
            # Check if the file name matches any pattern in required_files_list
            if any(fnmatch.fnmatch(file_name, pattern) for pattern in required_files_list):
                dtemp = np.loadtxt(file_path)
                h5data = gtemp.create_dataset("{0}".format(file_name),
                                              data=dtemp,
                                              compression="gzip",
                                              compression_opts=9)
                # save header
                ftemp = open(file_path, "r")
                header_text = str(ftemp.readline())
                ftemp.close()
                if header_text.startswith("#"):
                    h5data.attrs.create("header", np.string_(header_text))

                # Remove the original file
                remove(file_path)

        hf.close()
        shutil.move("{}.h5".format(results_name), final_results_folder)

    else:
        print("{} is broken, skipped".format(final_results_folder), flush=True)
    return (status)


def main(event_id):
    final_results_folder = "results"
    zip_spvn_results_into_hdf5(final_results_folder, event_id)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py event_id")
    else:
        event_id = sys.argv[1]
        main(event_id)
