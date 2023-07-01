#!/usr/bin/env python3
"""
    This script contains all the default parameters in the code.
"""

from os import path, makedirs
import sys
import shutil
import argparse

# control parameters
control_dict = {
    'walltime': "10:00:00",         # walltime to run
    'n_threads': 10,                # number of threads
}

dilepton_dict = {
    'hydro_flag': 2,                # read in mode for hydro medium
                                    # 0: read in hdf5 file
                                    # 1: read in binary file output from MUSIC
                                    # 2: read in binary file output from new MUSIC (no grid)
                                    # 3: read in binary file output from new MUSIC (on grid)
    'hydro_nskip_tau': 1,           # read in hydro slice every hydro_nskip_tau 
                                    # steps from the medium file
                                    # (only works for hydro_flag = 1)


    ##### emission rate #####

    'dilepton_emission_rate': 2,    # 0: analytical emission rate; PRC. 93, 044902, 2016
                                    # 1: LO rate; (2.4) of 1910.09567
                                    # 2: resummed NLO; from GJ, extrapolation will NOT work properly for large alpha_s, large muB, 
                                    # or if q falls outside the grid! In these cases, just the boundary value is used.

    'alpha_s': 0.3,                 # strong coupling, works for dilepton_emission_rate = 2


    ##### hydro profile #####

    'Xmin': -15.0,                  # minimum points along x direction
    'dx': 0.3,                      # lattice spacing along x direction
    'Ymin': -15.0,                  # minimum points along y direction
    'dy': 0.3,                      # lattice spacing along y direction

    'neta': 10,                     # number of points in eta direction
    'eta_i': 0.0,                   # beginning value of eta slice
    'eta_f': 3.0,                   # end value of eta slice

    'ETAmax': 2.0,                  # fluid cells outside this boundary are not considered

    'T_sw_high': 0.140,             # high end of the switching temperature, above which is full QGP
    'T_sw_low': 0.135,              # low end of the switching temperature
    'T_dec': 0.105,                 # freeze out temperature (GeV), below which is out of interest


    ##### hydro settings #####

    'HydroinfoVisflag': 1,          # determine whether to read in the viscous evolution information
    'HydroinfoBuffersize': 500,     # set the buffer size for hydro evolution profile

    'include_baryondiff_deltaf': 0, # switch to include baryon diffusion corrections
    'include_shearvisc_deltaf': 0,  # switch to include shear viscous corrections
    'turn_off_transverse_flow': 0,  # flag to turn off transverse flow in the photon calculation
    'turn_on_muB': 1,               # flag to include muB dependence in photon rates

    'test_code_flag': 0,            # flag to test the code with constant T, muB etc.
    'T_test': 0.2,
    'muB_test': 0.3,
    'rhoB_eplusp_test': 2.0,
    'inv_eplusp_test': 1.0,


    ##### photon kinematics #####

    'np': 20,                       # number of points for photon momentum
    'nphi': 20,                     # number of points for angles of photons momenta
    'nrapidity': 1,                 # number of points for photon rapidity, odd number
    'nm': 3,                        # number of points for dilepton invariant mass

    'photon_q_i': 0.001,            # the smallest photon momentum to be calculated
    'photon_q_f': 4.0,              # the largest photon momentum to be calculated
    'photon_phi_q_i': 0.0,          # the smallest angle of photon momentum
    'photon_phi_q_f': 6.2831853,    # the largest angle of photon momentum
    'photon_y_i': 0.0,              # the smallest photon rapidity
    'photon_y_f': 0.0,              # the largest photon rapidity
    'dilepton_mass_i': 0.1,         # the smallest dilepton invariant mass
    'dilepton_mass_f': 3.6,         # the largest dilepton invariant mass

    'norder': 3,                    # calculate photon vn to norder

    ##### differential setting #####

    'differential_flag': 1,         # determine whether to output differential photon yield and vn
                                    # 1: differential in T and tau
                                    # 2: differential in x and tau
                                    # 10: differeitial in all options above
    'tau_start': 1.0,               # emission start time (fm)
    'tau_end': 13.0,                # emission end time (fm)
    'n_tau_cut': 50,                # number of points in tau (range of tau is specified by tau_start and tau_end)
    'T_cuthigh': 0.25,              # maximum allowed emission T (GeV)
    'T_cutlow': 0.1,                # minimum allowed emission T (GeV)
    'nTcut': 30,                    # number of points in T (range of T is specified by T_cuthigh and T_cutlow)
    'dTau': 0.1,                    # lattice spacing along tau direction

    'calHGIdFlag': 0,               # Flag to decide whether to calculate individual HG channels

}

Parameters_list = [
    (dilepton_dict, "parameters.dat", 1)
]

path_list = [
    'model_parameters/'
]


def update_parameters_dict(par_dict_path):
    """This function update the parameters dictionaries with user's settings"""
    par_diretory = path.dirname(par_dict_path)
    sys.path.insert(0, par_diretory)
    print("\U0001F375  User's parameter file found in {}".format(par_diretory))

    parameters_dict = __import__(par_dict_path.split('.py')[0].split('/')[-1])
    
    dilepton_dict.update(parameters_dict.dilepton_dict)


def output_parameters_to_files(workfolder="."):
    """This function outputs parameters in dictionaries to files"""
    workfolder = path.abspath(workfolder)
    
    for idict, (parameters_dict, fname, itype) in enumerate(Parameters_list):
        #output_folder = path.join(workfolder, path_list[idict])
        output_folder = path.join(workfolder)
        if not path.exists(output_folder):
            makedirs(output_folder)
        f = open(path.join(output_folder, fname), "w")
        for key_name in parameters_dict:
            if itype == 1:
                f.write("{parameter_name} = {parameter_value}\n".format(
                    parameter_name=key_name,
                    parameter_value=parameters_dict[key_name]))
        f.close()

        print("\U0001F375  Output input parameter files to {}...".format(
                                                                output_folder))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='\U0000269B Welcome to dilepton parameter master',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-path', '--path', metavar='',
                        type=str, default='.',
                        help='output folder path')
    parser.add_argument('-par', '--par_dict', metavar='',
                        type=str, default='parameters_dict_user',
                        help='user-defined parameter dictionary filename')

    args = parser.parse_args()

    update_parameters_dict(path.abspath(args.par_dict))
    output_parameters_to_files(args.path)


