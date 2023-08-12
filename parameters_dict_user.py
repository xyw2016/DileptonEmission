#!/usr/bin/env python3
"""
    This script contains all the default parameters in the iEBE-MUSIC package.
"""

# control parameters
control_dict = {
    'walltime': "3:00:00",         # walltime to run
    'n_threads': 40,                # number of threads
    'n_memory_to_thread': 10,       # memory per thread
}

dilepton_dict = {
    ##### emission rate #####

    'dilepton_emission_rate': 2,    # 0: analytical emission rate; PRC. 93, 044902, 2016
                                    # 1: LO rate; (2.4) of 1910.09567
                                    # 2: resummed NLO; from GJ, extrapolation will NOT work properly for large alpha_s, large muB, 
                                    # or if q falls outside the grid! In these cases, just the boundary value is used.

    'alpha_s': 0.3,                 # strong coupling, works for dilepton_emission_rate = 2


    ##### hydro profile #####

    'ETAmax': 10.0,                  # fluid cells outside this boundary are not considered

    'T_sw_high': 0.105,             # high end of the switching temperature, above which is full QGP
    'T_sw_low': 0.103,              # low end of the switching temperature
    'T_dec': 0.120,                 # freeze out temperature (GeV), below which is out of interest


    ##### hydro settings #####

    'include_baryondiff_deltaf': 1, # switch to include baryon diffusion corrections
    'include_shearvisc_deltaf': 1,  # switch to include shear viscous corrections
    'turn_off_transverse_flow': 0,  # flag to turn off transverse flow in the photon calculation
    'turn_on_muB': 1,               # flag to include muB dependence in photon rates

    'test_code_flag': 0,            # flag to test the code with constant T, muB etc.
    'T_test': 0.2,
    'muB_test': 0.3,
    'rhoB_eplusp_test': 2.0,
    'inv_eplusp_test': 1.0,


    ##### photon kinematics #####

    'np': 20,                       # number of points for photon momentum
    'nphi': 40,                     # number of points for angles of photons momenta
    'nrapidity': 7,                 # number of points for photon rapidity, odd number
    'nm': 10,                       # number of points for dilepton invariant mass

    'photon_q_i': 0.0,              # the smallest photon momentum to be calculated
    'photon_q_f': 6.0,              # the largest photon momentum to be calculated
    'photon_phi_q_i': 0.0,          # the smallest angle of photon momentum
    'photon_phi_q_f': 6.2831853,    # the largest angle of photon momentum
    'photon_y_i': -1.0,              # the smallest photon rapidity
    'photon_y_f': 1.0,              # the largest photon rapidity
    'dilepton_mass_i': 0.1,         # the smallest dilepton invariant mass
    'dilepton_mass_f': 3.0,         # the largest dilepton invariant mass

    'use_logarithmic_mass_grid': 1, # choose between equal steps and logarithmically spaced mass grid

    'norder': 3,                    # calculate photon vn to norder

    ##### differential setting #####

    'differential_flag': 1,         # determine whether to output differential photon yield and vn
                                    # 1: differential in T and tau
                                    # 2: differential in x and tau
                                    # 10: differeitial in all options above
    'n_tau_cut': 100,                # number of points in tau (range of tau is specified by tau_start and tau_end)
    'nTcut': 50,                    # number of points in T (range of T is specified by T_cuthigh and T_cutlow)
}
