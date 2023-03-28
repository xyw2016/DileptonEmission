// Copyright 2016 Chun Shen
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <omp.h>

#include "Hydroinfo_h5.h"
#include "ThermalPhoton.h"
#include "QGP_LO.h"
#include "tensor_trans.h"
#include "PhotonEmission.h"
#include "ParameterReader.h"
#include "Arsenal.h"

using namespace std;
using ARSENAL::createA2DMatrix;
using ARSENAL::createA3DMatrix;
using ARSENAL::createA4DMatrix;
using ARSENAL::createA5DMatrix;
using ARSENAL::createA6DMatrix;
using ARSENAL::deleteA2DMatrix;
using ARSENAL::deleteA3DMatrix;
using ARSENAL::deleteA4DMatrix;
using ARSENAL::deleteA5DMatrix;
using ARSENAL::deleteA6DMatrix;
using TENSORTRANSFORM::getTransverseflow_u_mu_low;
using TENSORTRANSFORM::lorentz_boost_matrix;

using PhysConsts::hbarC;
using PhysConsts::eps;

#define CODE_TEST 0

PhotonEmission::PhotonEmission(std::shared_ptr<ParameterReader> paraRdr_in) {
    paraRdr = paraRdr_in;
    output_path = "results/";

    hydro_flag = paraRdr->getVal("hydro_flag");
    differential_flag = paraRdr->getVal("differential_flag");
    turn_off_transverse_flow = paraRdr->getVal("turn_off_transverse_flow");
    turn_on_muB_ = static_cast<int>(paraRdr->getVal("turn_on_muB", 1));

    diff_flag = paraRdr->getVal("differential_flag");

    // omp parameters
    CORES = 1;

#ifdef _OPENMP
    CORES = omp_get_max_threads();
#endif

    set_hydroGridinfo();
    print_hydroGridinfo();

    // read the photon emission rate tables
    InitializePhotonEmissionRateTables();

    lambda = createA2DMatrix(4, 4, 0.);

    dNd2pTdphidy_eq = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidy_tot = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTd2M_eq = createA2DMatrix(nm, np, 0);
    dNd2pTd2M_tot = createA2DMatrix(nm, np, 0);
    dNd2Mdy_eq.resize(nm, 0);
    dNd2Mdy_tot.resize(nm, 0);

    vnpT_cos_eq = createA3DMatrix(norder, nm, np, 0.);
    vnpT_sin_eq = createA3DMatrix(norder, nm, np, 0.);
    vnpT_cos_tot = createA3DMatrix(norder, nm, np, 0.);
    vnpT_sin_tot = createA3DMatrix(norder, nm, np, 0.);

    vn_cos_eq = createA2DMatrix(norder, nm, 0.);
    vn_sin_eq = createA2DMatrix(norder, nm, 0.);
    vn_cos_tot = createA2DMatrix(norder, nm, 0.);
    vn_sin_tot = createA2DMatrix(norder, nm, 0.);

    if (diff_flag == 1) {
        dNd2pTdphidydTdtau_eq = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_diff = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_tot = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);

        dNd2pTdphidydTdtau_eq_all = createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
        dNd2pTdphidydTdtau_diff_all = createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
        dNd2pTdphidydTdtau_tot_all = createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
    }
}


PhotonEmission::~PhotonEmission() {
    deleteA2DMatrix(lambda, 4);

    deleteA4DMatrix(dNd2pTdphidy_eq, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_tot, nm, np, nphi);
    deleteA2DMatrix(dNd2pTd2M_eq, nm);
    deleteA2DMatrix(dNd2pTd2M_tot, nm);

    deleteA3DMatrix(vnpT_cos_eq, norder, nm);
    deleteA3DMatrix(vnpT_sin_eq, norder, nm);
    deleteA3DMatrix(vnpT_cos_tot, norder, nm);
    deleteA3DMatrix(vnpT_sin_tot, norder, nm);

    deleteA2DMatrix(vn_cos_eq, norder);
    deleteA2DMatrix(vn_sin_eq, norder);
    deleteA2DMatrix(vn_cos_tot, norder);
    deleteA2DMatrix(vn_sin_tot, norder);

    if (diff_flag == 1) {
        nTcut = paraRdr->getVal("nTcut");
        n_tau_cut = paraRdr->getVal("n_tau_cut");
        deleteA6DMatrix(dNd2pTdphidydTdtau_eq, nTcut, n_tau_cut, nm, np, nphi);
        deleteA6DMatrix(dNd2pTdphidydTdtau_diff, nTcut, n_tau_cut, nm, np, nphi);
        deleteA6DMatrix(dNd2pTdphidydTdtau_tot, nTcut, n_tau_cut, nm, np, nphi);

        deleteA3DMatrix(dNd2pTdphidydTdtau_eq_all, nTcut, n_tau_cut);
        deleteA3DMatrix(dNd2pTdphidydTdtau_diff_all, nTcut, n_tau_cut);
        deleteA3DMatrix(dNd2pTdphidydTdtau_tot_all, nTcut, n_tau_cut);
    }
}


void PhotonEmission::set_hydroGridinfo() {
    gridX0 = paraRdr->getVal("Xmin");
    gridY0 = paraRdr->getVal("Ymin");
    gridDx = paraRdr->getVal("dx");
    gridDy = paraRdr->getVal("dy");
    gridDtau = paraRdr->getVal("dTau");

    nm = paraRdr->getVal("nm");
    neta = paraRdr->getVal("neta");
    np = paraRdr->getVal("np");
    nphi = paraRdr->getVal("nphi");
    nrapidity = paraRdr->getVal("nrapidity");
    norder = paraRdr->getVal("norder");

    T_dec = paraRdr->getVal("T_dec");
    T_sw_high = paraRdr->getVal("T_sw_high");
    T_sw_low = paraRdr->getVal("T_sw_low");

    T_cuthigh = paraRdr->getVal("T_cuthigh");
    T_cutlow = paraRdr->getVal("T_cutlow");
    nTcut = paraRdr->getVal("nTcut");
    tau_cut_high = paraRdr->getVal("tau_end");
    tau_cut_low = paraRdr->getVal("tau_start");
    n_tau_cut = paraRdr->getVal("n_tau_cut");

    calHGIdFlag = paraRdr->getVal("CalHGIdFlag");
}


void PhotonEmission::print_hydroGridinfo() {
    cout << "----------------------------------------" << endl;
    cout << "-- Parameters list for photon emission:" << endl;
    cout << "----------------------------------------" << endl;
    cout << "tau_start =" << paraRdr->getVal("tau_start") << " fm/c." << endl;
    cout << "tau_end =" << paraRdr->getVal("tau_end") << " fm/c." << endl;
    cout << "dTau = " << gridDtau << " fm/c" << endl;
    cout << "X_min = " << gridX0 << " fm/c" << endl;
    cout << "dx = " << gridDx << " fm/c" << endl;
    cout << "Y_min = " << gridY0 << " fm/c" << endl;
    cout << "dy = " << gridDy << " fm/c" << endl;
    cout << endl;

    cout << "T_dec = " << T_dec << " GeV." << endl;
    cout << "T_sw = " << T_sw_low <<  " to " << T_sw_high << " GeV."<< endl;
    cout << endl;

    cout << "Photon momentum: " << paraRdr->getVal("photon_q_i")
         << " to " << paraRdr->getVal("photon_q_f") << " GeV, "
         << "n_q =" << np << endl;
    cout << "Photon momentum angles: " << paraRdr->getVal("photon_phi_q_i")
         << " to " << paraRdr->getVal("photon_phi_q_f")
         << ", n_phi=" << nphi << endl;
    cout << "Photon momentum rapidity: " << paraRdr->getVal("photon_y_i")
         << " to " << paraRdr->getVal("photon_y_f")
         << ", n_y =" << nrapidity << endl;
    cout << "Dilepton invariant mass: " << paraRdr->getVal("dilepton_mass_i")
         << " to " << paraRdr->getVal("dilepton_mass_f") << " GeV, "
         << "n_m =" << nm << endl;
    cout << "Calculate individual channels in Hadron Resonance Gas phase: ";
    if (calHGIdFlag == 0) {
       cout << " No! " << endl;
    } else {
       cout << " Yes!" << endl;
    }

    
}


void PhotonEmission::InitializePhotonEmissionRateTables() {
    double photonrate_tb_Emin = paraRdr->getVal("PhotonemRatetableInfo_Emin");
    double photonrate_tb_Tmin = paraRdr->getVal("PhotonemRatetableInfo_Tmin");
    double photonrate_tb_dE = paraRdr->getVal("PhotonemRatetableInfo_dE");
    double photonrate_tb_dT = paraRdr->getVal("PhotonemRatetableInfo_dT");

    dilepton_QGP_LO = std::unique_ptr<ThermalPhoton>(
            new QGP_LO(paraRdr, "QGP_LO_total"));
    // dilepton_QGP_LO->setupEmissionrateFromFile(
    //         photonrate_tb_Tmin, photonrate_tb_dT,
    //         photonrate_tb_Emin, photonrate_tb_dE, true, true);// true for vis corrections
}


void PhotonEmission::calPhotonemission_3d(void *hydroinfo_ptr_in) {
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;
    hydroinfo_MUSIC_ptr =
                reinterpret_cast<Hydroinfo_MUSIC*>(hydroinfo_ptr_in);

    // photon momentum in the lab frame
    double M_ll[nm]; // invariant mass array
    double p_q[np], phi_q[nphi], y_q[nrapidity];
    double sin_phiq[nphi], cos_phiq[nphi];
    double p_lab_local[4];
    double p_lab_Min[4]; // dilepton 4-momentum in Minkowski coordinates
    for (int k = 0; k < nrapidity; k++) {
        y_q[k] = dilepton_QGP_LO->getPhotonrapidity(k);
    }
    for (int l = 0; l < np; l++) {
        p_q[l] = dilepton_QGP_LO->getPhotonp(l);
    }
    for (int m = 0; m < nphi; m++) {
        phi_q[m] = dilepton_QGP_LO->getPhotonphi(m);
        sin_phiq[m] = sin(phi_q[m]);
        cos_phiq[m] = cos(phi_q[m]);
    }
    for (int j = 0; j < nm; j++) {
        M_ll[j] = dilepton_QGP_LO->getDileptonMass(j);
    }

    // get hydro grid information
    double tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
    double dtau = hydroinfo_MUSIC_ptr->get_hydro_dtau();
    double dx = hydroinfo_MUSIC_ptr->get_hydro_dx();
    double deta = hydroinfo_MUSIC_ptr->get_hydro_deta();
    double eta_max = hydroinfo_MUSIC_ptr->get_hydro_eta_max();
    double X_max = hydroinfo_MUSIC_ptr->get_hydro_x_max();
    double Nskip_x = hydroinfo_MUSIC_ptr->get_hydro_Nskip_x();
    double Nskip_eta = hydroinfo_MUSIC_ptr->get_hydro_Nskip_eta();
    double Nskip_tau = hydroinfo_MUSIC_ptr->get_hydro_Nskip_tau();
    double volume_base = Nskip_tau*dtau*Nskip_x*dx*Nskip_x*dx*Nskip_eta*deta;

    // output results in (T, tau)
    double dT_cut = (T_cuthigh - T_cutlow)/(nTcut - 1);
    double dtau_cut = (tau_cut_high - tau_cut_low)/(n_tau_cut - 1);

    double tau_now = 0.0;
    double flow_u_mu_low[4];
    double flow_u_mu_Min[4]; // fluid velocity in Minkowski
    fluidCell_3D_new fluidCellptr;
    double **lambda_munu = createA2DMatrix(4, 4, 0.);

    // number of fluid cells
    long int number_of_cells =(
                hydroinfo_MUSIC_ptr->get_number_of_fluid_cells_3d());
    cout << "number of cells:" << number_of_cells << endl;

    // multi-threads setup
    long FO_chunk = number_of_cells / CORES;
    long remainder = number_of_cells  -  CORES * FO_chunk;

    cout << "Number of cores : " << CORES << endl;
    cout << "Chunk size = " << FO_chunk << endl;
    cout << "Remainder cells = " << remainder << endl;

    if(remainder != 0) FO_chunk++;

    cout << "----------------------------------------" << endl;

    // arrays to store values across cores
    double *dNd2pTdphidy_eq_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
    double *dNd2pTdphidy_diff_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
    double *dNd2pTdphidy_tot_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));

    // main loop begins ...
    // loop over all fluid cells
    // subdivide bite size chunks of freezeout surface across cores
    int ncells = 0;
    #pragma omp parallel for //reduction(+:ncells) reduction(add_dNd2pTdphidydTdtau_eq:dNd2pTdphidydTdtau_eq_all)
    for(long n = 0; n < CORES; n++)
    {
        long endFO = FO_chunk;

        for(long int icell = 0; icell < endFO; icell++)  // cell index inside each chunk
        {
            if((icell == endFO - 1) && (remainder != 0) && (n > remainder - 1)) continue;

            long cell_id = n  +  icell * CORES;

            hydroinfo_MUSIC_ptr->get_hydro_cell_info_3d(cell_id, fluidCellptr);

            double tau_local = tau0 + fluidCellptr.itau*dtau;
            double eta_local = -eta_max + fluidCellptr.ieta*deta;

#ifndef _OPENMP
            if (fabs(tau_now - tau_local) > 1e-10) {
                tau_now = tau_local;
                cout << "Calculating tau = " << setw(4) << setprecision(3)
                     << tau_now << " fm/c..." << endl;
            }
#endif
            // volume element: tau*dtau*dx*dy*deta,
            double volume = tau_local*volume_base;

            const double ed_local = fluidCellptr.ed;
            const double pd_local = fluidCellptr.pressure;
            double temp_local = fluidCellptr.temperature;
            double muB_local = turn_on_muB_*fluidCellptr.muB;
            const double rhoB_local = turn_on_muB_*fluidCellptr.rhoB;
            double temp_inv = 1/temp_local;
            double rhoB_over_eplusp = rhoB_local/(ed_local+pd_local);

            // if muB is off, muB is set to a tiny value
            if(turn_on_muB_==0)
                muB_local = eps;

            // validation setup
            if(CODE_TEST==1){
                temp_local = 0.25;
                temp_inv = 1/temp_local;
                muB_local = 0.9;
                rhoB_over_eplusp = 8.0;
                volume = 1.0;
                eta_local = 0.0;
                turn_off_transverse_flow = 1;
            }

            // fluid cell is out of interest
            if (temp_local < T_dec || temp_local > T_cuthigh || temp_local < T_cutlow)
                continue;

            // indices for the output in (T, tau)
            int idx_T = (int)((temp_local - T_cutlow)/dT_cut + eps);
            int idx_tau = (int)((tau_local - tau_cut_low)/dtau_cut + eps);

            double ux, uy, ueta;
            if (turn_off_transverse_flow == 1) {
                ux = 0.0;
                uy = 0.0;
                ueta = 0.0;
            } else {
                ux = fluidCellptr.ux;
                uy = fluidCellptr.uy;
                ueta = fluidCellptr.ueta;
            }

            double utau = sqrt(1. + ux*ux + uy*uy + ueta*ueta);

            flow_u_mu_low[0] = utau;
            flow_u_mu_low[1] = -ux;
            flow_u_mu_low[2] = -uy;
            flow_u_mu_low[3] = -ueta;

            double cosh_eta = cosh(eta_local);
            double sinh_eta = sinh(eta_local);

            // 4 velocity in Minkowski in lab frame

            flow_u_mu_Min[0] = cosh_eta*utau + sinh_eta*ueta;
            flow_u_mu_Min[1] = ux;
            flow_u_mu_Min[2] = uy;
            flow_u_mu_Min[3] = sinh_eta*utau + cosh_eta*ueta;

            // Lorentz boost matrix from lab to local rest frame

            lorentz_boost_matrix(lambda_munu, flow_u_mu_Min[0], 
                flow_u_mu_Min[1], flow_u_mu_Min[2], flow_u_mu_Min[3]);

            double pi11 = fluidCellptr.pi11;
            double pi12 = fluidCellptr.pi12;
            double pi13 = fluidCellptr.pi13;
            double pi22 = fluidCellptr.pi22;
            double pi23 = fluidCellptr.pi23;
            // reconstruct all other components of the shear stress tensor
            double pi01 = (ux*pi11 + uy*pi12 + ueta*pi13)/utau;
            double pi02 = (ux*pi12 + uy*pi22 + ueta*pi23)/utau;
            double pi33 = (utau*(ux*pi01 + uy*pi02) - utau*utau*(pi11 + pi22)
                           + ueta*(ux*pi13 + uy*pi23))/(utau*utau - ueta*ueta);
            double pi00 = pi11 + pi22 + pi33;
            double pi03 = (ux*pi13 + uy*pi23 + ueta*pi33)/utau;

            double bulkPi_local = fluidCellptr.bulkPi;

            // qeta = tau*qeta, qi is qi/kappa_hat
            double qx = fluidCellptr.qx;
            double qy = fluidCellptr.qy;
            double qeta = fluidCellptr.qz;
            double qtau = (ux*qx + uy*qy + ueta*qeta)/utau;

            double prefactor_pimunu = 1./2.;
            double prefactor_diff = 1./(4.*pow(2*M_PI, 5)) * temp_inv/pow(hbarC, 4);

            // photon momentum loops
            for (int k = 0; k < nrapidity; k++) {
                double cosh_y = cosh(y_q[k]);
                double sinh_y = sinh(y_q[k]);
                double cosh_y_minus_eta = cosh(y_q[k] - eta_local);
                double sinh_y_minus_eta = sinh(y_q[k] - eta_local);
                int i0 = nphi * k;
                for (int m = 0; m < nphi; m++) {
                    int i1 = (m+i0) * np;
                    for (int l = 0; l < np; l++) {
                        int i2 = (l+i1) * nm;
                        for (int j = 0; j < nm; j++) {
                            int i3 = (j+i2);

                            double M_T = sqrt(p_q[l]*p_q[l]+M_ll[j]*M_ll[j]); // transverse mass

                            // p_q is p_T array
                            // Minkowski four momentum in lab frame
                            p_lab_Min[0] = M_T * cosh_y;
                            p_lab_Min[1] = p_q[l]*cos_phiq[m];
                            p_lab_Min[2] = p_q[l]*sin_phiq[m];
                            p_lab_Min[3] = M_T * sinh_y;

                            // from Minkowski to Milne in lab frame
                            p_lab_local[0] = M_T * cosh_y_minus_eta;
                            p_lab_local[1] = p_q[l]*cos_phiq[m];
                            p_lab_local[2] = p_q[l]*sin_phiq[m];
                            p_lab_local[3] = M_T * sinh_y_minus_eta;

                            // from Lab frame to LRF frame
                            double p_Min_lrf[4];
                            for (int j = 0; j < 4; j++) {
                                p_Min_lrf[j] = 0.;
                                for (int i = 0; i < 4; i++) {
                                    p_Min_lrf[j] += lambda_munu[j][i]*p_lab_Min[i];
                                }                            
                            }

                            double pvec_lrf, pvec3; // spatial magnitude, |vec p| and |vec p|^5
                            pvec_lrf = sqrt(p_Min_lrf[1]*p_Min_lrf[1] + p_Min_lrf[2]*p_Min_lrf[2] + p_Min_lrf[3]*p_Min_lrf[3]);
                            pvec3 = pow(pvec_lrf, 3);

                            double Eq_localrest_temp = 0.0e0;
                            for (int local_i = 0; local_i < 4; local_i++) {
                                Eq_localrest_temp += (flow_u_mu_low[local_i]
                                                      *p_lab_local[local_i]);
                            }

                            double pi_photon = (
                                  p_lab_local[0]*p_lab_local[0]*pi00
                                - 2.*p_lab_local[0]*p_lab_local[1]*pi01
                                - 2.*p_lab_local[0]*p_lab_local[2]*pi02
                                - 2.*p_lab_local[0]*p_lab_local[3]*pi03
                                + p_lab_local[1]*p_lab_local[1]*pi11
                                + 2.*p_lab_local[1]*p_lab_local[2]*pi12
                                + 2.*p_lab_local[1]*p_lab_local[3]*pi13
                                + p_lab_local[2]*p_lab_local[2]*pi22
                                + 2.*p_lab_local[2]*p_lab_local[3]*pi23
                                + p_lab_local[3]*p_lab_local[3]*pi33);

                            // dot product of diffusion and dilepton momentum, calculated in lab frame
                            // note that 1/kappa_hat included
                            double diff_dot_p = qtau*p_lab_local[0] - qx*p_lab_local[1] 
                                                - qy*p_lab_local[2] - qeta*p_lab_local[3];
                            // validation setup
                            if(CODE_TEST==1){
                                diff_dot_p = 1.0;
                            }
                            
                            double Eq_localrest_Tb = Eq_localrest_temp;
                            double pi_photon_Tb = pi_photon*prefactor_pimunu;
                            double bulkPi_Tb = bulkPi_local;
                            double diff_Tb = prefactor_diff * diff_dot_p * M_ll[j] * M_ll[j] / pvec3;

                            double dNd2pTdphidy_cell_eq = 0.;
                            double dNd2pTdphidy_cell_diff = 0.;
                            double dNd2pTdphidy_cell_tot = 0.;

                            // begin to calculate thermal photon emission
                            if (temp_local > T_sw_high) {
                                // QGP emission
                                double QGP_fraction = 1.0;
                                dilepton_QGP_LO->calThermalPhotonemission_3d(
                                    Eq_localrest_Tb, M_ll[j], pi_photon_Tb, bulkPi_Tb, diff_Tb,
                                    temp_local, muB_local, rhoB_over_eplusp, volume, QGP_fraction,
                                    dNd2pTdphidy_cell_eq, dNd2pTdphidy_cell_diff, dNd2pTdphidy_cell_tot);

                                // total fluid cells of QGP emission
                                ncells++;
                            }

                            // add contributions from QGP and Hadronic matter, etc
                            dNd2pTdphidy_eq_all[n + CORES * i3] += dNd2pTdphidy_cell_eq;
                            dNd2pTdphidy_diff_all[n + CORES * i3] += dNd2pTdphidy_cell_diff;
                            dNd2pTdphidy_tot_all[n + CORES * i3] += dNd2pTdphidy_cell_tot;


                            // These 3 lines also work for 1 thread
                            // dNd2pTdphidydTdtau_eq[idx_T][idx_tau][j][l][m][k] += dNd2pTdphidy_cell_eq;
                            // dNd2pTdphidydTdtau_diff[idx_T][idx_tau][j][l][m][k] += dNd2pTdphidy_cell_diff;
                            // dNd2pTdphidydTdtau_tot[idx_T][idx_tau][j][l][m][k] += dNd2pTdphidy_cell_tot;


                            if (differential_flag == 1) {

                                dNd2pTdphidydTdtau_eq_all[idx_T][idx_tau][n+CORES*i3] += dNd2pTdphidy_cell_eq;;
                                dNd2pTdphidydTdtau_diff_all[idx_T][idx_tau][n+CORES*i3] += dNd2pTdphidy_cell_diff;
                                dNd2pTdphidydTdtau_tot_all[idx_T][idx_tau][n+CORES*i3] += dNd2pTdphidy_cell_tot;
                            }

                        } // M_ll
                    } // p_T
                } // phi_p
            } // y
        } // fluid cell
    }// cores

    #pragma omp parallel for collapse(4)
    for (int k = 0; k < nrapidity; k++) {
        for (int m = 0; m < nphi; m++) {
            for (int l = 0; l < np; l++) {
                for (int j = 0; j < nm; j++) {

                    int i3 = (j+(l+(m+nphi * k) * np) * nm);

                    double dN_pTdpTdphidy_eq_tmp = 0.0; // reduction variable
                    double dN_pTdpTdphidy_diff_tmp = 0.0;
                    double dN_pTdpTdphidy_tot_tmp = 0.0;

                    #pragma omp simd reduction(+:dN_pTdpTdphidy_eq_tmp,dN_pTdpTdphidy_diff_tmp,dN_pTdpTdphidy_tot_tmp)
                    for(long n = 0; n < CORES; n++)
                    {
                        dN_pTdpTdphidy_eq_tmp += dNd2pTdphidy_eq_all[n+CORES*i3];
                        dN_pTdpTdphidy_diff_tmp += dNd2pTdphidy_diff_all[n+CORES*i3];
                        dN_pTdpTdphidy_tot_tmp += dNd2pTdphidy_tot_all[n+CORES*i3];
                    } // sum over the cores

                    dNd2pTdphidy_eq[j][l][m][k] = dN_pTdpTdphidy_eq_tmp;
                    //dNd2pTdphidy_diff[j][l][m][k] = dN_pTdpTdphidy_diff_tmp;
                    dNd2pTdphidy_tot[j][l][m][k] = dN_pTdpTdphidy_tot_tmp;

                    // distribution in (T, tau)
                    if (differential_flag == 1) {
                        for (int iT = 0; iT < nTcut; iT++) {
                             for (int it = 0; it < n_tau_cut; it++) {

                                double dN_pTdpTdphidydTdtau_eq_tmp = 0.0;
                                double dN_pTdpTdphidydTdtau_diff_tmp = 0.0;
                                double dN_pTdpTdphidydTdtau_tot_tmp = 0.0;

                                #pragma omp simd reduction(+:dN_pTdpTdphidydTdtau_eq_tmp,dN_pTdpTdphidydTdtau_diff_tmp,dN_pTdpTdphidydTdtau_tot_tmp)
                                for(long n = 0; n < CORES; n++) {
                                    dN_pTdpTdphidydTdtau_eq_tmp += dNd2pTdphidydTdtau_eq_all[iT][it][n+CORES*i3];
                                    dN_pTdpTdphidydTdtau_diff_tmp += dNd2pTdphidydTdtau_diff_all[iT][it][n+CORES*i3];
                                    dN_pTdpTdphidydTdtau_tot_tmp += dNd2pTdphidydTdtau_tot_all[iT][it][n+CORES*i3];
                                }

                                dNd2pTdphidydTdtau_eq[iT][it][j][l][m][k] = dN_pTdpTdphidydTdtau_eq_tmp;
                                dNd2pTdphidydTdtau_diff[iT][it][j][l][m][k] = dN_pTdpTdphidydTdtau_diff_tmp;
                                dNd2pTdphidydTdtau_tot[iT][it][j][l][m][k] = dN_pTdpTdphidydTdtau_tot_tmp;

                            }
                        }
                    }

                } // M_ll
            } // p_T
        } // phi_p
    } // y


    // Total number of elements that satisfy the condition
    int total_count = ncells;
    printf("Cells above T_sw_high=%d...\n", total_count);

    // free memory
    free(dNd2pTdphidy_eq_all);
    free(dNd2pTdphidy_diff_all);
    free(dNd2pTdphidy_tot_all);
}


void PhotonEmission::calPhoton_SpvnpT_individualchannel() {
    // dilepton_QGP_LO->calPhoton_SpvnpT_shell();
    if (differential_flag == 1) {
        dilepton_QGP_LO->calPhoton_SpMatrix_dTdtau(dNd2pTdphidydTdtau_eq, 
            dNd2pTdphidydTdtau_tot, dNd2pTdphidydTdtau_diff);
        dilepton_QGP_LO->calPhoton_SpvnpT_dTdtau();
    }
}


void PhotonEmission::outputPhotonSpvn_individualchannel() {
    // dilepton_QGP_LO->outputPhoton_SpvnpT_shell(output_path);
    if (differential_flag == 1) {
        dilepton_QGP_LO->outputPhoton_SpvnpT_dTdtau(output_path);
        dilepton_QGP_LO->outputPhoton_spectra_dTdtau(output_path);
    }
}


void PhotonEmission::calPhoton_total_Spvn() {
   // #pragma omp parallel for collapse(4)
    for (int m = 0; m < nm; m++) {
        for (int i = 0; i < np; i++) {
            double p = dilepton_QGP_LO->getPhotonp(i);
            double pweight = dilepton_QGP_LO->getPhoton_pweight(i);
            for (int j = 0; j < nphi; j++) {
                double phi = dilepton_QGP_LO->getPhotonphi(j);
                double phiweight = dilepton_QGP_LO->getPhoton_phiweight(j);
                double dy = dilepton_QGP_LO->get_dy();
                double weight = phiweight*dy;
                for (int k = 0; k < nrapidity; k++) {
                    // integrate over rapidity and azimuthal angle of momentum
                    dNd2pTd2M_eq[m][i] += dNd2pTdphidy_eq[m][i][j][k]*weight;
                    dNd2pTd2M_tot[m][i] += dNd2pTdphidy_tot[m][i][j][k]*weight;
                    for (int order = 0; order < norder; order++) {
                        vnpT_cos_eq[order][m][i] += (
                                dNd2pTdphidy_eq[m][i][j][k]*cos(order*phi)*weight);
                        vnpT_cos_tot[order][m][i] += (
                                dNd2pTdphidy_tot[m][i][j][k]*cos(order*phi)*weight);
                        vnpT_sin_eq[order][m][i] += (
                                dNd2pTdphidy_eq[m][i][j][k]*sin(order*phi)*weight);
                        vnpT_sin_tot[order][m][i] += (
                                dNd2pTdphidy_tot[m][i][j][k]*sin(order*phi)*weight);
                    }
                }
            }
            dNd2Mdy_eq[m] += dNd2pTd2M_eq[m][i]*p*pweight; // integral over p
            dNd2Mdy_tot[m] += dNd2pTd2M_tot[m][i]*p*pweight;

            for (int order=0; order < norder; order++) {
                // sum over pT
                vn_cos_eq[order][m] += vnpT_cos_eq[order][m][i]*p*pweight;
                vn_sin_eq[order][m] += vnpT_sin_eq[order][m][i]*p*pweight;
                vn_cos_tot[order][m] += vnpT_cos_tot[order][m][i]*p*pweight;
                vn_sin_tot[order][m] += vnpT_sin_tot[order][m][i]*p*pweight;

                vnpT_cos_eq[order][m][i] = vnpT_cos_eq[order][m][i]/dNd2pTd2M_eq[m][i];
                vnpT_cos_tot[order][m][i] = vnpT_cos_tot[order][m][i]/dNd2pTd2M_tot[m][i];
                vnpT_sin_eq[order][m][i] = vnpT_sin_eq[order][m][i]/dNd2pTd2M_eq[m][i];
                vnpT_sin_tot[order][m][i] = vnpT_sin_tot[order][m][i]/dNd2pTd2M_tot[m][i];
            }
            dNd2pTd2M_eq[m][i] = dNd2pTd2M_eq[m][i]/(2*M_PI); // dN/(2pi pTdpT)
            dNd2pTd2M_tot[m][i] = dNd2pTd2M_tot[m][i]/(2*M_PI);
        }
    }
    
    for (int m = 0; m < nm; m++) {
        for (int order = 1; order < norder ; order++) {
            vn_cos_eq[order][m] = vn_cos_eq[order][m]/dNd2Mdy_eq[m];
            vn_sin_eq[order][m] = vn_sin_eq[order][m]/dNd2Mdy_eq[m];
            vn_cos_tot[order][m] = vn_cos_tot[order][m]/dNd2Mdy_tot[m];
            vn_sin_tot[order][m] = vn_sin_tot[order][m]/dNd2Mdy_tot[m];
        }
    }
}


void PhotonEmission::outputPhoton_total_SpMatrix_and_SpvnpT() {
    ostringstream filename_stream_eq_SpMatrix;
    ostringstream filename_stream_eq_Spvn;
    ostringstream filename_stream_SpMatrix;
    ostringstream filename_stream_Spvn;
    ostringstream filename_stream_inte_eq_Spvn;
    ostringstream filename_stream_inte_Spvn;

    string filename = "photon_total";

    filename_stream_eq_SpMatrix << output_path << filename
                                << "_eq_SpMatrix.dat";
    filename_stream_eq_Spvn << output_path << filename << "_eq_Spvn.dat";
    filename_stream_SpMatrix << output_path << filename << "_SpMatrix.dat";
    filename_stream_Spvn << output_path << filename << "_Spvn.dat";
    filename_stream_inte_eq_Spvn << output_path << filename
                                 << "_eq_Spvn_inte.dat";
    filename_stream_inte_Spvn << output_path << filename << "_Spvn_inte.dat";

    ofstream fphoton_eq_SpMatrix(filename_stream_eq_SpMatrix.str().c_str());
    ofstream fphoton_eq_Spvn(filename_stream_eq_Spvn.str().c_str());
    ofstream fphotonSpMatrix(filename_stream_SpMatrix.str().c_str());
    ofstream fphotonSpvn(filename_stream_Spvn.str().c_str());
    ofstream fphotoninte_eq_Spvn(filename_stream_inte_eq_Spvn.str().c_str());
    ofstream fphotoninteSpvn(filename_stream_inte_Spvn.str().c_str());

    for (int m = 0; m < nm; m++) {
        for (int i=0; i < nphi; i++) {
            double phi = dilepton_QGP_LO->getPhotonphi(i);
            fphoton_eq_SpMatrix << phi << "  ";
            fphotonSpMatrix << phi << "  ";
            for (int j = 0; j < np; j++) {
                double temp_eq = 0.0;
                double temp_tot = 0.0;
                double dy = dilepton_QGP_LO->get_dy();
                for (int k = 0; k < nrapidity; k++) {
                    temp_eq += dNd2pTdphidy_eq[m][j][i][k]*dy;
                    temp_tot += dNd2pTdphidy_tot[m][j][i][k]*dy;
                }
                fphoton_eq_SpMatrix << scientific << setprecision(6) << setw(16)
                                    << temp_eq << "  ";
                fphotonSpMatrix << scientific << setprecision(6) << setw(16)
                                << temp_tot << "  ";
            }
            fphoton_eq_SpMatrix << endl;
            fphotonSpMatrix << endl;
        }
    }

    // pT differential
    for (int m = 0; m < nm; m++) {
        for (int i = 0; i < np; i++) {
            double M_ll = dilepton_QGP_LO->getDileptonMass(m);
            fphoton_eq_Spvn << scientific << setprecision(6) << setw(16)
                            << M_ll << "  ";
            fphotonSpvn << scientific << setprecision(6) << setw(16)
                        << M_ll << "  ";
            double pT = dilepton_QGP_LO->getPhotonp(i);
            fphoton_eq_Spvn << scientific << setprecision(6) << setw(16)
                            << pT << "  " << dNd2pTd2M_eq[m][i] << "  ";
            fphotonSpvn << scientific << setprecision(6) << setw(16)
                        << pT << "  " << dNd2pTd2M_tot[m][i] << "  ";
            for (int order=1; order < norder; order++) {
                fphoton_eq_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vnpT_cos_eq[order][m][i] << "  "
                                << vnpT_sin_eq[order][m][i] << "  "
                                << sqrt(pow(vnpT_cos_eq[order][m][i], 2)
                                        + pow(vnpT_sin_eq[order][m][i], 2)) << "  ";
                fphotonSpvn << scientific << setprecision(6) << setw(16)
                            << order << "   " << vnpT_cos_tot[order][m][i] << "  "
                            << vnpT_sin_tot[order][m][i] << "  "
                            << sqrt(pow(vnpT_cos_tot[order][m][i], 2)
                                    + pow(vnpT_sin_tot[order][m][i], 2)) << "  ";
            }
            fphoton_eq_Spvn << endl;
            fphotonSpvn << endl;
        }
    }

    // pT integrated
    for (int m = 0; m < nm; m++) {
        double M_ll = dilepton_QGP_LO->getDileptonMass(m);
        fphotoninte_eq_Spvn << scientific << setprecision(6) << setw(16)
                        << M_ll << "  " << dNd2Mdy_eq[m]*M_ll << "  ";
        fphotoninteSpvn << scientific << setprecision(6) << setw(16)
                    << M_ll << "  " << dNd2Mdy_tot[m]*M_ll << "  ";

        for (int order = 0; order < norder; order++) {
            fphotoninte_eq_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vn_cos_eq[order][m] << "   "
                                << vn_sin_eq[order][m] << "   "
                                << sqrt(pow(vn_cos_eq[order][m], 2)
                                        + pow(vn_sin_eq[order][m], 2)) << "  ";
            fphotoninteSpvn << scientific << setprecision(6) << setw(16)
                            << order << "   " << vn_cos_tot[order][m] << "   "
                            << vn_sin_tot[order][m] << "   "
                            << sqrt(pow(vn_cos_tot[order][m], 2)
                                    + pow(vn_sin_tot[order][m], 2)) << "  ";
        }
        fphotoninte_eq_Spvn << endl;
        fphotoninteSpvn << endl;
    }
}
