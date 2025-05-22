// Copyright 2016 Chun Shen
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
//#include <omp.h>
#include <unordered_set>
#include <vector>
#ifndef _OPENMP
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#else
#include <omp.h>
#endif

#include "Arsenal.h"
#include "Hydroinfo_h5.h"
#include "ParameterReader.h"
#include "PhotonEmission.h"
#include "QGP_LO.h"
#include "QGP_LO_analytic.h"
#include "QGP_NLO.h"
#include "HadronGas_rho_omega_phi.h"
#include "ThermalPhoton.h"
#include "tensor_trans.h"

using namespace std;
using ARSENAL::createA1DMatrix;
using ARSENAL::createA2DMatrix;
using ARSENAL::createA3DMatrix;
using ARSENAL::createA4DMatrix;
using ARSENAL::createA5DMatrix;
using ARSENAL::createA6DMatrix;
using ARSENAL::deleteA1DMatrix;
using ARSENAL::deleteA2DMatrix;
using ARSENAL::deleteA3DMatrix;
using ARSENAL::deleteA4DMatrix;
using ARSENAL::deleteA5DMatrix;
using ARSENAL::deleteA6DMatrix;
using TENSORTRANSFORM::getTransverseflow_u_mu_low;
using TENSORTRANSFORM::lorentz_boost_matrix;

using PhysConsts::eps;
using PhysConsts::hbarC;

PhotonEmission::PhotonEmission(std::shared_ptr<ParameterReader> paraRdr_in) {
  paraRdr = paraRdr_in;
  output_path = "results/";
  hydro_flag = paraRdr->getVal("hydro_flag");
  differential_flag = paraRdr->getVal("differential_flag");
  turn_off_transverse_flow = paraRdr->getVal("turn_off_transverse_flow");
  turn_on_muB_ = static_cast<int>(paraRdr->getVal("turn_on_muB", 1));
  test_code_flag = static_cast<int>(paraRdr->getVal("test_code_flag", 0));
  emission_rate_flag = paraRdr->getVal("dilepton_emission_rate");

  // omp parameters
  CORES = 1;

#ifdef _OPENMP
  CORES = omp_get_max_threads();
#endif

  set_hydroGridinfo();
  print_hydroGridinfo();

  // read the photon emission rate tables
  std::cout << " Read emission rate! " << std::endl;
  InitializePhotonEmissionRateTables();
  std::cout << " Read emission rate! end" << std::endl;

   // dN/MdMdPT2dphidy
   dNd2pTdphidy_eq_lo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_tot_lo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_pol_lambda_theta_lo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_pol_lambda_norm_lo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_pol_lambda_phi_lo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);


   
   // dN/dMdPT2dphidy
   dNd2pTd2M_eq_lo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_tot_lo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_pol_lambda_theta_lo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_pol_lambda_norm_lo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_pol_lambda_phi_lo = createA2DMatrix(nm, np, 0);


   dNd2pTd2Mdy_eq_lo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_tot_lo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_pol_lambda_theta_lo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_pol_lambda_norm_lo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_pol_lambda_phi_lo = createA3DMatrix(nm, np, nrapidity, 0);


   dNd2Mdy_eq_lo = createA1DMatrix(nm, 0); 
   dNd2Mdy_tot_lo = createA1DMatrix(nm, 0); 
   dNd2Mdy_pol_lambda_theta_lo = createA1DMatrix(nm, 0); 
   dNd2Mdy_pol_lambda_norm_lo = createA1DMatrix(nm, 0); 
   dNd2Mdy_pol_lambda_phi_lo = createA1DMatrix(nm, 0); 


   vnpT_cos_eq_lo = createA3DMatrix(norder, nm, np, 0.);
   vnpT_sin_eq_lo = createA3DMatrix(norder, nm, np, 0.);
   vnpT_cos_tot_lo = createA3DMatrix(norder, nm, np, 0.);
   vnpT_sin_tot_lo = createA3DMatrix(norder, nm, np, 0.);


   vnMpTy_cos_eq_lo = createA4DMatrix(norder, nm, np, nrapidity, 0.);
   vnMpTy_sin_eq_lo = createA4DMatrix(norder, nm, np, nrapidity, 0.);
   vnMpTy_cos_tot_lo = createA4DMatrix(norder, nm, np, nrapidity, 0.);
   vnMpTy_sin_tot_lo = createA4DMatrix(norder, nm, np, nrapidity, 0.);


   vn_cos_eq_lo = createA2DMatrix(norder, nm, 0.);
   vn_sin_eq_lo = createA2DMatrix(norder, nm, 0.);
   vn_cos_tot_lo = createA2DMatrix(norder, nm, 0.);
   vn_sin_tot_lo = createA2DMatrix(norder, nm, 0.);

  
   // NLO
   dNd2pTdphidy_eq_nlo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_tot_nlo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_pol_lambda_theta_nlo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_pol_lambda_norm_nlo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
   dNd2pTdphidy_pol_lambda_phi_nlo = createA4DMatrix(nm, np, nphi, nrapidity, 0.);


   
   // dN/dMdPT2dphidy
   dNd2pTd2M_eq_nlo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_tot_nlo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_pol_lambda_theta_nlo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_pol_lambda_norm_nlo = createA2DMatrix(nm, np, 0);
   dNd2pTd2M_pol_lambda_phi_nlo = createA2DMatrix(nm, np, 0);


   dNd2pTd2Mdy_eq_nlo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_tot_nlo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_pol_lambda_theta_nlo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_pol_lambda_norm_nlo = createA3DMatrix(nm, np, nrapidity, 0);
   dNd2pTd2Mdy_pol_lambda_phi_nlo = createA3DMatrix(nm, np, nrapidity, 0);


   dNd2Mdy_eq_nlo = createA1DMatrix(nm, 0); 
   dNd2Mdy_tot_nlo = createA1DMatrix(nm, 0); 
   dNd2Mdy_pol_lambda_theta_nlo = createA1DMatrix(nm, 0); 
   dNd2Mdy_pol_lambda_norm_nlo = createA1DMatrix(nm, 0); 
   dNd2Mdy_pol_lambda_phi_nlo = createA1DMatrix(nm, 0); 


   vnpT_cos_eq_nlo = createA3DMatrix(norder, nm, np, 0.);
   vnpT_sin_eq_nlo = createA3DMatrix(norder, nm, np, 0.);
   vnpT_cos_tot_nlo = createA3DMatrix(norder, nm, np, 0.);
   vnpT_sin_tot_nlo = createA3DMatrix(norder, nm, np, 0.);


   vnMpTy_cos_eq_nlo = createA4DMatrix(norder, nm, np, nrapidity, 0.);
   vnMpTy_sin_eq_nlo = createA4DMatrix(norder, nm, np, nrapidity, 0.);
   vnMpTy_cos_tot_nlo = createA4DMatrix(norder, nm, np, nrapidity, 0.);
   vnMpTy_sin_tot_nlo = createA4DMatrix(norder, nm, np, nrapidity, 0.);


   vn_cos_eq_nlo = createA2DMatrix(norder, nm, 0.);
   vn_sin_eq_nlo = createA2DMatrix(norder, nm, 0.);
   vn_cos_tot_nlo = createA2DMatrix(norder, nm, 0.);
   vn_sin_tot_nlo = createA2DMatrix(norder, nm, 0.);


  // if (differential_flag == 1) {

  //     std::cout<<"Start differential table" <<std::endl;

  //     dNd2pTdphidydTdtau_eq = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi,
  //     nrapidity, 0.); dNd2pTdphidydTdtau_visc = createA6DMatrix(nTcut,
  //     n_tau_cut, nm, np, nphi, nrapidity, 0.); dNd2pTdphidydTdtau_diff =
  //     createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
  //     dNd2pTdphidydTdtau_tot = createA6DMatrix(nTcut, n_tau_cut, nm, np,
  //     nphi, nrapidity, 0.);

  //     dNd2pTdphidydTdtau_eq_all = createA3DMatrix(nTcut, n_tau_cut,
  //     CORES*nrapidity*np*nphi*nm, 0.); dNd2pTdphidydTdtau_visc_all =
  //     createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
  //     dNd2pTdphidydTdtau_diff_all = createA3DMatrix(nTcut, n_tau_cut,
  //     CORES*nrapidity*np*nphi*nm, 0.); dNd2pTdphidydTdtau_tot_all =
  //     createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
  //     std::cout<<"Initialize differential table" <<std::endl;
  // }
}

PhotonEmission::~PhotonEmission() {


  deleteA4DMatrix(dNd2pTdphidy_eq_lo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_tot_lo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_pol_lambda_theta_lo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_pol_lambda_norm_lo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_pol_lambda_phi_lo, nm, np, nphi);


  deleteA2DMatrix(dNd2pTd2M_eq_lo, nm);
  deleteA2DMatrix(dNd2pTd2M_tot_lo, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_theta_lo, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_norm_lo, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_phi_lo, nm);

  deleteA3DMatrix(dNd2pTd2Mdy_eq_lo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_tot_lo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_theta_lo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_norm_lo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_phi_lo, nm, np);

  deleteA3DMatrix(vnpT_cos_eq_lo, norder, nm);
  deleteA3DMatrix(vnpT_sin_eq_lo, norder, nm);
  deleteA3DMatrix(vnpT_cos_tot_lo, norder, nm);
  deleteA3DMatrix(vnpT_sin_tot_lo, norder, nm);

  deleteA4DMatrix(vnMpTy_cos_eq_lo, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_eq_lo, norder, nm, np);


  deleteA4DMatrix(vnMpTy_cos_tot_lo, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_tot_lo, norder, nm, np);

  deleteA2DMatrix(vn_cos_eq_lo, norder);
  deleteA2DMatrix(vn_sin_eq_lo, norder);
  deleteA2DMatrix(vn_cos_tot_lo, norder);
  deleteA2DMatrix(vn_sin_tot_lo, norder);


  deleteA1DMatrix(dNd2Mdy_eq_lo); 
  deleteA1DMatrix(dNd2Mdy_tot_lo); 
  deleteA1DMatrix(dNd2Mdy_pol_lambda_theta_lo);
  deleteA1DMatrix(dNd2Mdy_pol_lambda_norm_lo); 
  deleteA1DMatrix(dNd2Mdy_pol_lambda_phi_lo); 

  //NLO
  deleteA4DMatrix(dNd2pTdphidy_eq_nlo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_tot_nlo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_pol_lambda_theta_nlo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_pol_lambda_norm_nlo, nm, np, nphi);
  deleteA4DMatrix(dNd2pTdphidy_pol_lambda_phi_nlo, nm, np, nphi);


  deleteA2DMatrix(dNd2pTd2M_eq_nlo, nm);
  deleteA2DMatrix(dNd2pTd2M_tot_nlo, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_theta_nlo, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_norm_nlo, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_phi_nlo, nm);

  deleteA3DMatrix(dNd2pTd2Mdy_eq_nlo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_tot_nlo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_theta_nlo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_norm_nlo, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_phi_nlo, nm, np);

  deleteA3DMatrix(vnpT_cos_eq_nlo, norder, nm);
  deleteA3DMatrix(vnpT_sin_eq_nlo, norder, nm);
  deleteA3DMatrix(vnpT_cos_tot_nlo, norder, nm);
  deleteA3DMatrix(vnpT_sin_tot_nlo, norder, nm);

  deleteA4DMatrix(vnMpTy_cos_eq_nlo, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_eq_nlo, norder, nm, np);


  deleteA4DMatrix(vnMpTy_cos_tot_nlo, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_tot_nlo, norder, nm, np);

  deleteA2DMatrix(vn_cos_eq_nlo, norder);
  deleteA2DMatrix(vn_sin_eq_nlo, norder);
  deleteA2DMatrix(vn_cos_tot_nlo, norder);
  deleteA2DMatrix(vn_sin_tot_nlo, norder);


  deleteA1DMatrix(dNd2Mdy_eq_nlo); 
  deleteA1DMatrix(dNd2Mdy_tot_nlo); 
  deleteA1DMatrix(dNd2Mdy_pol_lambda_theta_nlo);
  deleteA1DMatrix(dNd2Mdy_pol_lambda_norm_nlo); 
  deleteA1DMatrix(dNd2Mdy_pol_lambda_phi_nlo); 



  // if (differential_flag == 1) {
  //     nTcut = paraRdr->getVal("nTcut");
  //     n_tau_cut = paraRdr->getVal("n_tau_cut");
  //     deleteA6DMatrix(dNd2pTdphidydTdtau_eq, nTcut, n_tau_cut, nm, np, nphi);
  //     deleteA6DMatrix(dNd2pTdphidydTdtau_visc, nTcut, n_tau_cut, nm, np,
  //     nphi); deleteA6DMatrix(dNd2pTdphidydTdtau_diff, nTcut, n_tau_cut, nm,
  //     np, nphi); deleteA6DMatrix(dNd2pTdphidydTdtau_tot, nTcut, n_tau_cut,
  //     nm, np, nphi);

  //     deleteA3DMatrix(dNd2pTdphidydTdtau_eq_all, nTcut, n_tau_cut);
  //     deleteA3DMatrix(dNd2pTdphidydTdtau_visc_all, nTcut, n_tau_cut);
  //     deleteA3DMatrix(dNd2pTdphidydTdtau_diff_all, nTcut, n_tau_cut);
  //     deleteA3DMatrix(dNd2pTdphidydTdtau_tot_all, nTcut, n_tau_cut);
  //     std::cout<<"Delete differential table" <<std::endl;
  // }
}

void PhotonEmission::set_hydroGridinfo() {
  gridX0 = paraRdr->getVal("Xmin");
  gridY0 = paraRdr->getVal("Ymin");
  gridDx = paraRdr->getVal("dx");
  gridDy = paraRdr->getVal("dy");
  gridDtau = paraRdr->getVal("dTau");
  neta = paraRdr->getVal("neta");
  ETAmax = paraRdr->getVal("ETAmax");

  nm = paraRdr->getVal("nm");
  np = paraRdr->getVal("np");
  nphi = paraRdr->getVal("nphi");
  nrapidity = paraRdr->getVal("nrapidity");
  norder = paraRdr->getVal("norder");

  T_dec = paraRdr->getVal("T_dec");
  T_sw_high = paraRdr->getVal("T_sw_high");
  T_sw_low = paraRdr->getVal("T_sw_low");

  nTcut = paraRdr->getVal("nTcut");
  n_tau_cut = paraRdr->getVal("n_tau_cut");

  T_test = paraRdr->getVal("T_test");
  muB_test = paraRdr->getVal("muB_test");
  rhoB_eplusp_test = paraRdr->getVal("rhoB_eplusp_test");
  inv_eplusp_test = paraRdr->getVal("inv_eplusp_test");

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
  cout << "T_sw = " << T_sw_low << " to " << T_sw_high << " GeV." << endl;
  cout << endl;

  cout << "dilepton_emission_rate = " << emission_rate_flag << endl;

  cout << "Photon momentum: " << paraRdr->getVal("photon_q_i") << " to "
       << paraRdr->getVal("photon_q_f") << " GeV, "
       << "n_q =" << np << endl;
  cout << "Photon momentum angles: " << paraRdr->getVal("photon_phi_q_i")
       << " to " << paraRdr->getVal("photon_phi_q_f") << ", n_phi=" << nphi
       << endl;
  cout << "Photon momentum rapidity: " << paraRdr->getVal("photon_y_i")
       << " to " << paraRdr->getVal("photon_y_f") << ", n_y =" << nrapidity
       << endl;
  cout << "Dilepton invariant mass: " << paraRdr->getVal("dilepton_mass_i")
       << " to " << paraRdr->getVal("dilepton_mass_f") << " GeV, "
       << "n_m =" << nm << endl;
  cout << "Calculate individual channels in Hadron Resonance Gas phase: ";
  if (calHGIdFlag == 0) {
    cout << " No! " << endl;
  } else {
    cout << " Yes!" << endl;
  }

  if (test_code_flag == 1) {
    cout << endl;
    cout << "Test code using: " << endl;
    cout << "T_test = " << T_test << " GeV." << endl;
    cout << "muB_test = " << muB_test << " GeV." << endl;
    cout << "rhoB_eplusp_test = " << rhoB_eplusp_test << " fm^-1." << endl;
  }
}

void PhotonEmission::InitializePhotonEmissionRateTables() {

  // if(emission_rate_flag==0){
  //     dilepton_QGP_thermal = std::unique_ptr<ThermalPhoton>(
  //         new QGP_LO_analytic(paraRdr, "QGP_LO_analytic_total"));
  // }

  dilepton_QGP_thermal_LO =
      std::unique_ptr<ThermalPhoton>(new QGP_LO(paraRdr, "QGP_LO_total"));
  HadronGas_rho_meson =
      std::unique_ptr<ThermalPhoton>(new HadronGas_rho(paraRdr, "HadronGas_rho"));
  dilepton_QGP_thermal =
      std::unique_ptr<ThermalPhoton>(new QGP_NLO(paraRdr, "QGP_NLO_total"));
  dilepton_QGP_thermal->readEmissionrateFromFile(
      true); // true to read in the NLO emission table
}

double PhotonEmission::suppression_factor(double tau, double T) {
  // double lambda= 1.249; //tau_chem = 1.5
  // double lambda= 0.833; //tau_chem = 1.0
  double sig_lambda = paraRdr->getVal("sig_lambda");
  double suppress_order = paraRdr->getVal("suppress_order");
  double etaovers = paraRdr->getVal("etaovers");
  // cout<< sig_lambda<<" "<<suppress_order<<endl;
  double A = 1.91882;
  double T_fm = T / 0.19733;
  double pi_tem = 3.1415926;
  if (fabs(sig_lambda) < 0.0001) {

    double tauR = etaovers * 4. * pi_tem / T_fm;
    double factor0 = 1 - exp(-A * tau / tauR);
    return pow(factor0, suppress_order);
  } else {

    return pow(1 - exp(-A * tau / sig_lambda), suppress_order);
  }
}

void PhotonEmission::calPhotonemission_2d(void *hydroinfo_ptr_in,int
hydro_mode) {

    // hydro data read in main.cpp
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;
    hydroinfo_MUSIC_ptr =
                reinterpret_cast<Hydroinfo_MUSIC*>(hydroinfo_ptr_in);

    int hydro_flag = paraRdr->getVal("hydro_flag");

    int Hydro_2D = paraRdr->getVal("Hydro_2D");
    double Hydro_etas_i = paraRdr->getVal("Hydro_etas_i");
    double Hydro_etas_f = paraRdr->getVal("Hydro_etas_f");
    int Hydro_netas_n = paraRdr->getVal("Hydro_netas_n");
    double d_Hydro_detas = (Hydro_etas_f -
    Hydro_etas_i)/((Hydro_netas_n-1)*1.0);

    // photon momentum in the lab frame
    double M_ll[nm]; // invariant mass array
    double p_q[np], phi_q[nphi], y_q[nrapidity];
    double sin_phiq[nphi], cos_phiq[nphi];

    for (int k = 0; k < nrapidity; k++) {
        y_q[k] = dilepton_QGP_thermal->getPhotonrapidity(k);
    }
    Dy = dilepton_QGP_thermal->get_Dy(); // total rapidity range y_f - y_i

    for (int l = 0; l < np; l++) {
        p_q[l] = dilepton_QGP_thermal->getPhotonp(l);
    }
    for (int m = 0; m < nphi; m++) {
        phi_q[m] = dilepton_QGP_thermal->getPhotonphi(m);
        sin_phiq[m] = sin(phi_q[m]);
        cos_phiq[m] = cos(phi_q[m]);
    }
    for (int j = 0; j < nm; j++) {
        M_ll[j] = dilepton_QGP_thermal->getDileptonMass(j);
    }

    // get hydro grid information
    double dtau = hydroinfo_MUSIC_ptr->get_hydro_dtau();
    double dx = hydroinfo_MUSIC_ptr->get_hydro_dx();
    double deta = hydroinfo_MUSIC_ptr->get_hydro_deta();
    double eta_max = hydroinfo_MUSIC_ptr->get_hydro_eta_max();
    double Nskip_x = hydroinfo_MUSIC_ptr->get_hydro_Nskip_x();
    double Nskip_eta = hydroinfo_MUSIC_ptr->get_hydro_Nskip_eta();
    double Nskip_tau = hydroinfo_MUSIC_ptr->get_hydro_Nskip_tau();

    if (Hydro_netas_n ==1 && Hydro_2D){
         d_Hydro_detas = deta;
    }
    if(Hydro_2D){
        cout<<" calculate dilepton with 2D Hydro "<<endl;
        cout<<" Hydro_etas_i "<< Hydro_etas_i <<endl;
        cout<<" Hydro_etas_f "<< Hydro_etas_f <<endl;
        cout<<" Hydro_netas_n "<< Hydro_netas_n <<endl;
        cout<<" d_Hydro_detas "<< d_Hydro_detas <<endl;

    }

    tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
    tau_max = hydroinfo_MUSIC_ptr->get_hydro_tau_max();

    // number of fluid cells
    long int number_of_cells =(
                hydroinfo_MUSIC_ptr->get_number_of_fluid_cells_3d());
    cout << "number of cells:" << number_of_cells << endl;
    CORES = 1;
    long n = 0;
    // multi-threads setup
    long FO_chunk = number_of_cells / CORES;
    long remainder = number_of_cells  -  CORES * FO_chunk;

    cout << "Number of cores : " << CORES << endl;
    cout << "Chunk size = " << FO_chunk << endl;
    cout << "Remainder cells = " << remainder << endl;

    if(remainder != 0) FO_chunk++;

    cout << "----------------------------------------" << endl;

    // arrays to store values across cores
    // double *dNd2pTdphidy_eq_all = (double*)calloc(nrapidity*np*nphi*nm,
    // sizeof(double)); double *dNd2pTdphidy_eqT_all =
    // (double*)calloc(nrapidity*np*nphi*nm, sizeof(double)); double
    // *dNd2pTdphidy_eqL_all = (double*)calloc(nrapidity*np*nphi*nm,
    // sizeof(double)); double *dNd2pTdphidy_visc_all =
    // (double*)calloc(nrapidity*np*nphi*nm, sizeof(double)); double
    // *dNd2pTdphidy_diff_all = (double*)calloc(nrapidity*np*nphi*nm,
    // sizeof(double)); double *dNd2pTdphidy_tot_all =
    // (double*)calloc(nrapidity*np*nphi*nm, sizeof(double)); double
    // *dNd2pTdphidy_pol_lambda_theta_all =
    // (double*)calloc(nrapidity*np*nphi*nm, sizeof(double)); double
    // *dNd2pTdphidy_pol_lambda_norm_all = (double*)calloc(nrapidity*np*nphi*nm,
    // sizeof(double)); double *dNd2pTdphidy_pol_lambda_phi_all =
    // (double*)calloc(nrapidity*np*nphi*nm, sizeof(double));

    // main loop begins ...
    // loop over all fluid cells
    // subdivide bite size chunks of freezeout surface across cores
    int ncells = 0;
    double tau_now = 0.0;

    if (Hydro_2D == 0){
        Hydro_netas_n = 1;
    }

    gridX0 = paraRdr->getVal("Xmin");
    gridY0 = paraRdr->getVal("Ymin");
    gridDx = paraRdr->getVal("dx");
    gridDy = paraRdr->getVal("dy");
    gridTau0 = paraRdr->getVal("tau_start");
    gridTauf = paraRdr->getVal("tau_end");
    gridDtau = paraRdr->getVal("dTau");

    gridNx = 2*fabs(gridX0)/gridDx + 1;
    gridNy = 2*fabs(gridY0)/gridDy + 1;

    gridTauf = tau_max;

    gridTau0 = tau0;

    gridDtau = dtau;
    std::cout << "tau0 "<< gridTau0 <<std::endl;
    std::cout << "tauf "<< gridTauf <<std::endl;
    std::cout << "dtau "<< gridDtau <<std::endl;

    int nFrame = static_cast<int>((gridTauf - gridTau0)/gridDtau + 1e-15) +1;

    fluidCell *fluidCellptr = new fluidCell;

    for(int ietas = 0; ietas < Hydro_netas_n ; ietas++){

        double eta_local = Hydro_etas_i + ietas*d_Hydro_detas;

        double volume_base0 = gridDx*gridDy*gridDtau*d_Hydro_detas;
        double preprefactor = 1.0;
        if(eta_local < 0){
            continue;
        }
        if(eta_local > 0){
            preprefactor =2.0;
        }

        double vweight = 1;
        if(Hydro_netas_n > 1 && (ietas ==0 || ietas == Hydro_netas_n -1)){
            vweight = 0.5;
        }
        vweight = vweight*preprefactor;
        volume_base0 = volume_base0*vweight;
        cout << "Updated Volume base = " << volume_base0 << endl;
        cout << "eta_local = " << eta_local << endl;
        cout << "vweight = " << vweight << endl;
        cout << "d_Hydro_detas = " << d_Hydro_detas << endl;

        for (int frameId = 0; frameId < nFrame; frameId++) {
            double tau_local = gridTau0 + frameId*gridDtau;
            std::cout<<"tau: "<< tau_local <<std::endl;
            for (int xi = 0; xi < gridNx; xi++) {
                double x_local = gridX0 + xi*gridDx;
                for (int yj = 0; yj < gridNy; yj++) {
                    double y_local = gridY0 + yj*gridDy;

                    hydroinfo_MUSIC_ptr->getHydroValues(
                        x_local, y_local, 0.0, tau_local, fluidCellptr);

            // fluid velocity in Minkowski
            double flow_u_mu_low[4];
            double flow_u_mu_Min[4];

            // dilepton 4-momentum in Minkowski coordinates
            double p_lab_local[4];
            double p_lab_Min[4];

            // volume element: tau*dtau*dx*dy*deta,
            double volume = tau_local*volume_base0;

            double ed_local = fluidCellptr->ed;
            double pd_local = fluidCellptr->pressure;
            double temp_local = fluidCellptr->temperature;
            double temp_inv = 1/temp_local;

            // note that when turn_on_muB_=0, muB and rhoB are set to zero as follows 
            double muB_local = fluidCellptr->muB; 
            double rhoB_local = fluidCellptr->rhoB; 
            double inv_eplusp = 1./(ed_local+pd_local);
            double rhoB_over_eplusp = rhoB_local*inv_eplusp;

            double T_sw = 0.166 - 0.139 * pow(muB_local, 2) - 0.053 *
            pow(muB_local, 4); // Cleymans et al

            // validation setup, using some constant values to test the code
            if(test_code_flag==1){
                temp_local = T_test;
                temp_inv = 1/temp_local;
                muB_local = muB_test;
                rhoB_over_eplusp = rhoB_eplusp_test;
                inv_eplusp = inv_eplusp_test;
                volume = 1.0;
                eta_local = 0.0;
                turn_off_transverse_flow = 1;
            }

            // fluid cell is out of interest, temperature below T_dec or eta
            //larger than ETAmax 
            if (hydro_flag==2 && (temp_local < T_dec ||
            eta_local > ETAmax))
                continue;
            if (hydro_flag==22 && (temp_local < T_dec || eta_local > ETAmax))
                continue;

            if (differential_flag == 1){
                if (temp_local > T_cuthigh||temp_local <
                T_cutlow||tau_local>tau_cut_high||tau_local<tau_cut_low){
                  printf("Warning: local temperature or proper time out of "
                    "[T_cutlow, T_cuthigh] or [tau_cut_low, tau_cut_high].\n");
                  printf("The boundaries should be enlarged to encolse all hydro "
                    "cells...\n");
                  continue;
                }
            }

	    if (temp_local > T_sw) {
                // total fluid cells of QGP emission
                // #pragma omp atomic update
                ncells++;
            }

            // ueta = tau*ueta
            double ux, uy, ueta;
            if (turn_off_transverse_flow == 1) {
                ux = 0.0;
                uy = 0.0;
                ueta = 0.0;
            } else {
                ux = fluidCellptr->ux;
                uy = fluidCellptr->uy;
                ueta = fluidCellptr->ueta;
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
            std::vector<std::vector<double>> lambda_munu(4,
            std::vector<double>(4, 0.0)); lorentz_boost_matrix(lambda_munu,
            flow_u_mu_Min[0],
                flow_u_mu_Min[1], flow_u_mu_Min[2], flow_u_mu_Min[3]);

            // Wij is reduced variables Wij/(e+P), Wieta = tau*Wieta
            double pi11 = fluidCellptr->pi[1][1];
            double pi12 = fluidCellptr->pi[1][2];
            double pi13 = fluidCellptr->pi[1][3];
            double pi22 = fluidCellptr->pi[2][2];
            double pi23 = fluidCellptr->pi[2][3];
            // reconstruct all other components of the shear stress tensor
            double pi01 = (ux*pi11 + uy*pi12 + ueta*pi13)/utau;
            double pi02 = (ux*pi12 + uy*pi22 + ueta*pi23)/utau;
            double pi33 = (utau*(ux*pi01 + uy*pi02) - utau*utau*(pi11 + pi22)
                           + ueta*(ux*pi13 + uy*pi23))/(utau*utau -
                           ueta*ueta);
            double pi00 = pi11 + pi22 + pi33;
            double pi03 = (ux*pi13 + uy*pi23 + ueta*pi33)/utau;

            double bulkPi_local = fluidCellptr->bulkPi;

            // qeta = tau*qeta, qi is qi/kappa_hat
            // todo check music output
            double qx = fluidCellptr->qmu[0];
            double qy = fluidCellptr->qmu[1];
            double qeta = fluidCellptr->qmu[3];
            double qtau = (ux*qx + uy*qy + ueta*qeta)/utau;

            // prefactors in the dissipative correction
            double Cq = 0.97;
            //https://journals.aps.org/prc/pdf/10.1103/PhysRevC.89.034904
            //under eq.6 
            double prefactor_diff = 1./(4.*pow(2*M_PI, 5)) *
            temp_inv/pow(hbarC, 4); 
            double prefactor_visc = Cq*prefactor_diff;

            double spsfactor = 1.0;
            double cell_tau = tau_local;

            if(hydro_mode == 22 ){
                spsfactor = suppression_factor(cell_tau,temp_local);
            }

            // photon momentum loops
            // #pragma omp parallel for collapse(4) private(p_lab_Min,
            //p_lab_local) 
            for (int k = 0; k < nrapidity; k++) {
                for (int m = 0; m < nphi; m++) {
                    for (int l = 0; l < np; l++) {
                        for (int j = 0; j < nm; j++) {

                            int i3 = (j+(l+(m+nphi * k) * np) * nm);
                            // p_q is p_T magnitude array
                            double M_T = sqrt(p_q[l]*p_q[l]+M_ll[j]*M_ll[j]);
                            // transverse mass

                            //y_q[k]=0.0;
                            //cos_phiq[m]=1.0;sin_phiq[m]=0.0;
                            double cosh_y = cosh(y_q[k]);
                            double sinh_y = sinh(y_q[k]);

                            double cosh_y_minus_eta = cosh(y_q[k] -
                            eta_local); 
                            double sinh_y_minus_eta = sinh(y_q[k]
                            - eta_local);
                             // from Minkowski to Milne in lab frame
                            p_lab_local[0] = M_T * cosh_y_minus_eta;
                            p_lab_local[1] = p_q[l]*cos_phiq[m];
                            p_lab_local[2] = p_q[l]*sin_phiq[m];
                            p_lab_local[3] = M_T * sinh_y_minus_eta;  //tau*p^eta

                            // Minkowski four momentum in lab frame
                            p_lab_Min[0] = M_T * cosh_y;
                            p_lab_Min[1] = p_q[l]*cos_phiq[m];
                            p_lab_Min[2] = p_q[l]*sin_phiq[m];
                            p_lab_Min[3] = M_T * sinh_y;

                            // Minkowski four momentum from Lab frame to LRF frame 
                            double p_Min_lrf[4]; 
                            for (int j = 0; j < 4;j++) {
                                p_Min_lrf[j] = 0.;
                                for (int i = 0; i < 4; i++) {
                                    p_Min_lrf[j] +=
                                    lambda_munu[j][i]*p_lab_Min[i];
                                }
                            }

                          double pvec_lrf, pvec3,
                          pvec5; // Minkowski spatial magnitude |vec p| and |vec p|^3
                          pvec_lrf = sqrt(p_Min_lrf[1] * p_Min_lrf[1] +
                                        p_Min_lrf[2] * p_Min_lrf[2] +
                                        p_Min_lrf[3] * p_Min_lrf[3]);
                          pvec3 = pow(pvec_lrf, 3);
                          pvec5 = pow(pvec_lrf, 5);

                            // pi^\mu\nu p_\mu p_\nu, calculated in lab frame
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
                            // note that 1/kappa_hat included in q^\mu
                            double diff_dot_p = qtau*p_lab_local[0] -
                              qx*p_lab_local[1] - qy*p_lab_local[2] -
                              qeta*p_lab_local[3];
                            // validation setup
                            if(test_code_flag==1){
                                pi_photon = 0.0;
                                diff_dot_p = 0.0;
                            }

                            double Eq_localrest_Tb = p_Min_lrf[0];

                            double M2 = M_ll[j] * M_ll[j];
                            double visc_fac = prefactor_visc * pi_photon * M2
                            / pvec5; 
                            double diff_fac = prefactor_diff *
                            diff_dot_p * M2 / pvec3; 
                            double bulkPi_fac =
                            bulkPi_local;

                            

                            // begin to calculate thermal photon emission
                            //if (hydro_flag==2 && temp_local > T_sw) {
                                // QGP emission
                            double QGP_fraction = 1.0;
                            QGP_fraction = QGP_fraction*spsfactor;
                            muB_local = turn_on_muB_ * muB_local;
                            rhoB_over_eplusp = turn_on_muB_ * rhoB_over_eplusp;
          
                            dilepton_QGP_thermal->calThermalPhotonemission_3d(
                                    p_lab_Min, flow_u_mu_Min, Eq_localrest_Tb, M_ll[j],
                                    visc_fac, bulkPi_fac, diff_fac, temp_local, muB_local,
                                    inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,n,i3);
              
                            dilepton_QGP_thermal_LO->calThermalPhotonemission_3d(
                                    p_lab_Min, flow_u_mu_Min, Eq_localrest_Tb, M_ll[j],
                                    visc_fac, bulkPi_fac, diff_fac, temp_local, muB_local,
                                    inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,n,i3);
                            if (calHGIdFlag)
                            {
                              HadronGas_rho_meson->calThermalPhotonemission_3d(
                                      p_lab_Min, flow_u_mu_Min, Eq_localrest_Tb, M_ll[j],
                                      visc_fac, bulkPi_fac, diff_fac, temp_local, muB_local,
                                      inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,n,i3);
                            }

                            //}

                      

                        } // M_ll
                    } // p_T
                } // phi_p
            } // y
        }  //gridx
        } //gridy
    } //tau

    } //eta

    dilepton_QGP_thermal->reduce_multile_core();
    dilepton_QGP_thermal_LO->reduce_multile_core();
    HadronGas_rho_meson->reduce_multile_core();

    
    // Total number of elements that satisfy the condition
    // int total_count = ncells;
    printf("Cells above T_sw_high=%d...\n", ncells);

}

void PhotonEmission::calPhotonemission_3d(void *hydroinfo_ptr_in,
                                          int hydro_mode) {

  // hydro data read in main.cpp
  Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;
  hydroinfo_MUSIC_ptr = reinterpret_cast<Hydroinfo_MUSIC *>(hydroinfo_ptr_in);

  int hydro_flag = paraRdr->getVal("hydro_flag");

  int Hydro_2D = paraRdr->getVal("Hydro_2D");
  double Hydro_etas_i = paraRdr->getVal("Hydro_etas_i");
  double Hydro_etas_f = paraRdr->getVal("Hydro_etas_f");
  int Hydro_netas_n = paraRdr->getVal("Hydro_netas_n");
  double d_Hydro_detas =
      (Hydro_etas_f - Hydro_etas_i) / ((Hydro_netas_n - 1) * 1.0);

  // photon momentum in the lab frame
  double M_ll[nm]; // invariant mass array
  double p_q[np], phi_q[nphi], y_q[nrapidity];
  double sin_phiq[nphi], cos_phiq[nphi];

  for (int k = 0; k < nrapidity; k++) {
    y_q[k] = dilepton_QGP_thermal->getPhotonrapidity(k);
  }
  Dy = dilepton_QGP_thermal->get_Dy(); // total rapidity range y_f - y_i

  for (int l = 0; l < np; l++) {
    p_q[l] = dilepton_QGP_thermal->getPhotonp(l);
  }
  for (int m = 0; m < nphi; m++) {
    phi_q[m] = dilepton_QGP_thermal->getPhotonphi(m);
    sin_phiq[m] = sin(phi_q[m]);
    cos_phiq[m] = cos(phi_q[m]);
  }
  for (int j = 0; j < nm; j++) {
    M_ll[j] = dilepton_QGP_thermal->getDileptonMass(j);
  }

  // get hydro grid information
  double dtau = hydroinfo_MUSIC_ptr->get_hydro_dtau();
  double dx = hydroinfo_MUSIC_ptr->get_hydro_dx();
  double deta = hydroinfo_MUSIC_ptr->get_hydro_deta();
  double eta_max = hydroinfo_MUSIC_ptr->get_hydro_eta_max();
  double Nskip_x = hydroinfo_MUSIC_ptr->get_hydro_Nskip_x();
  double Nskip_eta = hydroinfo_MUSIC_ptr->get_hydro_Nskip_eta();
  double Nskip_tau = hydroinfo_MUSIC_ptr->get_hydro_Nskip_tau();

  double volume_base =
      Nskip_tau * dtau * Nskip_x * dx * Nskip_x * dx * Nskip_eta * deta;

  if (Hydro_netas_n == 1 && Hydro_2D) {
    d_Hydro_detas = deta;
  }
  if (Hydro_2D) {
    cout << " calculate dilepton with 2D Hydro " << endl;
    cout << " Hydro_etas_i " << Hydro_etas_i << endl;
    cout << " Hydro_etas_f " << Hydro_etas_f << endl;
    cout << " Hydro_netas_n " << Hydro_netas_n << endl;
    cout << " d_Hydro_detas " << d_Hydro_detas << endl;
  }

  tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
  tau_max = hydroinfo_MUSIC_ptr->get_hydro_tau_max();
  tau_cut_low = tau0 - dtau;
  tau_cut_high = tau_max + dtau;
  T_cutlow = hydroinfo_MUSIC_ptr->get_hydro_T_min();
  T_cuthigh = hydroinfo_MUSIC_ptr->get_hydro_T_max();

  // output results in (T, tau)
  double dT_cut = (T_cuthigh - T_cutlow) / (nTcut - 1);
  double dtau_cut = (tau_cut_high - tau_cut_low) / (n_tau_cut - 1);

  std::cout << "Differential table T range: " << T_cuthigh << " GeV"
            << " - " << T_cutlow << " GeV"
            << " dT:" << dT_cut << " GeV" << std::endl;
  std::cout << "Differential table tau range: " << tau_cut_high << " fm"
            << " - " << tau_cut_low << " fm"
            << " dtau:" << dtau_cut << " fm" << std::endl;

  if (differential_flag == 1) {
    if (tau_max > tau_cut_high || tau0 < tau_cut_low) {
      printf("Warning: proper time out of [tau_cut_low, tau_cut_high].\n");
      printf(
          "The boundaries should be enlarged to encolse all hydro cells...\n");
    }
  }

  // number of fluid cells
  long int number_of_cells =
      (hydroinfo_MUSIC_ptr->get_number_of_fluid_cells_3d());
  cout << "number of cells:" << number_of_cells << endl;

  // multi-threads setup
  long FO_chunk = number_of_cells / CORES;
  long remainder = number_of_cells - CORES * FO_chunk;

  cout << "Number of cores : " << CORES << endl;
  cout << "Chunk size = " << FO_chunk << endl;
  cout << "Remainder cells = " << remainder << endl;
  cout << "Volume base = " << volume_base << endl;

  if (remainder != 0)
    FO_chunk++;

  cout << "----------------------------------------" << endl;

  // arrays to store values across cores
  // double *dNd2pTdphidy_eq_all = (double*)calloc(CORES * nrapidity*np*nphi*nm,
  // sizeof(double)); double *dNd2pTdphidy_eqT_all = (double*)calloc(CORES *
  // nrapidity*np*nphi*nm, sizeof(double)); double *dNd2pTdphidy_eqL_all =
  // (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double)); double
  // *dNd2pTdphidy_visc_all = (double*)calloc(CORES * nrapidity*np*nphi*nm,
  // sizeof(double)); double *dNd2pTdphidy_diff_all = (double*)calloc(CORES *
  // nrapidity*np*nphi*nm, sizeof(double)); double *dNd2pTdphidy_tot_all =
  // (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double)); double
  // *dNd2pTdphidy_pol_lambda_theta_all = (double*)calloc(CORES *
  // nrapidity*np*nphi*nm, sizeof(double)); double
  // *dNd2pTdphidy_pol_lambda_norm_all = (double*)calloc(CORES *
  // nrapidity*np*nphi*nm, sizeof(double)); double
  // *dNd2pTdphidy_pol_lambda_phi_all = (double*)calloc(CORES *
  // nrapidity*np*nphi*nm, sizeof(double));

  // main loop begins ...
  // loop over all fluid cells
  // subdivide bite size chunks of freezeout surface across cores
  int ncells = 0;
  double tau_now = 0.0;

  if (Hydro_2D == 0) {
    Hydro_netas_n = 1;
  }

  for (int ietas = 0; ietas < Hydro_netas_n; ietas++) {

    double ietas_local = Hydro_etas_i + ietas * d_Hydro_detas;

    double volume_base0 = volume_base;
    if (Hydro_2D) {
      volume_base0 = volume_base / (Nskip_eta * deta);
      double vweight = 1;
      if (Hydro_netas_n > 1 && (ietas == 0 || ietas == Hydro_netas_n - 1)) {
        vweight = 0.5;
      }
      volume_base0 = volume_base0 * vweight * d_Hydro_detas;
      cout << "Updated Volume base = " << volume_base0 << endl;
      cout << "eta_local = " << ietas_local << endl;
      cout << "vweight = " << vweight << endl;
      cout << "d_Hydro_detas = " << d_Hydro_detas << endl;
    }

#pragma omp parallel for reduction(+ : ncells)
    for (long n = 0; n < CORES; n++) {
      long endFO = FO_chunk;

      int processedCells = 0;           // Counter for processed cells
      int printInterval = 0.25 * endFO; // Print at 25% intervals

      //for (long int icell = 0; icell <  endFO;
    //        icell++) // cell index inside each chunk
        for (long int icell = 0; icell < 10000;
            icell++) // cell index inside each chunk
       
           {
        if ((icell == endFO - 1) && (remainder != 0) && (n > remainder - 1))
          continue;

          

        long cell_id = n + icell * CORES;
        std::cout<< cell_id << std::endl;
        // fluid velocity in Minkowski
        double flow_u_mu_low[4];
        double flow_u_mu_Min[4];

        // dilepton 4-momentum in Minkowski coordinates
        double p_lab_local[4];
        double p_lab_Min[4];

        fluidCell_3D_new fluidCellptr;

        hydroinfo_MUSIC_ptr->get_hydro_cell_info_3d(cell_id, fluidCellptr);

        double tau_local = tau0 + fluidCellptr.itau * dtau;
        double eta_local = -eta_max + fluidCellptr.ieta * deta;

        if (Hydro_2D) {
          eta_local = ietas_local;
        }

#ifdef _OPENMP
        processedCells++; // Increment the counter for processed cells

        if (processedCells % printInterval == 0) {
          double progress = (processedCells * 100) / FO_chunk;
          int threadNum = omp_get_thread_num();
          std::cout << "Thread " << threadNum << ": Processing " << progress
                    << "% of cells" << std::endl;
        }
#else

        if (fabs(tau_now - tau_local) > 1e-10) {
          tau_now = tau_local;
          cout << "Calculating tau = " << setw(4) << setprecision(3) << tau_now
               << " fm/c..." << endl;
        }
#endif

        // volume element: tau*dtau*dx*dy*deta,
        double volume = tau_local * volume_base0;

        double ed_local = fluidCellptr.ed;
        double pd_local = fluidCellptr.pressure;
        double temp_local = fluidCellptr.temperature;
        double temp_inv = 1 / temp_local;

        // note that when turn_on_muB_=0, muB and rhoB are set to zero as
        // follows
        double muB_local = fluidCellptr.muB;
        double rhoB_local = fluidCellptr.rhoB;
        double inv_eplusp = 1. / (ed_local + pd_local);
        double rhoB_over_eplusp = rhoB_local * inv_eplusp;

        double T_sw = 0.166 - 0.139 * pow(muB_local, 2) -
                      0.053 * pow(muB_local, 4); // Cleymans et al

        // validation setup, using some constant values to test the code
        if (test_code_flag == 1) {
          temp_local = T_test;
          temp_inv = 1 / temp_local;
          muB_local = muB_test;
          rhoB_over_eplusp = rhoB_eplusp_test;
          inv_eplusp = inv_eplusp_test;
          volume = 1.0;
          eta_local = 0.0;
          turn_off_transverse_flow = 1;
        }

        // fluid cell is out of interest, temperature below T_dec or eta larger
        // than ETAmax
        if (hydro_flag == 2 && (temp_local < T_dec || eta_local > ETAmax))
          continue;
        if (hydro_flag == 22 && (temp_local < T_dec || eta_local > ETAmax))
          continue;

        if (differential_flag == 1) {
          if (temp_local > T_cuthigh || temp_local < T_cutlow ||
              tau_local > tau_cut_high || tau_local < tau_cut_low) {
            printf("Warning: local temperature or proper time out of "
                   "[T_cutlow, T_cuthigh] or [tau_cut_low, tau_cut_high].\n");
            printf("The boundaries should be enlarged to encolse all hydro "
                   "cells...\n");
            continue;
          }
        }

        if (temp_local > T_sw) {
          // total fluid cells of QGP emission
          // #pragma omp atomic update
          ncells++;
        }

        // indices for the output in (T, tau)
        int idx_T = (int)std::floor((temp_local - T_cutlow) / dT_cut + eps);
        int idx_tau =
            (int)std::floor((tau_local - tau_cut_low) / dtau_cut + eps);

        // ueta = tau*ueta
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

        double utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);

        flow_u_mu_low[0] = utau;
        flow_u_mu_low[1] = -ux;
        flow_u_mu_low[2] = -uy;
        flow_u_mu_low[3] = -ueta;

        double cosh_eta = cosh(eta_local);
        double sinh_eta = sinh(eta_local);

        // 4 velocity in Minkowski in lab frame
        flow_u_mu_Min[0] = cosh_eta * utau + sinh_eta * ueta;
        flow_u_mu_Min[1] = ux;
        flow_u_mu_Min[2] = uy;
        flow_u_mu_Min[3] = sinh_eta * utau + cosh_eta * ueta;

        // Lorentz boost matrix from lab to local rest frame
        std::vector<std::vector<double>> lambda_munu(
            4, std::vector<double>(4, 0.0));
        lorentz_boost_matrix(lambda_munu, flow_u_mu_Min[0], flow_u_mu_Min[1],
                             flow_u_mu_Min[2], flow_u_mu_Min[3]);

        // Wij is reduced variables Wij/(e+P), Wieta = tau*Wieta
        double pi11 = fluidCellptr.pi11;
        double pi12 = fluidCellptr.pi12;
        double pi13 = fluidCellptr.pi13;
        double pi22 = fluidCellptr.pi22;
        double pi23 = fluidCellptr.pi23;
        // reconstruct all other components of the shear stress tensor
        double pi01 = (ux * pi11 + uy * pi12 + ueta * pi13) / utau;
        double pi02 = (ux * pi12 + uy * pi22 + ueta * pi23) / utau;
        double pi33 =
            (utau * (ux * pi01 + uy * pi02) - utau * utau * (pi11 + pi22) +
             ueta * (ux * pi13 + uy * pi23)) /
            (utau * utau - ueta * ueta);
        double pi00 = pi11 + pi22 + pi33;
        double pi03 = (ux * pi13 + uy * pi23 + ueta * pi33) / utau;

        double bulkPi_local = fluidCellptr.bulkPi;

        // qeta = tau*qeta, qi is qi/kappa_hat
        double qx = fluidCellptr.qx;
        double qy = fluidCellptr.qy;
        double qeta = fluidCellptr.qz;
        double qtau = (ux * qx + uy * qy + ueta * qeta) / utau;

        // prefactors in the dissipative correction
        double Cq = 0.97; // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.89.034904
                          // under eq.6
        double prefactor_diff =
            1. / (4. * pow(2 * M_PI, 5)) * temp_inv / pow(hbarC, 4);
        double prefactor_visc = Cq * prefactor_diff;

        double spsfactor = 1.0;
        double cell_tau = tau_local;

        if (hydro_mode == 22) {
          spsfactor = suppression_factor(cell_tau, temp_local);
        }

        // photon momentum loops
        // #pragma omp parallel for collapse(4) private(p_lab_Min, p_lab_local)
        for (int k = 0; k < nrapidity; k++) {
          for (int m = 0; m < nphi; m++) {
            for (int l = 0; l < np; l++) {
              for (int j = 0; j < nm; j++) {

                int i3 = (j + (l + (m + nphi * k) * np) * nm);
                // p_q is p_T magnitude array
                double M_T = sqrt(p_q[l] * p_q[l] +
                                  M_ll[j] * M_ll[j]); // transverse mass

                // y_q[k]=0.0;
                // cos_phiq[m]=1.0;sin_phiq[m]=0.0;
                double cosh_y = cosh(y_q[k]);
                double sinh_y = sinh(y_q[k]);

                double cosh_y_minus_eta = cosh(y_q[k] - eta_local);
                double sinh_y_minus_eta = sinh(y_q[k] - eta_local);
                // from Minkowski to Milne in lab frame
                p_lab_local[0] = M_T * cosh_y_minus_eta;
                p_lab_local[1] = p_q[l] * cos_phiq[m];
                p_lab_local[2] = p_q[l] * sin_phiq[m];
                p_lab_local[3] = M_T * sinh_y_minus_eta; // tau*p^eta

                // Minkowski four momentum in lab frame
                p_lab_Min[0] = M_T * cosh_y;
                p_lab_Min[1] = p_q[l] * cos_phiq[m];
                p_lab_Min[2] = p_q[l] * sin_phiq[m];
                p_lab_Min[3] = M_T * sinh_y;

                // Minkowski four momentum from Lab frame to LRF frame
                double p_Min_lrf[4];
                for (int j = 0; j < 4; j++) {
                  p_Min_lrf[j] = 0.;
                  // std::cout << j <<" wxy 2 "<<std::endl;

                  for (int i = 0; i < 4; i++) {
                    p_Min_lrf[j] += lambda_munu[j][i] * p_lab_Min[i];
                  }
                }
                // std::cout << j <<" wxy 1 "<<std::endl;

                double pvec_lrf, pvec3,
                    pvec5; // Minkowski spatial magnitude |vec p| and |vec p|^3
                pvec_lrf = sqrt(p_Min_lrf[1] * p_Min_lrf[1] +
                                p_Min_lrf[2] * p_Min_lrf[2] +
                                p_Min_lrf[3] * p_Min_lrf[3]);
                pvec3 = pow(pvec_lrf, 3);
                pvec5 = pow(pvec_lrf, 5);

                // pi^\mu\nu p_\mu p_\nu, calculated in lab frame
                double pi_photon =
                    (p_lab_local[0] * p_lab_local[0] * pi00 -
                     2. * p_lab_local[0] * p_lab_local[1] * pi01 -
                     2. * p_lab_local[0] * p_lab_local[2] * pi02 -
                     2. * p_lab_local[0] * p_lab_local[3] * pi03 +
                     p_lab_local[1] * p_lab_local[1] * pi11 +
                     2. * p_lab_local[1] * p_lab_local[2] * pi12 +
                     2. * p_lab_local[1] * p_lab_local[3] * pi13 +
                     p_lab_local[2] * p_lab_local[2] * pi22 +
                     2. * p_lab_local[2] * p_lab_local[3] * pi23 +
                     p_lab_local[3] * p_lab_local[3] * pi33);

                // dot product of diffusion and dilepton momentum, calculated in
                // lab frame note that 1/kappa_hat included in q^\mu
                double diff_dot_p = qtau * p_lab_local[0] -
                                    qx * p_lab_local[1] - qy * p_lab_local[2] -
                                    qeta * p_lab_local[3];
                // validation setup
                if (test_code_flag == 1) {
                  pi_photon = 0.0;
                  diff_dot_p = 0.0;
                }

                double Eq_localrest_Tb = p_Min_lrf[0];

                double M2 = M_ll[j] * M_ll[j];
                double visc_fac = prefactor_visc * pi_photon * M2 / pvec5;
                double diff_fac = prefactor_diff * diff_dot_p * M2 / pvec3;
                double bulkPi_fac = bulkPi_local;

                // begin to calculate thermal photon emission
                //if (hydro_flag == 2 && temp_local > T_sw) {
                  // QGP emission
                  double QGP_fraction = 1.0;
                  QGP_fraction = QGP_fraction*spsfactor;
                  muB_local = turn_on_muB_ * muB_local;
                  rhoB_over_eplusp = turn_on_muB_ * rhoB_over_eplusp;

                  dilepton_QGP_thermal->calThermalPhotonemission_3d(
                      p_lab_Min, flow_u_mu_Min, Eq_localrest_Tb, M_ll[j],
                      visc_fac, bulkPi_fac, diff_fac, temp_local, muB_local,
                      inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,n,i3);

                  dilepton_QGP_thermal_LO->calThermalPhotonemission_3d(
                      p_lab_Min, flow_u_mu_Min, Eq_localrest_Tb, M_ll[j],
                      visc_fac, bulkPi_fac, diff_fac, temp_local, muB_local,
                      inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,n,i3);
                //}
                //else{
                  
                  //double QGP_fraction = 1.0;
                  if (calHGIdFlag)
                  {
                  HadronGas_rho_meson->calThermalPhotonemission_3d(
                    p_lab_Min, flow_u_mu_Min, Eq_localrest_Tb, M_ll[j],
                    visc_fac, bulkPi_fac, diff_fac, temp_local, muB_local,
                    inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,n,i3);
                  }

                //}

                // add contributions from QGP and Hadronic matter, etc
                // these are dN/(MdM pTdpTdphi dy)


                // if (differential_flag == 1) {
                //     dNd2pTdphidydTdtau_eq_all[idx_T][idx_tau][n+CORES*i3] +=
                //     dNd2pTdphidy_cell_eq*spsfactor;
                //     dNd2pTdphidydTdtau_visc_all[idx_T][idx_tau][n+CORES*i3]
                //     += dNd2pTdphidy_cell_visc*spsfactor;
                //     dNd2pTdphidydTdtau_diff_all[idx_T][idx_tau][n+CORES*i3]
                //     += dNd2pTdphidy_cell_diff*spsfactor;
                //     dNd2pTdphidydTdtau_tot_all[idx_T][idx_tau][n+CORES*i3] +=
                //     dNd2pTdphidy_cell_tot*spsfactor;
                // }

              } // M_ll
            }   // p_T
          }     // phi_p
        }       // y
      }         // fluid cell
    }           // cores
  }

  dilepton_QGP_thermal->reduce_multile_core();
  dilepton_QGP_thermal_LO->reduce_multile_core();
  HadronGas_rho_meson->reduce_multile_core();
  
  // Total number of elements that satisfy the condition
  // int total_count = ncells;
  printf("Cells above T_sw_high=%d...\n", ncells);

}

void PhotonEmission::calPhoton_SpvnpT_individualchannel() {
  dilepton_QGP_thermal->calPhoton_SpvnpT_shell();
  dilepton_QGP_thermal_LO->calPhoton_SpvnpT_shell();
  HadronGas_rho_meson->calPhoton_SpvnpT_shell();

  if (differential_flag == 1) {
    dilepton_QGP_thermal->calPhoton_SpMatrix_dTdtau(
        dNd2pTdphidydTdtau_eq, dNd2pTdphidydTdtau_visc, dNd2pTdphidydTdtau_diff,
        dNd2pTdphidydTdtau_tot);
    dilepton_QGP_thermal->calPhoton_Spectra_dTdtau();
    dilepton_QGP_thermal->calPhoton_Spvn_dTdtau();
  }
}

void PhotonEmission::outputPhotonSpvn_individualchannel(std::string type_str ) {
  dilepton_QGP_thermal->outputPhoton_SpvnpT_shell(output_path,type_str);
  dilepton_QGP_thermal_LO->outputPhoton_SpvnpT_shell(output_path,type_str);
  HadronGas_rho_meson->outputPhoton_SpvnpT_shell(output_path,type_str);

  if (differential_flag == 1) {
    dilepton_QGP_thermal->outputPhoton_Spectra_dTdtau(
        output_path, T_cuthigh, T_cutlow, tau_cut_high, tau_cut_low);

    dilepton_QGP_thermal->outputPhoton_Spvn_dTdtau(
        output_path, T_cuthigh, T_cutlow, tau_cut_high, tau_cut_low);

    dilepton_QGP_thermal->outputPhoton_Spectra_full_diff(
        output_path, T_cuthigh, T_cutlow, tau_cut_high, tau_cut_low);
  }
}

void PhotonEmission::calPhoton_total_Spvn() {

    for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
        for (int j = 0; j < nphi; j++) {
          for (int k = 0; k < nrapidity; k++) {
            dNd2pTdphidy_eq_lo[m][i][j][k] = dilepton_QGP_thermal_LO->get_dNd2pTdphidy_eq(m,i,j,k)
                                           + HadronGas_rho_meson->get_dNd2pTdphidy_eq(m,i,j,k) ;
            
            dNd2pTdphidy_tot_lo[m][i][j][k] = dilepton_QGP_thermal_LO->get_dNd2pTdphidy_tot(m,i,j,k)
                                            + HadronGas_rho_meson->get_dNd2pTdphidy_tot(m,i,j,k) ;
            dNd2pTdphidy_pol_lambda_theta_lo[m][i][j][k] = dilepton_QGP_thermal_LO->get_dNd2pTdphidy_pol_lambda_theta(m,i,j,k) + HadronGas_rho_meson->get_dNd2pTdphidy_pol_lambda_theta(m,i,j,k) ;

            dNd2pTdphidy_pol_lambda_norm_lo[m][i][j][k] = dilepton_QGP_thermal_LO->get_dNd2pTdphidy_pol_lambda_norm(m,i,j,k) + HadronGas_rho_meson->get_dNd2pTdphidy_pol_lambda_norm(m,i,j,k) ;

            dNd2pTdphidy_pol_lambda_phi_lo[m][i][j][k] = dilepton_QGP_thermal_LO->get_dNd2pTdphidy_pol_lambda_phi(m,i,j,k) + HadronGas_rho_meson->get_dNd2pTdphidy_pol_lambda_phi(m,i,j,k) ;


            dNd2pTdphidy_eq_nlo[m][i][j][k] = dilepton_QGP_thermal->get_dNd2pTdphidy_eq(m,i,j,k)
            + HadronGas_rho_meson->get_dNd2pTdphidy_eq(m,i,j,k) ;

            dNd2pTdphidy_tot_nlo[m][i][j][k] = dilepton_QGP_thermal->get_dNd2pTdphidy_tot(m,i,j,k)
             + HadronGas_rho_meson->get_dNd2pTdphidy_tot(m,i,j,k) ;
            dNd2pTdphidy_pol_lambda_theta_nlo[m][i][j][k] = dilepton_QGP_thermal->get_dNd2pTdphidy_pol_lambda_theta(m,i,j,k) + HadronGas_rho_meson->get_dNd2pTdphidy_pol_lambda_theta(m,i,j,k) ;

            dNd2pTdphidy_pol_lambda_norm_nlo[m][i][j][k] = dilepton_QGP_thermal->get_dNd2pTdphidy_pol_lambda_norm(m,i,j,k) + HadronGas_rho_meson->get_dNd2pTdphidy_pol_lambda_norm(m,i,j,k) ;

            dNd2pTdphidy_pol_lambda_phi_nlo[m][i][j][k] = dilepton_QGP_thermal->get_dNd2pTdphidy_pol_lambda_phi(m,i,j,k) + HadronGas_rho_meson->get_dNd2pTdphidy_pol_lambda_phi(m,i,j,k) ;


        }
      }
    }
  }

  calPhoton_total_Spvn_shell();

}


void PhotonEmission::calPhoton_total_Spvn_shell() {

    calPhoton_SpvnpT(dNd2pTdphidy_eq_lo, vnpT_cos_eq_lo,vnpT_sin_eq_lo,
                    vn_cos_eq_lo, vn_sin_eq_lo,
                    dNd2Mdy_eq_lo, dNd2pTd2M_eq_lo,
                    dNd2pTd2Mdy_eq_lo,vnMpTy_cos_eq_lo,vnMpTy_sin_eq_lo);
    calPhoton_SpvnpT(dNd2pTdphidy_tot_lo, vnpT_cos_tot_lo,vnpT_sin_tot_lo,
                      vn_cos_tot_lo, vn_sin_tot_lo,
                      dNd2Mdy_tot_lo, dNd2pTd2M_tot_lo,
                      dNd2pTd2Mdy_tot_lo,vnMpTy_cos_tot_lo,vnMpTy_sin_tot_lo);

    calPhoton_SpvnpT(dNd2pTdphidy_eq_nlo, vnpT_cos_eq_nlo,vnpT_sin_eq_nlo,
                        vn_cos_eq_nlo, vn_sin_eq_nlo,
                        dNd2Mdy_eq_nlo, dNd2pTd2M_eq_nlo,
                        dNd2pTd2Mdy_eq_nlo,vnMpTy_cos_eq_nlo,vnMpTy_sin_eq_nlo);
    calPhoton_SpvnpT(dNd2pTdphidy_tot_nlo, vnpT_cos_tot_nlo,vnpT_sin_tot_nlo,
                        vn_cos_tot_nlo, vn_sin_tot_nlo,
                        dNd2Mdy_tot_nlo, dNd2pTd2M_tot_nlo,
                        dNd2pTd2Mdy_tot_nlo,vnMpTy_cos_tot_nlo,vnMpTy_sin_tot_nlo);

    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_theta_lo,dNd2pTd2M_pol_lambda_theta_lo,dNd2Mdy_pol_lambda_theta_lo,dNd2pTd2Mdy_pol_lambda_theta_lo);
    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_norm_lo,dNd2pTd2M_pol_lambda_norm_lo,dNd2Mdy_pol_lambda_norm_lo,dNd2pTd2Mdy_pol_lambda_norm_lo);
    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_phi_lo,dNd2pTd2M_pol_lambda_phi_lo,dNd2Mdy_pol_lambda_phi_lo,dNd2pTd2Mdy_pol_lambda_phi_lo);


    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_theta_nlo,dNd2pTd2M_pol_lambda_theta_nlo,dNd2Mdy_pol_lambda_theta_nlo,dNd2pTd2Mdy_pol_lambda_theta_nlo);
    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_norm_nlo,dNd2pTd2M_pol_lambda_norm_nlo,dNd2Mdy_pol_lambda_norm_nlo,dNd2pTd2Mdy_pol_lambda_norm_nlo);
    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_phi_nlo,dNd2pTd2M_pol_lambda_phi_nlo,dNd2Mdy_pol_lambda_phi_nlo,dNd2pTd2Mdy_pol_lambda_phi_nlo);


}

void PhotonEmission::calPhoton_SpvnpT_pol(double ****dNd2pTdphidy_pol_lambda_theta, double **dNd2pTd2M_pol_lambda_theta, double*dNd2Mdy_pol_lambda_theta, double*** dNd2pTd2Mdy_pol_lambda_theta ){

  double dy = dilepton_QGP_thermal->get_dy();
  for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
          double p_local = dilepton_QGP_thermal->getPhotonp(i); // p_T
          double pweight_local = dilepton_QGP_thermal->getPhoton_pweight(i);
          for (int j = 0; j < nphi; j++) {
              double phi_local = dilepton_QGP_thermal->getPhotonphi(j); // phi
              double phi_weight_local = dilepton_QGP_thermal->getPhoton_phiweight(j);
              
              for (int k = 0; k < nrapidity; k++) {
                double y_weight_local = dilepton_QGP_thermal->getPhoton_yweight(k); // y
                // below dy/Dy to make everything rapidity density, Dy is the rapidity
                // range
                double weight = phi_weight_local * y_weight_local * dy / Dy;
                // integrate over rapidity, azimuthal angle of momentum
                dNd2pTd2M_pol_lambda_theta[m][i] += dNd2pTdphidy_pol_lambda_theta[m][i][j][k]*weight;
              }
          }
          // pT integrated spectra, dN/(MdMdy)
          dNd2Mdy_pol_lambda_theta[m] += dNd2pTd2M_pol_lambda_theta[m][i]* p_local * pweight_local;
      }
  }



  for (int k = 0; k < nrapidity; k++){
      double y_weight_local =dilepton_QGP_thermal-> getPhoton_yweight(k); // y
      for (int m = 0; m < nm; m++) {
          for (int i = 0; i < np; i++) {
              double p_local = dilepton_QGP_thermal->getPhotonp(i); // p_T
              double pweight_local = dilepton_QGP_thermal->getPhoton_pweight(i);
              
              for (int j = 0; j < nphi; j++) {
                  double phi_local = dilepton_QGP_thermal->getPhotonphi(j); // phi
                  double phi_weight_local = dilepton_QGP_thermal->getPhoton_phiweight(j);
                  dNd2pTd2Mdy_pol_lambda_theta[m][i][k] += dNd2pTdphidy_pol_lambda_theta[m][i][j][k]*phi_weight_local;
              
              }
          }
      }
  }

  
}


void PhotonEmission::calPhoton_SpvnpT( double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
                                       double **vn_cos_eq, double **vn_sin_eq,
                                       double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
                                      double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq) {
    double dy = dilepton_QGP_thermal->get_dy();
    for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
        double p_local = dilepton_QGP_thermal->getPhotonp(i); // p_T
        double pweight_local = dilepton_QGP_thermal->getPhoton_pweight(i);
        for (int j = 0; j < nphi; j++) {
          double phi_local = dilepton_QGP_thermal->getPhotonphi(j); // phi
          double phi_weight_local = dilepton_QGP_thermal->getPhoton_phiweight(j);
          for (int k = 0; k < nrapidity; k++) {
            double y_weight_local = dilepton_QGP_thermal->getPhoton_yweight(k); // y
            // below dy/Dy to make everything rapidity density, Dy is the rapidity
            // range
            double weight = phi_weight_local * y_weight_local * dy / Dy;
            // integrate over rapidity, azimuthal angle of momentum
            dNd2pTd2M_eq[m][i] += dNd2pTdphidy_eq[m][i][j][k] * weight;
          
  
            for (int order = 0; order < norder; order++) {
              vnpT_cos_eq[order][m][i] +=
                  (dNd2pTdphidy_eq[m][i][j][k] * weight * cos(order * phi_local));
              vnpT_sin_eq[order][m][i] +=
                  (dNd2pTdphidy_eq[m][i][j][k] * weight * sin(order * phi_local));
            
            }
          }
        }
  
        for (int order = 0; order < norder; order++) {
          // m integrated flow, numerator
          vn_cos_eq[order][m] += vnpT_cos_eq[order][m][i] * p_local * pweight_local;
          vn_sin_eq[order][m] += vnpT_sin_eq[order][m][i] * p_local * pweight_local;
        
  
          // // pT differential flow
          vnpT_cos_eq[order][m][i] =
              vnpT_cos_eq[order][m][i] / dNd2pTd2M_eq[m][i];    
          vnpT_sin_eq[order][m][i] =
              vnpT_sin_eq[order][m][i] / dNd2pTd2M_eq[m][i];
        }
  
        // pT integrated spectra, dN/(MdMdy)
        dNd2Mdy_eq[m] += dNd2pTd2M_eq[m][i] * p_local * pweight_local;
        
  
        // pT differential spectra, dN/(2pi pTdpT MdM dy)
        dNd2pTd2M_eq[m][i] = dNd2pTd2M_eq[m][i] / (2 * M_PI);
      
      }
    }
  
    // pT integrated flow
    for (int m = 0; m < nm; m++) {
      for (int order = 1; order < norder; order++) {
        vn_cos_eq[order][m] = vn_cos_eq[order][m] / dNd2Mdy_eq[m];
        vn_sin_eq[order][m] = vn_sin_eq[order][m] / dNd2Mdy_eq[m];
      }
    }
  
    for (int k = 0; k < nrapidity; k++) {
      double y_weight_local = dilepton_QGP_thermal->getPhoton_yweight(k); // y
      for (int m = 0; m < nm; m++) {
        for (int i = 0; i < np; i++) {
          double p_local = dilepton_QGP_thermal->getPhotonp(i); // p_T
          double pweight_local = dilepton_QGP_thermal->getPhoton_pweight(i);
  
          for (int j = 0; j < nphi; j++) {
            double phi_local = dilepton_QGP_thermal->getPhotonphi(j); // phi
            double phi_weight_local = dilepton_QGP_thermal->getPhoton_phiweight(j);
  
            dNd2pTd2Mdy_eq[m][i][k] += dNd2pTdphidy_eq[m][i][j][k] * phi_weight_local;
      
  
            for (int order = 0; order < norder; order++) {
              vnMpTy_cos_eq[order][m][i][k] +=
                  (dNd2pTdphidy_eq[m][i][j][k] * phi_weight_local * cos(order * phi_local));
            
              vnMpTy_sin_eq[order][m][i][k] +=
                  (dNd2pTdphidy_eq[m][i][j][k] * phi_weight_local * sin(order * phi_local));
              
            }
          }
  
          for (int order = 0; order < norder; order++) {
  
            vnMpTy_cos_eq[order][m][i][k] =
                vnMpTy_cos_eq[order][m][i][k] / dNd2pTd2Mdy_eq[m][i][k];
            vnMpTy_sin_eq[order][m][i][k] =
                vnMpTy_sin_eq[order][m][i][k] / dNd2pTd2Mdy_eq[m][i][k];
          }
        }
      }
    }


}

void PhotonEmission::outputPhoton_SpvnpT(string path, string type_str,
  double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
  double **vn_cos_eq, double **vn_sin_eq,
  double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
  double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq) {

    ostringstream filename_stream_eq_SpMatrix;
    ostringstream filename_stream_eq_Spvn;
    ostringstream filename_stream_eq_inte_Spvn;
    ostringstream filename_stream_eq_SpMatrix_dy;
    string file_end = ".dat";
    filename_stream_eq_SpMatrix << path <<  type_str <<"_SpMatrix"
    << file_end;
    filename_stream_eq_Spvn << path  << type_str << "_Spvn" << file_end;
    filename_stream_eq_inte_Spvn << path  << type_str << "_Spvn_inte"
    << file_end;
    filename_stream_eq_SpMatrix_dy << path  << type_str << "_Spvn_MpTy"
     << file_end;

    ofstream fphoton_eq_SpMatrix_dy(filename_stream_eq_SpMatrix_dy.str().c_str());
    ofstream fphoton_eq_SpMatrix(filename_stream_eq_SpMatrix.str().c_str());
    ofstream fphoton_eq_Spvn(filename_stream_eq_Spvn.str().c_str());
    ofstream fphoton_eq_inte_Spvn(filename_stream_eq_inte_Spvn.str().c_str());

    double dy = dilepton_QGP_thermal->get_dy();
    for (int m = 0; m < nm; m++) {
      for (int i = 0; i < nphi; i++) {
        double phi_local = dilepton_QGP_thermal->getPhotonphi(i);
        fphoton_eq_SpMatrix << phi_local << "  ";
        for (int j = 0; j < np; j++) {
          double temp_eq = 0.0; 
          for (int k = 0; k < nrapidity; k++) {
            double y_weight_local = dilepton_QGP_thermal->getPhoton_yweight(k);
            
            double weight = y_weight_local * dy / Dy;
            temp_eq += dNd2pTdphidy_eq[m][j][i][k] * weight;
          }
          fphoton_eq_SpMatrix << scientific << setprecision(6) << setw(16)
                              << temp_eq << "  ";
        }
        fphoton_eq_SpMatrix << endl;
      }
    }


    for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
        for (int k = 0; k < nrapidity; k++) {
          double mll_local = dilepton_QGP_thermal->getDileptonMass(m); // Mll
          double p_local = dilepton_QGP_thermal->getPhotonp(i);        // p_T
          double y_local = dilepton_QGP_thermal->getPhotonrapidity(k);
  
          fphoton_eq_SpMatrix_dy << scientific << setprecision(6) << setw(16)
                                 << mll_local << "  " << p_local << " " << y_local
                                 << " " << dNd2pTd2Mdy_eq[m][i][k] << "  ";
          for (int order = 1; order < norder; order++) {
  
            fphoton_eq_SpMatrix_dy << scientific << setprecision(6) << setw(16)
                                   << order << "  "
                                   << vnMpTy_cos_eq[order][m][i][k] << " "
                                   << vnMpTy_sin_eq[order][m][i][k] << " "
                                   << sqrt(pow(vnMpTy_cos_eq[order][m][i][k], 2) +
                                           pow(vnMpTy_sin_eq[order][m][i][k], 2))
                                   << "  ";
          }
  
          fphoton_eq_SpMatrix_dy << endl;
        }
      }
    }

      // pT differential, dN/(2pi pTdpT MdM dy) and vn(M, pT)
for (int m = 0; m < nm; m++) {
  for (int i = 0; i < np; i++) {
    double M_ll_local = dilepton_QGP_thermal->getDileptonMass(m);
    double pT_local = dilepton_QGP_thermal->getPhotonp(i);
    fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) << M_ll_local
                    << "  " << pT_local << "  " << dNd2pTd2M_eq[m][i] << "  ";


    for (int order = 1; order < norder; order++) {
      fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) << order
                      << "   " << vnpT_cos_eq[order][m][i] << "  "
                      << vnpT_sin_eq[order][m][i] << "  "
                      << sqrt(pow(vnpT_cos_eq[order][m][i], 2) +
                              pow(vnpT_sin_eq[order][m][i], 2))
                      << "  ";
    
    }
    fphoton_eq_Spvn << endl;
  }
}

  // pT integrated
for (int m = 0; m < nm; m++) {
double M_ll_local = dilepton_QGP_thermal->getDileptonMass(m);
// get dN/dMdy from dN/MdMdy
fphoton_eq_inte_Spvn << scientific << setprecision(6) << setw(16) << M_ll_local
                   << "  " << dNd2Mdy_eq[m] << "  ";

for (int order = 0; order < norder; order++) {
fphoton_eq_inte_Spvn << scientific << setprecision(6) << setw(16) << order
                     << "   " << vn_cos_eq[order][m] << "   "
                     << vn_sin_eq[order][m] << "   "
                     << sqrt(pow(vn_cos_eq[order][m], 2) +
                             pow(vn_sin_eq[order][m], 2))
                     << "  ";
}
fphoton_eq_inte_Spvn << endl;
}

}



void PhotonEmission::outputPhoton_SpvnpT_pol(std::string path, std::string type_str, double **** dNd2pTdphidy_pol_lambda_theta,double *** dNd2pTd2Mdy_pol_lambda_theta, double** dNd2pTd2M_pol_lambda_theta,double* dNd2Mdy_pol_lambda_theta, double* dNd2Mdy_pol_lambda_norm) {

  ostringstream filename_stream_pol_lambda_theta_SpMatrix;
  ostringstream filename_stream_pol_lambda_theta_Spvn;
  ostringstream filename_stream_pol_lambda_theta_inte_Spvn;
  

  ostringstream filename_stream_pol_lambda_theta_SpMatrix_dy;
  

  string file_end = ".dat";
  filename_stream_pol_lambda_theta_SpMatrix << path <<  type_str <<"_SpMatrix" << file_end;
  filename_stream_pol_lambda_theta_Spvn << path << type_str << "_Spvn" << file_end;
  filename_stream_pol_lambda_theta_inte_Spvn << path << type_str << "_Spvn_inte"<< file_end;
  filename_stream_pol_lambda_theta_SpMatrix_dy << path << type_str << "_Spvn_MpTy"<< file_end;



  ofstream fphoton_pol_lambda_theta_SpMatrix_dy(filename_stream_pol_lambda_theta_SpMatrix_dy.str().c_str());
  ofstream fphoton_pol_lambda_theta_SpMatrix(filename_stream_pol_lambda_theta_SpMatrix.str().c_str());
  ofstream fphoton_pol_lambda_theta_Spvn(filename_stream_pol_lambda_theta_Spvn.str().c_str());
  ofstream fphoton_pol_lambda_theta_inte_Spvn(filename_stream_pol_lambda_theta_inte_Spvn.str().c_str());

  double dy = dilepton_QGP_thermal->get_dy();
  for (int m = 0; m < nm; m++) {
      for (int i=0; i < nphi; i++) {
          double phi_local = dilepton_QGP_thermal->getPhotonphi(i);
          fphoton_pol_lambda_theta_SpMatrix << phi_local << "  ";
         
          for (int j = 0; j < np; j++) {
              double temp_pol_lambda_theta = 0.0;

              for (int k = 0; k < nrapidity; k++) {
                  double y_weight_local = dilepton_QGP_thermal->getPhoton_yweight(k);
                  // below dy/Dy to make everything rapidity density, Dy is the rapidity range
                  double weight = y_weight_local*dy/Dy;
               
                  temp_pol_lambda_theta += dNd2pTdphidy_pol_lambda_theta[m][j][i][k]*weight;
              }
              fphoton_pol_lambda_theta_SpMatrix << scientific << setprecision(6) << setw(16)
                                  << temp_pol_lambda_theta << "  ";            
          }
          fphoton_pol_lambda_theta_SpMatrix << endl;
      }
  }

  
  for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
              for (int k = 0; k < nrapidity; k++){
                  double mll_local = dilepton_QGP_thermal->getDileptonMass(m); //Mll
                  double p_local = dilepton_QGP_thermal->getPhotonp(i); //p_T
                  double y_local = dilepton_QGP_thermal->getPhotonrapidity(k);
                  
                  
                  fphoton_pol_lambda_theta_SpMatrix_dy<< scientific << setprecision(6) << setw(16) 
                              << mll_local << "  "<<p_local<<" "<<y_local<<" "<< dNd2pTd2Mdy_pol_lambda_theta[m][i][k] << "  ";
                  fphoton_pol_lambda_theta_SpMatrix_dy<< endl;
              
              }
      }
  }


  // pT differential, dN/(2pi pTdpT MdM dy) and vn(M, pT)
  for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
          double M_ll_local = dilepton_QGP_thermal->getDileptonMass(m);
          double pT_local = dilepton_QGP_thermal->getPhotonp(i);
         
          fphoton_pol_lambda_theta_Spvn << scientific << setprecision(6) << setw(16)
                          << M_ll_local << "  "<< pT_local << "  " << dNd2pTd2M_pol_lambda_theta[m][i] << "  ";
          fphoton_pol_lambda_theta_Spvn << endl;
          
      }
  }
  

  if (type_str.find("pol_lambda_norm") != std::string::npos)
  {
  // pT integrated
  for (int m = 0; m < nm; m++) {
      double M_ll_local = dilepton_QGP_thermal->getDileptonMass(m);
      // get dN/dMdy from dN/MdMdy
      fphoton_pol_lambda_theta_inte_Spvn << scientific << setprecision(6) << setw(16)
                      << M_ll_local << "  " << dNd2Mdy_pol_lambda_theta[m]/dNd2Mdy_pol_lambda_norm[m] << "  ";
      fphoton_pol_lambda_theta_inte_Spvn << endl;

  }
  }


}





void PhotonEmission::outputPhoton_total_SpMatrix_and_SpvnpT(std::string type_str) {

  ostringstream file_name_label;  
  file_name_label << "all_channel_eq_lo_" << type_str ;

  outputPhoton_SpvnpT(output_path, file_name_label.str(),dNd2pTdphidy_eq_lo, vnpT_cos_eq_lo,vnpT_sin_eq_lo,
    vn_cos_eq_lo, vn_sin_eq_lo,
    dNd2Mdy_eq_lo, dNd2pTd2M_eq_lo,
    dNd2pTd2Mdy_eq_lo,vnMpTy_cos_eq_lo,vnMpTy_sin_eq_lo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_tot_lo_" << type_str ;

  outputPhoton_SpvnpT(output_path, file_name_label.str(),dNd2pTdphidy_tot_lo, vnpT_cos_tot_lo,vnpT_sin_tot_lo,
      vn_cos_tot_lo, vn_sin_tot_lo,
      dNd2Mdy_tot_lo, dNd2pTd2M_tot_lo,
      dNd2pTd2Mdy_tot_lo,vnMpTy_cos_tot_lo,vnMpTy_sin_tot_lo);
  
  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_eq_nlo_" << type_str ;

  outputPhoton_SpvnpT(output_path, file_name_label.str(),dNd2pTdphidy_eq_nlo, vnpT_cos_eq_nlo,vnpT_sin_eq_nlo,
        vn_cos_eq_nlo, vn_sin_eq_nlo,
        dNd2Mdy_eq_nlo, dNd2pTd2M_eq_nlo,
        dNd2pTd2Mdy_eq_nlo,vnMpTy_cos_eq_nlo,vnMpTy_sin_eq_nlo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_tot_nlo_" << type_str ;
    
  outputPhoton_SpvnpT(output_path, file_name_label.str(),dNd2pTdphidy_tot_nlo, vnpT_cos_tot_nlo,vnpT_sin_tot_nlo,
          vn_cos_tot_nlo, vn_sin_tot_nlo,
          dNd2Mdy_tot_nlo, dNd2pTd2M_tot_nlo,
          dNd2pTd2Mdy_tot_nlo,vnMpTy_cos_tot_lo,vnMpTy_sin_tot_nlo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_pol_lambda_theta_lo_" << type_str ;

  outputPhoton_SpvnpT_pol(output_path, file_name_label.str(), dNd2pTdphidy_pol_lambda_theta_lo, dNd2pTd2Mdy_pol_lambda_theta_lo, dNd2pTd2M_pol_lambda_theta_lo, dNd2Mdy_pol_lambda_theta_lo,  dNd2Mdy_pol_lambda_norm_lo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_pol_lambda_phi_lo_" << type_str ;
  
  outputPhoton_SpvnpT_pol(output_path, file_name_label.str(), dNd2pTdphidy_pol_lambda_phi_lo, dNd2pTd2Mdy_pol_lambda_phi_lo, dNd2pTd2M_pol_lambda_phi_lo, dNd2Mdy_pol_lambda_phi_lo,  dNd2Mdy_pol_lambda_norm_lo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_pol_lambda_norm_lo_" << type_str ;
  
  outputPhoton_SpvnpT_pol(output_path,file_name_label.str(), dNd2pTdphidy_pol_lambda_norm_lo, dNd2pTd2Mdy_pol_lambda_norm_lo, dNd2pTd2M_pol_lambda_norm_lo, dNd2Mdy_pol_lambda_norm_lo,  dNd2Mdy_pol_lambda_norm_lo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_pol_lambda_theta_nlo_" << type_str ;


  outputPhoton_SpvnpT_pol(output_path, file_name_label.str(), dNd2pTdphidy_pol_lambda_theta_nlo, dNd2pTd2Mdy_pol_lambda_theta_nlo, dNd2pTd2M_pol_lambda_theta_nlo, dNd2Mdy_pol_lambda_theta_nlo,  dNd2Mdy_pol_lambda_norm_nlo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_pol_lambda_phi_nlo_" << type_str ;
  
  outputPhoton_SpvnpT_pol(output_path, file_name_label.str(), dNd2pTdphidy_pol_lambda_phi_nlo, dNd2pTd2Mdy_pol_lambda_phi_nlo, dNd2pTd2M_pol_lambda_phi_nlo, dNd2Mdy_pol_lambda_phi_nlo,  dNd2Mdy_pol_lambda_norm_nlo);

  file_name_label.str("");
  file_name_label.clear(); 
  file_name_label << "all_channel_pol_lambda_norm_nlo_" << type_str ;
  
  outputPhoton_SpvnpT_pol(output_path, file_name_label.str(), dNd2pTdphidy_pol_lambda_norm_lo, dNd2pTd2Mdy_pol_lambda_norm_lo, dNd2pTd2M_pol_lambda_norm_nlo, dNd2Mdy_pol_lambda_norm_nlo,  dNd2Mdy_pol_lambda_norm_nlo);

  

}


void PhotonEmission::calPhoton_total_Spvn_sum(const PhotonEmission &spvn_tem) {
  // integrate over rapidity, azimuthal angle of momentum, and pT
  // calculate dN/(2pi dydM)

  for (int m = 0; m < nm; m++) {

    dNd2Mdy_eq_lo[m] = 0.0;
    dNd2Mdy_tot_lo[m] = 0.0;
    dNd2Mdy_pol_lambda_theta_lo[m] = 0.0;
    dNd2Mdy_pol_lambda_norm_lo[m] = 0.0;
    dNd2Mdy_pol_lambda_phi_lo[m] = 0.0;

    dNd2Mdy_eq_nlo[m] = 0.0;
    dNd2Mdy_tot_nlo[m] = 0.0;
    dNd2Mdy_pol_lambda_theta_nlo[m] = 0.0;
    dNd2Mdy_pol_lambda_norm_nlo[m] = 0.0;
    dNd2Mdy_pol_lambda_phi_nlo[m] = 0.0;

    for (int order = 0; order < norder; order++) {
      vn_cos_eq_lo[order][m] = 0.0;
      vn_sin_eq_lo[order][m] = 0.0;
      vn_cos_tot_lo[order][m] = 0.0;
      vn_sin_tot_lo[order][m] = 0.0;

      vn_cos_eq_nlo[order][m] = 0.0;
      vn_sin_eq_nlo[order][m] = 0.0;
      vn_cos_tot_nlo[order][m] = 0.0;
      vn_sin_tot_nlo[order][m] = 0.0;
    }

    for (int i = 0; i < np; i++) {
      dNd2pTd2M_eq_lo[m][i] = 0.0;
      dNd2pTd2M_tot_lo[m][i] = 0.0;
      dNd2pTd2M_pol_lambda_theta_lo[m][i] = 0.0;
      dNd2pTd2M_pol_lambda_norm_lo[m][i] = 0.0;
      dNd2pTd2M_pol_lambda_phi_lo[m][i] = 0.0;

      dNd2pTd2M_eq_nlo[m][i] = 0.0;
      dNd2pTd2M_tot_nlo[m][i] = 0.0;
      dNd2pTd2M_pol_lambda_theta_nlo[m][i] = 0.0;
      dNd2pTd2M_pol_lambda_norm_nlo[m][i] = 0.0;
      dNd2pTd2M_pol_lambda_phi_nlo[m][i] = 0.0;

      for (int order = 0; order < norder; order++) {
        vnpT_cos_eq_lo[order][m][i] = 0.0;
        vnpT_cos_tot_lo[order][m][i] = 0.0;
        vnpT_sin_eq_lo[order][m][i] = 0.0;
        vnpT_sin_tot_lo[order][m][i] = 0.0;

        vnpT_cos_eq_nlo[order][m][i] = 0.0;
        vnpT_cos_tot_nlo[order][m][i] = 0.0;
        vnpT_sin_eq_nlo[order][m][i] = 0.0;
        vnpT_sin_tot_nlo[order][m][i] = 0.0;

      }

      for (int j = 0; j < nphi; j++) {
        for (int k = 0; k < nrapidity; k++) {
          dNd2pTdphidy_eq_lo[m][i][j][k] += spvn_tem.dNd2pTdphidy_eq_lo[m][i][j][k];
          dNd2pTdphidy_tot_lo[m][i][j][k] += spvn_tem.dNd2pTdphidy_tot_lo[m][i][j][k];
          dNd2pTdphidy_pol_lambda_theta_lo[m][i][j][k] +=
              spvn_tem.dNd2pTdphidy_pol_lambda_theta_lo[m][i][j][k];
          dNd2pTdphidy_pol_lambda_norm_lo[m][i][j][k] +=
              spvn_tem.dNd2pTdphidy_pol_lambda_norm_lo[m][i][j][k];
          dNd2pTdphidy_pol_lambda_phi_lo[m][i][j][k] +=
              spvn_tem.dNd2pTdphidy_pol_lambda_phi_lo[m][i][j][k];

          dNd2pTdphidy_eq_nlo[m][i][j][k] += spvn_tem.dNd2pTdphidy_eq_nlo[m][i][j][k];
          dNd2pTdphidy_tot_nlo[m][i][j][k] += spvn_tem.dNd2pTdphidy_tot_nlo[m][i][j][k];
          dNd2pTdphidy_pol_lambda_theta_nlo[m][i][j][k] +=
                  spvn_tem.dNd2pTdphidy_pol_lambda_theta_nlo[m][i][j][k];
          dNd2pTdphidy_pol_lambda_norm_nlo[m][i][j][k] +=
                  spvn_tem.dNd2pTdphidy_pol_lambda_norm_nlo[m][i][j][k];
          dNd2pTdphidy_pol_lambda_phi_nlo[m][i][j][k] +=
                  spvn_tem.dNd2pTdphidy_pol_lambda_phi_nlo[m][i][j][k];
        }
      }
    }
  }
  calPhoton_total_Spvn_shell();
  //calPhoton_total_Spvn();
}

// void PhotonEmission::outputPhoton_total_SpMatrix_and_SpvnpT(int hydro_mode) {
//   ostringstream filename_stream_eq_SpMatrix;
//   ostringstream filename_stream_eq_Spvn;
//   ostringstream filename_stream_eq_inte_Spvn;

//   ostringstream filename_stream_eq_TL_SpMatrix;
//   ostringstream filename_stream_eq_TL_Spvn;
//   ostringstream filename_stream_eq_TL_inte_Spvn;

//   ostringstream filename_stream_visc_SpMatrix;
//   ostringstream filename_stream_visc_Spvn;
//   ostringstream filename_stream_visc_inte_Spvn;

//   ostringstream filename_stream_diff_SpMatrix;
//   ostringstream filename_stream_diff_Spvn;
//   ostringstream filename_stream_diff_inte_Spvn;

//   ostringstream filename_stream_tot_SpMatrix;
//   ostringstream filename_stream_tot_Spvn;
//   ostringstream filename_stream_tot_inte_Spvn;

//   ostringstream filename_stream_pol_lambda_theta_SpMatrix;
//   ostringstream filename_stream_pol_lambda_theta_Spvn;
//   ostringstream filename_stream_pol_lambda_theta_inte_Spvn;

//   ostringstream filename_stream_pol_lambda_norm_SpMatrix;
//   ostringstream filename_stream_pol_lambda_norm_Spvn;
//   ostringstream filename_stream_pol_lambda_norm_inte_Spvn;

//   ostringstream filename_stream_pol_lambda_phi_SpMatrix;
//   ostringstream filename_stream_pol_lambda_phi_Spvn;
//   ostringstream filename_stream_pol_lambda_phi_inte_Spvn;

//   // rapdity
//   ostringstream filename_stream_eq_SpMatrix_dy;
//   ostringstream filename_stream_eq_TL_SpMatrix_dy;
//   ostringstream filename_stream_visc_SpMatrix_dy;
//   ostringstream filename_stream_diff_SpMatrix_dy;
//   ostringstream filename_stream_tot_SpMatrix_dy;
//   ostringstream filename_stream_pol_lambda_theta_SpMatrix_dy;
//   ostringstream filename_stream_pol_lambda_norm_SpMatrix_dy;
//   ostringstream filename_stream_pol_lambda_phi_SpMatrix_dy;

//   string filename = " ";
//   bool flag_hydro = paraRdr->getVal("flag_hydro");
//   bool flag_prehydro = paraRdr->getVal("flag_prehydro");

//   double sig_lambda = 0.0;
//   double suppress_order = 0.0;

//   if (hydro_mode == 22 && flag_prehydro) {
//     filename = "Pre_hydro_dilepton";
//     sig_lambda = paraRdr->getVal("sig_lambda");
//     suppress_order = paraRdr->getVal("suppress_order");
//   }
//   if (hydro_mode == 12 && flag_hydro) {
//     filename = "Hydro_dilepton";
//   }

//   if (hydro_mode == -1 && flag_hydro && flag_prehydro) {
//     filename = "Total_dilepton";
//     sig_lambda = paraRdr->getVal("sig_lambda");
//     suppress_order = paraRdr->getVal("suppress_order");
//   }

//   ostringstream stream_sig;
//   stream_sig << fixed << setprecision(3) << sig_lambda;
//   string str_sig = stream_sig.str();

//   string num_sf = to_string(static_cast<int>(round(suppress_order)));

//   string file_end = "_" + str_sig + "_" + num_sf + ".dat";
//   if (filename == "Hydro_dilepton") {
//     file_end = ".dat";
//   }

//   filename_stream_eq_SpMatrix << output_path << filename << "_eq_SpMatrix"
//                               << file_end;
//   filename_stream_eq_Spvn << output_path << filename << "_eq_Spvn" << file_end;
//   filename_stream_eq_inte_Spvn << output_path << filename << "_eq_Spvn_inte"
//                                << file_end;

//   filename_stream_eq_TL_SpMatrix << output_path << filename << "_eq_TL_SpMatrix"
//                                  << file_end;
//   filename_stream_eq_TL_Spvn << output_path << filename << "_eq_TL_Sp"
//                              << file_end;
//   filename_stream_eq_TL_inte_Spvn << output_path << filename << "_eq_TL_Sp_inte"
//                                   << file_end;

//   filename_stream_visc_SpMatrix << output_path << filename << "_visc_SpMatrix"
//                                 << file_end;
//   filename_stream_visc_Spvn << output_path << filename << "_visc_Spvn"
//                             << file_end;
//   filename_stream_visc_inte_Spvn << output_path << filename << "_visc_Spvn_inte"
//                                  << file_end;

//   filename_stream_diff_SpMatrix << output_path << filename << "_diff_SpMatrix"
//                                 << file_end;
//   filename_stream_diff_Spvn << output_path << filename << "_diff_Spvn"
//                             << file_end;
//   filename_stream_diff_inte_Spvn << output_path << filename << "_diff_Spvn_inte"
//                                  << file_end;

//   filename_stream_tot_SpMatrix << output_path << filename << "_tot_SpMatrix"
//                                << file_end;
//   filename_stream_tot_Spvn << output_path << filename << "_tot_Spvn"
//                            << file_end;
//   filename_stream_tot_inte_Spvn << output_path << filename << "_tot_Spvn_inte"
//                                 << file_end;

//   filename_stream_pol_lambda_theta_SpMatrix
//       << output_path << filename << "_pol_lambda_theta_SpMatrix" << file_end;
//   filename_stream_pol_lambda_theta_Spvn << output_path << filename
//                                         << "_pol_lambda_theta_Spvn" << file_end;
//   filename_stream_pol_lambda_theta_inte_Spvn
//       << output_path << filename << "_pol_lambda_theta_Spvn_inte" << file_end;

//   filename_stream_pol_lambda_norm_SpMatrix
//       << output_path << filename << "_pol_lambda_norm_SpMatrix" << file_end;
//   filename_stream_pol_lambda_norm_Spvn << output_path << filename
//                                        << "_pol_lambda_norm_Spvn" << file_end;
//   filename_stream_pol_lambda_norm_inte_Spvn
//       << output_path << filename << "_pol_lambda_norm_Spvn_inte" << file_end;

//   filename_stream_pol_lambda_phi_SpMatrix
//       << output_path << filename << "_pol_lambda_phi_SpMatrix" << file_end;
//   filename_stream_pol_lambda_phi_Spvn << output_path << filename
//                                       << "_pol_lambda_phi_Spvn" << file_end;
//   filename_stream_pol_lambda_phi_inte_Spvn
//       << output_path << filename << "_pol_lambda_phi_Spvn_inte" << file_end;

//   // rapidiy
//   filename_stream_eq_SpMatrix_dy << output_path << filename << "_eq_Spvn_MpTy"
//                                  << file_end;
//   filename_stream_eq_TL_SpMatrix_dy << output_path << filename
//                                     << "_eq_TL_Spvn_MpTy" << file_end;
//   filename_stream_visc_SpMatrix_dy << output_path << filename
//                                    << "_visc_Spvn_MpTy" << file_end;
//   filename_stream_diff_SpMatrix_dy << output_path << filename
//                                    << "_diff_Spvn_MpTy" << file_end;
//   filename_stream_tot_SpMatrix_dy << output_path << filename << "_tot_Spvn_MpTy"
//                                   << file_end;
//   filename_stream_pol_lambda_theta_SpMatrix_dy
//       << output_path << filename << "_pol_lambda_theta_Spvn_MpTy" << file_end;
//   filename_stream_pol_lambda_norm_SpMatrix_dy
//       << output_path << filename << "_pol_lambda_norm_Spvn_MpTy" << file_end;
//   filename_stream_pol_lambda_phi_SpMatrix_dy
//       << output_path << filename << "_pol_lambda_phi_Spvn_MpTy" << file_end;

//   ofstream fphoton_eq_SpMatrix_dy(filename_stream_eq_SpMatrix_dy.str().c_str());
//   ofstream fphoton_eq_TL_SpMatrix_dy(
//       filename_stream_eq_TL_SpMatrix_dy.str().c_str());
//   ofstream fphoton_visc_SpMatrix_dy(
//       filename_stream_visc_SpMatrix_dy.str().c_str());
//   ofstream fphoton_diff_SpMatrix_dy(
//       filename_stream_diff_SpMatrix_dy.str().c_str());
//   ofstream fphoton_tot_SpMatrix_dy(
//       filename_stream_tot_SpMatrix_dy.str().c_str());
//   ofstream fphoton_pol_lambda_theta_SpMatrix_dy(
//       filename_stream_pol_lambda_theta_SpMatrix_dy.str().c_str());
//   ofstream fphoton_pol_lambda_norm_SpMatrix_dy(
//       filename_stream_pol_lambda_norm_SpMatrix_dy.str().c_str());
//   ofstream fphoton_pol_lambda_phi_SpMatrix_dy(
//       filename_stream_pol_lambda_phi_SpMatrix_dy.str().c_str());

//   ofstream fphoton_eq_SpMatrix(filename_stream_eq_SpMatrix.str().c_str());
//   ofstream fphoton_eq_Spvn(filename_stream_eq_Spvn.str().c_str());
//   ofstream fphoton_eq_inte_Spvn(filename_stream_eq_inte_Spvn.str().c_str());

//   ofstream fphoton_eq_TL_SpMatrix(filename_stream_eq_TL_SpMatrix.str().c_str());
//   ofstream fphoton_eq_TL_Spvn(filename_stream_eq_TL_Spvn.str().c_str());
//   ofstream fphoton_eq_TL_inte_Spvn(
//       filename_stream_eq_TL_inte_Spvn.str().c_str());

//   ofstream fphoton_visc_SpMatrix(filename_stream_visc_SpMatrix.str().c_str());
//   ofstream fphoton_visc_Spvn(filename_stream_visc_Spvn.str().c_str());
//   ofstream fphoton_visc_inte_Spvn(filename_stream_visc_inte_Spvn.str().c_str());

//   ofstream fphoton_diff_SpMatrix(filename_stream_diff_SpMatrix.str().c_str());
//   ofstream fphoton_diff_Spvn(filename_stream_diff_Spvn.str().c_str());
//   ofstream fphoton_diff_inte_Spvn(filename_stream_diff_inte_Spvn.str().c_str());

//   ofstream fphoton_tot_SpMatrix(filename_stream_tot_SpMatrix.str().c_str());
//   ofstream fphoton_tot_Spvn(filename_stream_tot_Spvn.str().c_str());
//   ofstream fphoton_tot_inte_Spvn(filename_stream_tot_inte_Spvn.str().c_str());

//   ofstream fphoton_pol_lambda_theta_SpMatrix(
//       filename_stream_pol_lambda_theta_SpMatrix.str().c_str());
//   ofstream fphoton_pol_lambda_theta_Spvn(
//       filename_stream_pol_lambda_theta_Spvn.str().c_str());
//   ofstream fphoton_pol_lambda_theta_inte_Spvn(
//       filename_stream_pol_lambda_theta_inte_Spvn.str().c_str());

//   ofstream fphoton_pol_lambda_norm_SpMatrix(
//       filename_stream_pol_lambda_norm_SpMatrix.str().c_str());
//   ofstream fphoton_pol_lambda_norm_Spvn(
//       filename_stream_pol_lambda_norm_Spvn.str().c_str());
//   ofstream fphoton_pol_lambda_norm_inte_Spvn(
//       filename_stream_pol_lambda_norm_inte_Spvn.str().c_str());

//   ofstream fphoton_pol_lambda_phi_SpMatrix(
//       filename_stream_pol_lambda_phi_SpMatrix.str().c_str());
//   ofstream fphoton_pol_lambda_phi_Spvn(
//       filename_stream_pol_lambda_phi_Spvn.str().c_str());
//   ofstream fphoton_pol_lambda_phi_inte_Spvn(
//       filename_stream_pol_lambda_phi_inte_Spvn.str().c_str());

//   double dy = dilepton_QGP_thermal->get_dy();
//   for (int m = 0; m < nm; m++) {
//     for (int i = 0; i < nphi; i++) {
//       double phi = dilepton_QGP_thermal->getPhotonphi(i);
//       fphoton_eq_SpMatrix << phi << "  ";
//       fphoton_eq_TL_SpMatrix << phi << "  ";
//       fphoton_visc_SpMatrix << phi << "  ";
//       fphoton_diff_SpMatrix << phi << "  ";
//       fphoton_tot_SpMatrix << phi << "  ";
//       fphoton_pol_lambda_theta_SpMatrix << phi << "  ";
//       fphoton_pol_lambda_norm_SpMatrix << phi << "  ";
//       fphoton_pol_lambda_phi_SpMatrix << phi << "  ";

//       for (int j = 0; j < np; j++) {
//         double temp_eq = 0.0;
//         double temp_eqT = 0.0;
//         double temp_eqL = 0.0;
//         double temp_visc = 0.0;
//         double temp_diff = 0.0;
//         double temp_tot = 0.0;
//         double temp_pol_lambda_theta = 0.0;
//         double temp_pol_lambda_norm = 0.0;
//         double temp_pol_lambda_phi = 0.0;

//         for (int k = 0; k < nrapidity; k++) {
//           double y_weight = dilepton_QGP_thermal->getPhoton_yweight(k);
//           // below dy/Dy to make everything rapidity density, Dy is the rapidity
//           // range
//           double weight = y_weight * dy / Dy;
//           temp_eq += dNd2pTdphidy_eq[m][j][i][k] * weight;
//           temp_eqT += dNd2pTdphidy_eqT[m][j][i][k] * weight;
//           temp_eqL += dNd2pTdphidy_eqL[m][j][i][k] * weight;
//           temp_visc += dNd2pTdphidy_visc[m][j][i][k] * weight;
//           temp_diff += dNd2pTdphidy_diff[m][j][i][k] * weight;
//           temp_tot += dNd2pTdphidy_tot[m][j][i][k] * weight;
//           temp_pol_lambda_theta +=
//               dNd2pTdphidy_pol_lambda_theta[m][j][i][k] * weight;
//           temp_pol_lambda_norm +=
//               dNd2pTdphidy_pol_lambda_norm[m][j][i][k] * weight;
//           temp_pol_lambda_phi +=
//               dNd2pTdphidy_pol_lambda_phi[m][j][i][k] * weight;
//         }
//         fphoton_eq_SpMatrix << scientific << setprecision(6) << setw(16)
//                             << temp_eq << "  ";
//         fphoton_eq_TL_SpMatrix << scientific << setprecision(6) << setw(16)
//                                << temp_eqT << "  " << temp_eqL << "  ";
//         fphoton_visc_SpMatrix << scientific << setprecision(6) << setw(16)
//                               << temp_visc << "  ";
//         fphoton_diff_SpMatrix << scientific << setprecision(6) << setw(16)
//                               << temp_diff << "  ";
//         fphoton_tot_SpMatrix << scientific << setprecision(6) << setw(16)
//                              << temp_tot << "  ";
//         fphoton_pol_lambda_theta_SpMatrix << scientific << setprecision(6)
//                                           << setw(16) << temp_pol_lambda_theta
//                                           << "  ";
//         fphoton_pol_lambda_norm_SpMatrix << scientific << setprecision(6)
//                                          << setw(16) << temp_pol_lambda_norm
//                                          << "  ";
//         fphoton_pol_lambda_phi_SpMatrix << scientific << setprecision(6)
//                                         << setw(16) << temp_pol_lambda_phi
//                                         << "  ";
//       }
//       fphoton_eq_SpMatrix << endl;
//       fphoton_eq_TL_SpMatrix << endl;
//       fphoton_visc_SpMatrix << endl;
//       fphoton_diff_SpMatrix << endl;
//       fphoton_tot_SpMatrix << endl;
//       fphoton_pol_lambda_theta_SpMatrix << endl;
//       fphoton_pol_lambda_norm_SpMatrix << endl;
//       fphoton_pol_lambda_phi_SpMatrix << endl;
//     }
//   }

//   for (int m = 0; m < nm; m++) {
//     for (int i = 0; i < np; i++) {
//       for (int k = 0; k < nrapidity; k++) {
//         double mll_local = dilepton_QGP_thermal->getDileptonMass(m); // Mll
//         double p_local = dilepton_QGP_thermal->getPhotonp(i);        // p_T
//         double y_local = dilepton_QGP_thermal->getPhotonrapidity(k);

//         fphoton_eq_SpMatrix_dy << scientific << setprecision(6) << setw(16)
//                                << mll_local << "  " << p_local << " " << y_local
//                                << " " << dNd2pTd2Mdy_eq[m][i][k] << "  ";
//         fphoton_eq_TL_SpMatrix_dy << scientific << setprecision(6) << setw(16)
//                                   << mll_local << "  " << p_local << " "
//                                   << y_local << " " << dNd2pTd2Mdy_eqT[m][i][k]
//                                   << "  " << dNd2pTd2Mdy_eqL[m][i][k] << " ";
//         fphoton_visc_SpMatrix_dy << scientific << setprecision(6) << setw(16)
//                                  << mll_local << "  " << p_local << " "
//                                  << y_local << " " << dNd2pTd2Mdy_visc[m][i][k]
//                                  << "  ";
//         fphoton_diff_SpMatrix_dy << scientific << setprecision(6) << setw(16)
//                                  << mll_local << "  " << p_local << " "
//                                  << y_local << " " << dNd2pTd2Mdy_diff[m][i][k]
//                                  << "  ";
//         fphoton_tot_SpMatrix_dy << scientific << setprecision(6) << setw(16)
//                                 << mll_local << "  " << p_local << " "
//                                 << y_local << " " << dNd2pTd2Mdy_tot[m][i][k]
//                                 << "  ";
//         fphoton_pol_lambda_theta_SpMatrix_dy
//             << scientific << setprecision(6) << setw(16) << mll_local << "  "
//             << p_local << " " << y_local << " "
//             << dNd2pTd2Mdy_pol_lambda_theta[m][i][k] << "  ";
//         fphoton_pol_lambda_norm_SpMatrix_dy
//             << scientific << setprecision(6) << setw(16) << mll_local << "  "
//             << p_local << " " << y_local << " "
//             << dNd2pTd2Mdy_pol_lambda_norm[m][i][k] << "  ";
//         fphoton_pol_lambda_phi_SpMatrix_dy
//             << scientific << setprecision(6) << setw(16) << mll_local << "  "
//             << p_local << " " << y_local << " "
//             << dNd2pTd2Mdy_pol_lambda_phi[m][i][k] << "  ";
//         for (int order = 1; order < norder; order++) {

//           fphoton_eq_SpMatrix_dy << scientific << setprecision(6) << setw(16)
//                                  << order << "  "
//                                  << vnMpTy_cos_eq[order][m][i][k] << " "
//                                  << vnMpTy_sin_eq[order][m][i][k] << " "
//                                  << sqrt(pow(vnMpTy_cos_eq[order][m][i][k], 2) +
//                                          pow(vnMpTy_sin_eq[order][m][i][k], 2))
//                                  << "  ";
//           fphoton_visc_SpMatrix_dy
//               << scientific << setprecision(6) << setw(16) << order << "  "
//               << vnMpTy_cos_visc[order][m][i][k] << " "
//               << vnMpTy_sin_visc[order][m][i][k] << " "
//               << sqrt(pow(vnMpTy_cos_visc[order][m][i][k], 2) +
//                       pow(vnMpTy_sin_visc[order][m][i][k], 2))
//               << "  ";
//           fphoton_diff_SpMatrix_dy
//               << scientific << setprecision(6) << setw(16) << order << "  "
//               << vnMpTy_cos_diff[order][m][i][k] << " "
//               << vnMpTy_sin_diff[order][m][i][k] << " "
//               << sqrt(pow(vnMpTy_cos_diff[order][m][i][k], 2) +
//                       pow(vnMpTy_sin_diff[order][m][i][k], 2))
//               << "  ";

//           fphoton_tot_SpMatrix_dy
//               << scientific << setprecision(6) << setw(16) << order << "  "
//               << vnMpTy_cos_tot[order][m][i][k] << " "
//               << vnMpTy_sin_tot[order][m][i][k] << " "
//               << sqrt(pow(vnMpTy_cos_tot[order][m][i][k], 2) +
//                       pow(vnMpTy_sin_tot[order][m][i][k], 2))
//               << "  ";
//         }

//         fphoton_eq_SpMatrix_dy << endl;
//         fphoton_eq_TL_SpMatrix_dy << endl;
//         fphoton_visc_SpMatrix_dy << endl;
//         fphoton_diff_SpMatrix_dy << endl;
//         fphoton_tot_SpMatrix_dy << endl;
//         fphoton_pol_lambda_theta_SpMatrix_dy << endl;
//         fphoton_pol_lambda_norm_SpMatrix_dy << endl;
//         fphoton_pol_lambda_phi_SpMatrix_dy << endl;
//       }
//     }
//   }

//   // pT differential, dN/(2pi pTdpT MdM dy) and vn(M, pT)
//   for (int m = 0; m < nm; m++) {
//     for (int i = 0; i < np; i++) {
//       double M_ll = dilepton_QGP_thermal->getDileptonMass(m);
//       double pT = dilepton_QGP_thermal->getPhotonp(i);
//       fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                       << "  " << pT << "  " << dNd2pTd2M_eq[m][i] << "  ";
//       fphoton_eq_TL_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                          << "  " << pT << "  " << dNd2pTd2M_eqT[m][i] << "  "
//                          << dNd2pTd2M_eqL[m][i];
//       fphoton_visc_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                         << "  " << pT << "  " << dNd2pTd2M_visc[m][i] << "  ";
//       fphoton_diff_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                         << "  " << pT << "  " << dNd2pTd2M_diff[m][i] << "  ";
//       fphoton_tot_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                        << "  " << pT << "  " << dNd2pTd2M_tot[m][i] << "  ";

//       fphoton_pol_lambda_theta_Spvn << scientific << setprecision(6) << setw(16)
//                                     << M_ll << "  " << pT << "  "
//                                     << dNd2pTd2M_pol_lambda_theta[m][i] << "  ";

//       fphoton_pol_lambda_norm_Spvn << scientific << setprecision(6) << setw(16)
//                                    << M_ll << "  " << pT << "  "
//                                    << dNd2pTd2M_pol_lambda_norm[m][i] << "  ";

//       fphoton_pol_lambda_phi_Spvn << scientific << setprecision(6) << setw(16)
//                                   << M_ll << "  " << pT << "  "
//                                   << dNd2pTd2M_pol_lambda_phi[m][i] << "  ";

//       for (int order = 1; order < norder; order++) {
//         fphoton_eq_Spvn << scientific << setprecision(6) << setw(16) << order
//                         << "   " << vnpT_cos_eq[order][m][i] << "  "
//                         << vnpT_sin_eq[order][m][i] << "  "
//                         << sqrt(pow(vnpT_cos_eq[order][m][i], 2) +
//                                 pow(vnpT_sin_eq[order][m][i], 2))
//                         << "  ";
//         fphoton_visc_Spvn << scientific << setprecision(6) << setw(16) << order
//                           << "   " << vnpT_cos_visc[order][m][i] << "  "
//                           << vnpT_sin_visc[order][m][i] << "  "
//                           << sqrt(pow(vnpT_cos_visc[order][m][i], 2) +
//                                   pow(vnpT_sin_visc[order][m][i], 2))
//                           << "  ";
//         fphoton_diff_Spvn << scientific << setprecision(6) << setw(16) << order
//                           << "   " << vnpT_cos_diff[order][m][i] << "  "
//                           << vnpT_sin_diff[order][m][i] << "  "
//                           << sqrt(pow(vnpT_cos_diff[order][m][i], 2) +
//                                   pow(vnpT_sin_diff[order][m][i], 2))
//                           << "  ";
//         fphoton_tot_Spvn << scientific << setprecision(6) << setw(16) << order
//                          << "   " << vnpT_cos_tot[order][m][i] << "  "
//                          << vnpT_sin_tot[order][m][i] << "  "
//                          << sqrt(pow(vnpT_cos_tot[order][m][i], 2) +
//                                  pow(vnpT_sin_tot[order][m][i], 2))
//                          << "  ";
//       }
//       fphoton_eq_Spvn << endl;
//       fphoton_eq_TL_Spvn << endl;
//       fphoton_visc_Spvn << endl;
//       fphoton_diff_Spvn << endl;
//       fphoton_tot_Spvn << endl;
//       fphoton_pol_lambda_theta_Spvn << endl;
//       fphoton_pol_lambda_norm_Spvn << endl;
//       fphoton_pol_lambda_phi_Spvn << endl;
//     }
//   }

//   // pT integrated
//   for (int m = 0; m < nm; m++) {
//     double M_ll = dilepton_QGP_thermal->getDileptonMass(m);
//     // get dN/dMdy from dN/MdMdy
//     fphoton_eq_inte_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                          << "  " << dNd2Mdy_eq[m] << "  ";
//     fphoton_eq_TL_inte_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                             << "  " << dNd2Mdy_eqT[m] << "  " << dNd2Mdy_eqL[m];
//     fphoton_visc_inte_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                            << "  " << dNd2Mdy_visc[m] << "  ";
//     fphoton_diff_inte_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                            << "  " << dNd2Mdy_diff[m] << "  ";
//     fphoton_tot_inte_Spvn << scientific << setprecision(6) << setw(16) << M_ll
//                           << "  " << dNd2Mdy_tot[m] << "  ";

//     fphoton_pol_lambda_theta_inte_Spvn
//         << scientific << setprecision(6) << setw(16) << M_ll << "  "
//         << dNd2Mdy_pol_lambda_theta[m] / dNd2Mdy_pol_lambda_norm[m] << "  ";
//     fphoton_pol_lambda_phi_inte_Spvn
//         << scientific << setprecision(6) << setw(16) << M_ll << "  "
//         << dNd2Mdy_pol_lambda_phi[m] / dNd2Mdy_pol_lambda_norm[m] << "  ";

//     fphoton_pol_lambda_norm_inte_Spvn << scientific << setprecision(6)
//                                       << setw(16) << M_ll << "  "
//                                       << dNd2Mdy_pol_lambda_norm[m] << "  ";

//     for (int order = 0; order < norder; order++) {
//       fphoton_eq_inte_Spvn << scientific << setprecision(6) << setw(16) << order
//                            << "   " << vn_cos_eq[order][m] << "   "
//                            << vn_sin_eq[order][m] << "   "
//                            << sqrt(pow(vn_cos_eq[order][m], 2) +
//                                    pow(vn_sin_eq[order][m], 2))
//                            << "  ";
//       fphoton_visc_inte_Spvn
//           << scientific << setprecision(6) << setw(16) << order << "   "
//           << vn_cos_visc[order][m] << "   " << vn_sin_visc[order][m] << "   "
//           << sqrt(pow(vn_cos_visc[order][m], 2) + pow(vn_sin_visc[order][m], 2))
//           << "  ";
//       fphoton_diff_inte_Spvn
//           << scientific << setprecision(6) << setw(16) << order << "   "
//           << vn_cos_diff[order][m] << "   " << vn_sin_diff[order][m] << "   "
//           << sqrt(pow(vn_cos_diff[order][m], 2) + pow(vn_sin_diff[order][m], 2))
//           << "  ";
//       fphoton_tot_inte_Spvn
//           << scientific << setprecision(6) << setw(16) << order << "   "
//           << vn_cos_tot[order][m] << "   " << vn_sin_tot[order][m] << "   "
//           << sqrt(pow(vn_cos_tot[order][m], 2) + pow(vn_sin_tot[order][m], 2))
//           << "  ";
//     }
//     fphoton_eq_inte_Spvn << endl;
//     fphoton_eq_TL_inte_Spvn << endl;
//     fphoton_visc_inte_Spvn << endl;
//     fphoton_diff_inte_Spvn << endl;
//     fphoton_tot_inte_Spvn << endl;
//     fphoton_pol_lambda_theta_inte_Spvn << endl;
//     fphoton_pol_lambda_norm_inte_Spvn << endl;
//     fphoton_pol_lambda_phi_inte_Spvn << endl;
//   }
// }
