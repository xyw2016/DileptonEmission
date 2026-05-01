/////////////////////////////////////////////////////////////////////////
//  To do in the future:
//      change the integration routines into gaussian integration
/////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
//#include <omp.h>
#ifndef _OPENMP
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#else
#include <omp.h>
#endif

#include "Arsenal.h"
#include "ParameterReader.h"
#include "QGP_NLO.h"
#include "Table2D.h"
#include "ThermalPhoton.h"
#include "data_struct.h"
#include "gauss_quadrature.h"
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
using ARSENAL::logarithmic_mass_grid;
using TENSORTRANSFORM::boost_matrix;


using PhysConsts::me;

ThermalPhoton::ThermalPhoton(std::shared_ptr<ParameterReader> paraRdr_in,
                             std::string emissionProcess)
    : grid_T(), grid_L() {

  paraRdr = paraRdr_in;
  emissionProcess_name = emissionProcess;

    // omp parameters
    CORES = 1;

    #ifdef _OPENMP
      CORES = omp_get_max_threads();
    #endif
  //std::cout << " wxy code: "<< CORES <<std::endl; 
  neta = paraRdr->getVal("neta");
  nm = paraRdr->getVal("nm");
  np = paraRdr->getVal("np");
  nphi = paraRdr->getVal("nphi");
  nrapidity = paraRdr->getVal("nrapidity");
  norder = paraRdr->getVal("norder");
  rate_path_ = "ph_rates/";
        
 


  CS_frame = paraRdr->getVal("CS_frame");

  bRateTable_ = false;
  bShearVisCorr_ = false;
  bBulkVisCorr_ = false;
  bDiffusionCorr_ = false;

  turn_on_muB_ = static_cast<int>(paraRdr->getVal("turn_on_muB", 1));
  include_diff_deltaf = paraRdr->getVal("include_baryondiff_deltaf");
  include_visc_deltaf = paraRdr->getVal("include_shearvisc_deltaf");

  alpha_s = paraRdr->getVal("alpha_s");

  // if muB is off, no need to include diffusion correction
  if (turn_on_muB_ == 0)
    include_diff_deltaf = 0;

  // initial variables for photon spectra
  double p_i = paraRdr->getVal("photon_q_i");
  double p_f = paraRdr->getVal("photon_q_f");
  double phi_i = paraRdr->getVal("photon_phi_q_i");
  double phi_f = paraRdr->getVal("photon_phi_q_f");
  double y_i = paraRdr->getVal("photon_y_i");
  double y_f = paraRdr->getVal("photon_y_f");

  double m_i = paraRdr->getVal("dilepton_mass_i");
  double m_f = paraRdr->getVal("dilepton_mass_f");

  // Choose between equal steps and logarithmically spaced mass grid
  use_logarithmic_mass_grid = paraRdr->getVal("use_logarithmic_mass_grid");

  p = new double[np];
  p_weight = new double[np];
  phi = new double[nphi];
  phi_weight = new double[nphi];

  gauss_quadrature(np, 1, 0.0, 0.0, p_i, p_f, p, p_weight);
  gauss_quadrature(nphi, 1, 0.0, 0.0, phi_i, phi_f, phi, phi_weight);

  // dilepton rapidity
  y_weight = new double[nrapidity];

  if (nrapidity > 1) {
    Dy = y_f - y_i;
    dy = Dy / (nrapidity - 1);
    trapezoidal_weights(nrapidity, y_weight);
  } else {
    Dy = 1.0;
    dy = 1.0;
    y_weight[0] = 1.0;
  }

  y.resize(nrapidity, 0);

  for (int i = 0; i < nrapidity; i++) {
    y[i] = y_i + i * dy;
  }

  // dilepton invariant mass
  if (use_logarithmic_mass_grid) {
    M = logarithmic_mass_grid(m_i, m_f, nm);
  } else {
    M.resize(nm, 0);
    dM = (m_f - m_i) / (nm - 1);
    for (int i = 0; i < nm; i++) {
      M[i] = m_i + i * dM;
    }
  }

  dNd2pTdphidy_eq_all = (double*)calloc(CORES * nrapidity*np*nphi*nm,
  sizeof(double)); 
  dNd2pTdphidy_eqT_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double)); 
  dNd2pTdphidy_eqL_all  =  (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
  dNd2pTdphidy_visc_all = (double*)calloc(CORES * nrapidity*np*nphi*nm,sizeof(double)); 
  dNd2pTdphidy_diff_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double)); 
  dNd2pTdphidy_tot_all  =  (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double)); 
  dNd2pTdphidy_pol_lambda_theta_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double)); 
  dNd2pTdphidy_pol_lambda_norm_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double)); 
  dNd2pTdphidy_pol_lambda_phi_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));


    // dN/MdMdPT2dphidy
  dNd2pTdphidy_eq = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
  dNd2pTdphidy_eqT = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
  dNd2pTdphidy_eqL = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
  dNd2pTdphidy_visc = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
  dNd2pTdphidy_diff = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
  dNd2pTdphidy_tot = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
  dNd2pTdphidy_pol_lambda_theta = createA4DMatrix(nm, np, nphi, nrapidity,
  0.); dNd2pTdphidy_pol_lambda_norm = createA4DMatrix(nm, np, nphi,
  nrapidity, 0.); dNd2pTdphidy_pol_lambda_phi = createA4DMatrix(nm, np, nphi,
  nrapidity, 0.);


    // // dN/dMdPT2dphidy
  dNd2pTd2M_eq = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_eqT = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_eqL = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_visc = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_diff = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_tot = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_pol_lambda_theta = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_pol_lambda_norm = createA2DMatrix(nm, np, 0);
  dNd2pTd2M_pol_lambda_phi = createA2DMatrix(nm, np, 0);

  dNd2pTd2Mdy_eq = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_eqT = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_eqL = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_visc = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_diff = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_tot = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_pol_lambda_theta = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_pol_lambda_norm = createA3DMatrix(nm, np, nrapidity, 0);
  dNd2pTd2Mdy_pol_lambda_phi = createA3DMatrix(nm, np, nrapidity, 0);

  dNd2Mdy_eq = createA1DMatrix(nm, 0); 
  dNd2Mdy_eqT = createA1DMatrix(nm, 0); 
  dNd2Mdy_eqL = createA1DMatrix(nm, 0); 
  dNd2Mdy_visc = createA1DMatrix(nm, 0); 
  dNd2Mdy_diff = createA1DMatrix(nm, 0); 
  dNd2Mdy_tot = createA1DMatrix(nm, 0); 
  dNd2Mdy_pol_lambda_theta = createA1DMatrix(nm, 0); 
  dNd2Mdy_pol_lambda_norm = createA1DMatrix(nm, 0); 
  dNd2Mdy_pol_lambda_phi = createA1DMatrix(nm, 0); 

  vnpT_cos_eq = createA3DMatrix(norder, nm, np, 0.);
  vnpT_sin_eq = createA3DMatrix(norder, nm, np, 0.);
  vnpT_cos_eqT = createA3DMatrix(norder, nm, np, 0.);
  vnpT_sin_eqT = createA3DMatrix(norder, nm, np, 0.);
  vnpT_cos_eqL = createA3DMatrix(norder, nm, np, 0.);
  vnpT_sin_eqL = createA3DMatrix(norder, nm, np, 0.);

  vnpT_cos_visc = createA3DMatrix(norder, nm, np, 0.);
  vnpT_sin_visc = createA3DMatrix(norder, nm, np, 0.);
  vnpT_cos_diff = createA3DMatrix(norder, nm, np, 0.);
  vnpT_sin_diff = createA3DMatrix(norder, nm, np, 0.);
  vnpT_cos_tot = createA3DMatrix(norder, nm, np, 0.);
  vnpT_sin_tot = createA3DMatrix(norder, nm, np, 0.);

  vnMpTy_cos_eq = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_sin_eq = createA4DMatrix(norder, nm, np, nrapidity, 0.);

  vnMpTy_cos_eqT = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_sin_eqT = createA4DMatrix(norder, nm, np, nrapidity, 0.);

  vnMpTy_cos_eqL = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_sin_eqL = createA4DMatrix(norder, nm, np, nrapidity, 0.);

  vnMpTy_cos_visc = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_sin_visc = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_cos_diff = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_sin_diff = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_cos_tot = createA4DMatrix(norder, nm, np, nrapidity, 0.);
  vnMpTy_sin_tot = createA4DMatrix(norder, nm, np, nrapidity, 0.);

  vn_cos_eq = createA2DMatrix(norder, nm, 0.);
  vn_sin_eq = createA2DMatrix(norder, nm, 0.);
  
  vn_cos_eqT = createA2DMatrix(norder, nm, 0.);
  vn_sin_eqT = createA2DMatrix(norder, nm, 0.);

  vn_cos_eqL = createA2DMatrix(norder, nm, 0.);
  vn_sin_eqL = createA2DMatrix(norder, nm, 0.);

  vn_cos_visc = createA2DMatrix(norder, nm, 0.);
  vn_sin_visc = createA2DMatrix(norder, nm, 0.);
  vn_cos_diff = createA2DMatrix(norder, nm, 0.);
  vn_sin_diff = createA2DMatrix(norder, nm, 0.);
  vn_cos_tot = createA2DMatrix(norder, nm, 0.);
  vn_sin_tot = createA2DMatrix(norder, nm, 0.);




  int differential_flag = paraRdr->getVal("differential_flag");

  if (differential_flag == 1 or differential_flag > 10) {
    nTcut = paraRdr->getVal("nTcut");
    n_tau_cut = paraRdr->getVal("n_tau_cut");

    dNd2pTdphidydTdtau_eq =
        createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidydTdtau_visc =
        createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidydTdtau_diff =
        createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidydTdtau_tot =
        createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);

    dNpTdpTdydTdtau_eq = createA4DMatrix(nTcut, n_tau_cut, nm, np, 0.);
    dNpTdpTdydTdtau_visc = createA4DMatrix(nTcut, n_tau_cut, nm, np, 0.);
    dNpTdpTdydTdtau_diff = createA4DMatrix(nTcut, n_tau_cut, nm, np, 0.);
    dNpTdpTdydTdtau_tot = createA4DMatrix(nTcut, n_tau_cut, nm, np, 0.);

    dNdydTdtau_eq = createA3DMatrix(nTcut, n_tau_cut, nm, 0.);
    dNdydTdtau_visc = createA3DMatrix(nTcut, n_tau_cut, nm, 0.);
    dNdydTdtau_diff = createA3DMatrix(nTcut, n_tau_cut, nm, 0.);
    dNdydTdtau_tot = createA3DMatrix(nTcut, n_tau_cut, nm, 0.);

    vndTdtau_cos_eq = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    vndTdtau_sin_eq = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    vndTdtau_cos_visc = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    vndTdtau_sin_visc = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    vndTdtau_cos_diff = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    vndTdtau_sin_diff = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    vndTdtau_cos_tot = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    vndTdtau_sin_tot = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);

    // in temperature or proper time
    dNpTdpTdydT_eq = createA3DMatrix(nTcut, nm, np, 0.);
    dNpTdpTdydT_visc = createA3DMatrix(nTcut, nm, np, 0.);
    dNpTdpTdydT_diff = createA3DMatrix(nTcut, nm, np, 0.);
    dNpTdpTdydT_tot = createA3DMatrix(nTcut, nm, np, 0.);

    dNpTdpTdydtau_eq = createA3DMatrix(n_tau_cut, nm, np, 0.);
    dNpTdpTdydtau_visc = createA3DMatrix(n_tau_cut, nm, np, 0.);
    dNpTdpTdydtau_diff = createA3DMatrix(n_tau_cut, nm, np, 0.);
    dNpTdpTdydtau_tot = createA3DMatrix(n_tau_cut, nm, np, 0.);

    dNdydT_eq = createA2DMatrix(nTcut, nm, 0.);
    dNdydT_visc = createA2DMatrix(nTcut, nm, 0.);
    dNdydT_diff = createA2DMatrix(nTcut, nm, 0.);
    dNdydT_tot = createA2DMatrix(nTcut, nm, 0.);

    dNdydtau_eq = createA2DMatrix(n_tau_cut, nm, 0.);
    dNdydtau_visc = createA2DMatrix(n_tau_cut, nm, 0.);
    dNdydtau_diff = createA2DMatrix(n_tau_cut, nm, 0.);
    dNdydtau_tot = createA2DMatrix(n_tau_cut, nm, 0.);
  }
}

ThermalPhoton::~ThermalPhoton() {

  delete[] p;
  delete[] p_weight;
  delete[] phi;
  delete[] phi_weight;
  delete[] y_weight;
   
  //xyw
  free(dNd2pTdphidy_eq_all);
  free(dNd2pTdphidy_eqT_all);
  free(dNd2pTdphidy_eqL_all);
  free(dNd2pTdphidy_visc_all);
  free(dNd2pTdphidy_diff_all);
  free(dNd2pTdphidy_tot_all);
  free(dNd2pTdphidy_pol_lambda_theta_all);
  free(dNd2pTdphidy_pol_lambda_norm_all);
  free(dNd2pTdphidy_pol_lambda_phi_all);

    deleteA4DMatrix(dNd2pTdphidy_eq, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_eqT, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_eqL, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_visc, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_diff, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_tot, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_pol_lambda_theta, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_pol_lambda_norm, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_pol_lambda_phi, nm, np, nphi);



  deleteA2DMatrix(dNd2pTd2M_eq, nm);
  deleteA2DMatrix(dNd2pTd2M_eqT, nm);
  deleteA2DMatrix(dNd2pTd2M_eqL, nm);
  deleteA2DMatrix(dNd2pTd2M_visc, nm);
  deleteA2DMatrix(dNd2pTd2M_diff, nm);
  deleteA2DMatrix(dNd2pTd2M_tot, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_theta, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_norm, nm);
  deleteA2DMatrix(dNd2pTd2M_pol_lambda_phi, nm);

  deleteA3DMatrix(dNd2pTd2Mdy_eq, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_eqT, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_eqL, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_visc, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_diff, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_tot, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_theta, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_norm, nm, np);
  deleteA3DMatrix(dNd2pTd2Mdy_pol_lambda_phi, nm, np);

  deleteA3DMatrix(vnpT_cos_eq, norder, nm);
  deleteA3DMatrix(vnpT_sin_eq, norder, nm);
  deleteA3DMatrix(vnpT_cos_eqT, norder, nm);
  deleteA3DMatrix(vnpT_sin_eqT, norder, nm);
  deleteA3DMatrix(vnpT_cos_eqL, norder, nm);
  deleteA3DMatrix(vnpT_sin_eqL, norder, nm);
  deleteA3DMatrix(vnpT_cos_visc, norder, nm);
  deleteA3DMatrix(vnpT_sin_visc, norder, nm);
  deleteA3DMatrix(vnpT_cos_diff, norder, nm);
  deleteA3DMatrix(vnpT_sin_diff, norder, nm);
  deleteA3DMatrix(vnpT_cos_tot, norder, nm);
  deleteA3DMatrix(vnpT_sin_tot, norder, nm);

  deleteA4DMatrix(vnMpTy_cos_eq, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_eq, norder, nm, np);

  deleteA4DMatrix(vnMpTy_cos_eqT, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_eqT, norder, nm, np);

  deleteA4DMatrix(vnMpTy_cos_eqL, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_eqL, norder, nm, np);

  deleteA4DMatrix(vnMpTy_cos_visc, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_visc, norder, nm, np);
  deleteA4DMatrix(vnMpTy_cos_diff, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_diff, norder, nm, np);
  deleteA4DMatrix(vnMpTy_cos_tot, norder, nm, np);
  deleteA4DMatrix(vnMpTy_sin_tot, norder, nm, np);

  deleteA2DMatrix(vn_cos_eq, norder);
  deleteA2DMatrix(vn_sin_eq, norder);

  deleteA2DMatrix(vn_cos_eqT, norder);
  deleteA2DMatrix(vn_sin_eqT, norder);

  deleteA2DMatrix(vn_cos_eqL, norder);
  deleteA2DMatrix(vn_sin_eqL, norder);

  deleteA2DMatrix(vn_cos_visc, norder);
  deleteA2DMatrix(vn_sin_visc, norder);
  deleteA2DMatrix(vn_cos_diff, norder);
  deleteA2DMatrix(vn_sin_diff, norder);
  deleteA2DMatrix(vn_cos_tot, norder);
  deleteA2DMatrix(vn_sin_tot, norder);


  deleteA1DMatrix(dNd2Mdy_eq); 
  deleteA1DMatrix(dNd2Mdy_eqT); 
  deleteA1DMatrix(dNd2Mdy_eqL); 
  deleteA1DMatrix(dNd2Mdy_visc); 
  deleteA1DMatrix(dNd2Mdy_diff); 
  deleteA1DMatrix(dNd2Mdy_tot); 
  deleteA1DMatrix(dNd2Mdy_pol_lambda_theta);
  deleteA1DMatrix(dNd2Mdy_pol_lambda_norm); 
  deleteA1DMatrix(dNd2Mdy_pol_lambda_phi); 



  int differential_flag = paraRdr->getVal("differential_flag");
  if (differential_flag == 1 or differential_flag > 10) {
    deleteA3DMatrix(dNdydTdtau_eq, nTcut, n_tau_cut);
    deleteA3DMatrix(dNdydTdtau_visc, nTcut, n_tau_cut);
    deleteA3DMatrix(dNdydTdtau_diff, nTcut, n_tau_cut);
    deleteA3DMatrix(dNdydTdtau_tot, nTcut, n_tau_cut);

    deleteA4DMatrix(vndTdtau_cos_eq, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(vndTdtau_sin_eq, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(vndTdtau_cos_visc, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(vndTdtau_sin_visc, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(vndTdtau_cos_diff, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(vndTdtau_sin_diff, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(vndTdtau_cos_tot, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(vndTdtau_sin_tot, nTcut, n_tau_cut, nm);

    deleteA6DMatrix(dNd2pTdphidydTdtau_eq, nTcut, n_tau_cut, nm, np, nphi);
    deleteA6DMatrix(dNd2pTdphidydTdtau_visc, nTcut, n_tau_cut, nm, np, nphi);
    deleteA6DMatrix(dNd2pTdphidydTdtau_diff, nTcut, n_tau_cut, nm, np, nphi);
    deleteA6DMatrix(dNd2pTdphidydTdtau_tot, nTcut, n_tau_cut, nm, np, nphi);

    deleteA4DMatrix(dNpTdpTdydTdtau_eq, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(dNpTdpTdydTdtau_visc, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(dNpTdpTdydTdtau_diff, nTcut, n_tau_cut, nm);
    deleteA4DMatrix(dNpTdpTdydTdtau_tot, nTcut, n_tau_cut, nm);

    deleteA3DMatrix(dNpTdpTdydT_eq, nTcut, nm);
    deleteA3DMatrix(dNpTdpTdydT_visc, nTcut, nm);
    deleteA3DMatrix(dNpTdpTdydT_diff, nTcut, nm);
    deleteA3DMatrix(dNpTdpTdydT_tot, nTcut, nm);

    deleteA3DMatrix(dNpTdpTdydtau_eq, n_tau_cut, nm);
    deleteA3DMatrix(dNpTdpTdydtau_visc, n_tau_cut, nm);
    deleteA3DMatrix(dNpTdpTdydtau_diff, n_tau_cut, nm);
    deleteA3DMatrix(dNpTdpTdydtau_tot, n_tau_cut, nm);

    deleteA2DMatrix(dNdydT_eq, nTcut);
    deleteA2DMatrix(dNdydT_visc, nTcut);
    deleteA2DMatrix(dNdydT_diff, nTcut);
    deleteA2DMatrix(dNdydT_tot, nTcut);

    deleteA2DMatrix(dNdydtau_eq, n_tau_cut);
    deleteA2DMatrix(dNdydtau_visc, n_tau_cut);
    deleteA2DMatrix(dNdydtau_diff, n_tau_cut);
    deleteA2DMatrix(dNdydtau_tot, n_tau_cut);
  }
  deleteA4DMatrix(grid_T.F, grid_T.nx, grid_T.ny, grid_T.nz);
  deleteA4DMatrix(grid_L.F, grid_L.nx, grid_L.ny, grid_L.nz);
}

void ThermalPhoton::readEmissionrateFromFile(bool bRateTable) {

  cout << "----------------------------------------" << endl;
  cout << "-- Read in emission rate table:" << endl;
  cout << "----------------------------------------" << endl;

  bRateTable_ = bRateTable;

  // read in equilibrium rate
  ostringstream eqrate_filename_stream;
  eqrate_filename_stream << rate_path_ << "rate_" << emissionProcess_name
                         << "_eqrate.dat";

  string fname = eqrate_filename_stream.str();

  if (emissionProcess_name == "QGP_NLO_total" ){

  

  initialize(fname, a_list, B_list, M_list, k_list, rhoT_list, rhoL_list);

  if (a_list.empty() || B_list.empty() || M_list.empty() || k_list.empty() ||
      rhoT_list.empty()) {
    // Handle the error here

    printf("Error: the vectors are still empty. Emission rate table was not "
           "read properly.\n");
  }

  // Update the tables
  grid_T = Table(a_list, B_list, M_list, k_list, rhoT_list);
  grid_L = Table(a_list, B_list, M_list, k_list, rhoL_list);

  // print some details:
  cout << "-> boundaries of the table ..." << endl;
  cout << " min alpha: " << grid_T.x_min << " , max alpha: " << grid_T.x_max
       << endl;
  cout << " min muB: " << grid_T.y_min << " , max muB: " << grid_T.y_max
       << endl;
  cout << " min M: " << grid_T.z_min << " , max M: " << grid_T.z_max << endl;
  cout << " min k: " << grid_T.w_min << " , max k: " << grid_T.w_max << endl;
  cout << " (muB, M and k are given in units of T!)" << endl << endl;
  }
}

void ThermalPhoton::analyticRates(double T, double muB, vector<double> &Eq,
                                  double *M_ll, std::vector<double> &eqrate_ptr,
                                  int nm, int np, int nphi, int nrapidity) {
  for (unsigned int i = 0; i < eqrate_ptr.size(); i++) {
    eqrate_ptr[i] = 1e-16;
  }
}

void ThermalPhoton::analyticRatesShearVis(double T, vector<double> &Eq,
                                          double *M_ll,
                                          std::vector<double> &visrate_ptr) {
  for (unsigned int i = 0; i < visrate_ptr.size(); i++) {
    visrate_ptr[i] = 0.;
  }
}

void ThermalPhoton::analyticRatesBulkVis(double T, vector<double> &Eq,
                                         double *M_ll,
                                         std::vector<double> &bulkvis_ptr) {
  for (unsigned int i = 0; i < bulkvis_ptr.size(); i++) {
    bulkvis_ptr[i] = 0.;
  }
}

void ThermalPhoton::FiniteBaryonRates(
    double T, double muB, double inv_eplusp, double rhoB_over_eplusp, double Eq,
    double M_ll, double &eqrate_ptr, double &eqrateT_ptr, double &eqrateL_ptr,
    double &viscrate_ptr, double &diffrate_ptr, int include_visc_deltaf,
    int include_diff_deltaf) {
        eqrate_ptr = 1.e-16;
        eqrateT_ptr = 1.e-16;
        eqrateL_ptr = 1.e-16;
        viscrate_ptr = 0.;
        diffrate_ptr = 0. ;
        //em_diffrate = 0.;
}

void ThermalPhoton::getPhotonemissionRate(
    double &Eq, double &M_ll, double &pi_factor, double &bulkPi_factor,
    double &diff_factor, double &T, double &muB, double &inv_eplusp,
    double &rhoB_over_eplusp, double& em_eqrate ,double& em_eqrateT,
    double &em_eqrateL, double& em_visrate,double& em_bulkvis, double& em_diffrate) {

  if (emissionProcess_name == "QGP_NLO_total") {
    // interpolate NLO equilibrium rate
    double k = sqrt(Eq * Eq - M_ll * M_ll);
    NLO_rate(grid_T, grid_L, Eq, k, alpha_s, muB, T, me, em_eqrate,
        em_eqrateT, em_eqrateL);
    em_diffrate = 0.;

  } else if (emissionProcess_name == "QGP_LO_total") {
    // use LO analytical form
    FiniteBaryonRates(T, muB, inv_eplusp, rhoB_over_eplusp, Eq, M_ll,
                      em_eqrate, em_eqrateT, em_eqrateL, em_visrate,
                      em_diffrate, include_visc_deltaf, include_diff_deltaf);
  } else if (emissionProcess_name == "HadronGas_rho") {
    double k = sqrt(Eq * Eq - M_ll * M_ll);
    getRateFromTable( Eq, T,  k, M_ll,
                      em_eqrate, em_eqrateT, em_eqrateL);
  } 
  else if (emissionProcess_name == "HadronGas_4pi") {
    double k = sqrt(Eq * Eq - M_ll * M_ll);
    getRateFromTable( Eq, T,  k, M_ll,
                      em_eqrate, em_eqrateT, em_eqrateL);
  }  
  else {
    std::cout << " error!: emissionProcess! " << std::endl;
  }
  em_visrate = pi_factor * em_visrate;
  em_diffrate = diff_factor * em_diffrate;
}



void ThermalPhoton::calThermalPhotonemission_3d(
    double (&p_lab_Min)[4], double (&flow_u_mu_Min)[4], double &Eq,
    double &M_ll, double &pi_zz, double &bulkPi, double &diff_factor, double &T,
    double &muB, double &inv_eplusp, double &rhoB_over_eplusp, double &volume,
    double &fraction, long& n, int& i3) {

  const double volfrac = volume * fraction;
 
  
  double em_eqrate_local = 0.0;
  double em_eqrateT_local= 0.0;
  double em_eqrateL_local= 0.0;
  double em_visrate_local= 0.0;
  double em_bulkvis_local= 0.0;
  double em_diffrate_local= 0.0;
  getPhotonemissionRate(Eq, M_ll, pi_zz, bulkPi, diff_factor, T, muB,
                        inv_eplusp, rhoB_over_eplusp,em_eqrate_local,em_eqrateT_local,
                        em_eqrateL_local,em_visrate_local,em_bulkvis_local,em_diffrate_local);
    

  double temp_eq_sum   = em_eqrate_local*volfrac;
  double temp_eqT_sum  = em_eqrateT_local*volfrac;
  double temp_eqL_sum  = em_eqrateL_local*volfrac;
  double temp_visc_sum = em_visrate_local*volfrac;
  double temp_bulkvis_sum = em_bulkvis_local*volfrac;
  double temp_diff_sum = em_diffrate_local*volfrac;

  
    
    
  dNd2pTdphidy_eq_all[n + CORES * i3] += temp_eq_sum;
  dNd2pTdphidy_eqT_all[n + CORES * i3] += temp_eqT_sum;
  dNd2pTdphidy_eqL_all[n + CORES * i3] += temp_eqL_sum;
  dNd2pTdphidy_visc_all[n +CORES * i3] += temp_eq_sum + temp_visc_sum;
  dNd2pTdphidy_diff_all[n + CORES * i3] += temp_eq_sum + temp_diff_sum;
  dNd2pTdphidy_tot_all[n + CORES * i3] += (temp_eq_sum + temp_visc_sum + temp_diff_sum);
  
//   if (n==0)
//   std::cout<<" check: "<< n <<" "<< i3 << " "<< CORES<<" "<< temp_eq_sum<<" "<< dNd2pTdphidy_eq_all[n + CORES * i3]<< std::endl;
    
  // if (n == 0 && i3<10000 && i3 > 9990){
  //     std::cout<<n<<" "<< i3 <<" "<< temp_eq_sum <<std::endl;
  // }

  // dilepton polarization lambda_theta
  double p_lab_Min_norm[4];
  double p_lab_min_vec_square = sqrt( p_lab_Min[1]*p_lab_Min[1]
                                      + p_lab_Min[2]*p_lab_Min[2]
                                      + p_lab_Min[3]*p_lab_Min[3]);

  p_lab_Min_norm[1]=p_lab_Min[1]/p_lab_min_vec_square;
  p_lab_Min_norm[2]=p_lab_Min[2]/p_lab_min_vec_square;
  p_lab_Min_norm[3]=p_lab_Min[3]/p_lab_min_vec_square;

  double u_dot_kvec_norm = flow_u_mu_Min[1]*p_lab_Min_norm[1]
                          +flow_u_mu_Min[2]*p_lab_Min_norm[2]
                          +flow_u_mu_Min[3]*p_lab_Min_norm[3];

  double u_dot_kvec = flow_u_mu_Min[1]*p_lab_Min[1]
                      +flow_u_mu_Min[2]*p_lab_Min[2]
                      +flow_u_mu_Min[3]*p_lab_Min[3];

  double factor_pol_0 = pow((p_lab_Min[0]*u_dot_kvec_norm -
  p_lab_min_vec_square*flow_u_mu_Min[0]),2)/
                        (pow((p_lab_Min[0]*flow_u_mu_Min[0]-u_dot_kvec),2)-M_ll*M_ll);
  factor_pol_0 = factor_pol_0 - 1./3.0;


  ///////////// CS_frame
  double v_dilepton_boost[3] = {p_lab_Min[1]/p_lab_Min[0],p_lab_Min[2]/p_lab_Min[0],p_lab_Min[3]/p_lab_Min[0]};
  int rows = 4, cols = 4;
  double** lambda_munu = createA2DMatrix(rows, cols, 0.);  
  boost_matrix(lambda_munu, v_dilepton_boost[0], 
        v_dilepton_boost[1], v_dilepton_boost[2]);  

  double pA_beam_lab[4] = {1,0,0,1};
  double pB_beam_lab[4] = {1,0,0,-1};
  double pA_beam_lrf[4] = {0,0,0,0};
  double pB_beam_lrf[4] = {0,0,0,0};



    double uflow_dilepton_lrf[4] = {0,0,0,0};
     
    for (int boostj = 0; boostj < 4; boostj++) {
        //double checkk = 0.0;   
        for (int boosti = 0; boosti < 4; boosti++) {
            pA_beam_lrf[boostj] += lambda_munu[boostj][boosti]*pA_beam_lab[boosti];
            pB_beam_lrf[boostj] += lambda_munu[boostj][boosti]*pB_beam_lab[boosti];
            uflow_dilepton_lrf[boostj] += lambda_munu[boostj][boosti]*flow_u_mu_Min[boosti];
        }
        
    }

    
   

    double unit_z[4] = {pA_beam_lrf[0]-pB_beam_lrf[0],pA_beam_lrf[1]-pB_beam_lrf[1],pA_beam_lrf[2]-pB_beam_lrf[2],pA_beam_lrf[3]-pB_beam_lrf[3]};
    
    double unit_z_sq = sqrt(unit_z[1]* unit_z[1] + unit_z[2]* unit_z[2] + unit_z[3]* unit_z[3]);

    unit_z[1] =  unit_z[1]/unit_z_sq;
    unit_z[2] =  unit_z[2]/unit_z_sq;
    unit_z[3] =  unit_z[3]/unit_z_sq;

   

    double uz_sq_lrf = pow(uflow_dilepton_lrf[1]*unit_z[1] + uflow_dilepton_lrf[2]*unit_z[2] +uflow_dilepton_lrf[3]*unit_z[3],2);
    double u_vec_sq = uflow_dilepton_lrf[1]*uflow_dilepton_lrf[1] + uflow_dilepton_lrf[2]*uflow_dilepton_lrf[2] +uflow_dilepton_lrf[3]*uflow_dilepton_lrf[3];
    
     double check_u1 = flow_u_mu_Min[1]+(p_lab_Min[0]/M_ll -1)*u_dot_kvec_norm*p_lab_Min_norm[1] - p_lab_Min_norm[1]*flow_u_mu_Min[0]/M_ll;

    
    if (CS_frame == 1 )
    {
    double factor_pol_0_tem = uz_sq_lrf/u_vec_sq;            
    factor_pol_0 = factor_pol_0_tem - 1./3.0;
    }

    

    ///////////// CS_frame




  double factor_pol_1 = me*me/(M_ll*M_ll);
  double rho_delta= temp_eqT_sum - temp_eqL_sum;
  double rho_V = temp_eqT_sum*2 + temp_eqL_sum;

   


  double lambda_theta1 = factor_pol_0*(1.0 - 4.0*factor_pol_1)*rho_delta;
  double lambda_theta2 = ( 4.0*(1+2.0*factor_pol_1)*rho_V/3.0 - lambda_theta1
  );

  double lambda_theta = 3.0*lambda_theta1/lambda_theta2;

  if ( fabs(lambda_theta2) < 1e-18){
      lambda_theta = 0.0;
  }

  double dNd2pTdphidy_cell_lambda_theta =
  lambda_theta*temp_eq_sum/(1.0+lambda_theta/3.);
  double dNd2pTdphidy_cell_lambda_norm = temp_eq_sum/(1.0+lambda_theta/3.);

  

  if(std::isnan(dNd2pTdphidy_cell_lambda_theta))
  {
      std::cout<< "ERROR: lambda_theta is NAN !!!!!! "<<std::endl;
      std::cout<< temp_eqT_sum << " "<< temp_eqL_sum <<" "<<lambda_theta1<<" "<<factor_pol_0<<" "<<lambda_theta2<< " "<<dNd2pTdphidy_cell_lambda_theta<<std::endl;
  }

  dNd2pTdphidy_pol_lambda_theta_all[n + CORES * i3] += dNd2pTdphidy_cell_lambda_theta;
  dNd2pTdphidy_pol_lambda_norm_all[n + CORES * i3] += dNd2pTdphidy_cell_lambda_norm;


  // lambda_phi

  double p_beam[4] = {0,0,0,1};
  double phat_dot_khat = p_beam[3]*p_lab_Min_norm[3];
  double uflow_dot_khat = u_dot_kvec_norm;
  double phat_dot_uflow = p_beam[3]*flow_u_mu_Min[3];

  double phat_cross_khat_dot_uflow = p_lab_Min_norm[1]*flow_u_mu_Min[2]
                         - p_lab_Min_norm[2]*flow_u_mu_Min[1];

  double lambda_phi_ux_sq =
  pow(phat_dot_khat,2)*pow(uflow_dot_khat,2)+pow(phat_dot_uflow,2)
                          - 2*phat_dot_uflow*phat_dot_khat*uflow_dot_khat;

  lambda_phi_ux_sq = lambda_phi_ux_sq/(1- phat_dot_khat* phat_dot_khat);

  double lambda_phi_uy_sq =
  phat_cross_khat_dot_uflow*phat_cross_khat_dot_uflow/(1- phat_dot_khat*
  phat_dot_khat);

  double uflow_lrf_sq =
  pow((p_lab_Min[0]*flow_u_mu_Min[0]-u_dot_kvec)/M_ll,2) - 1.0;

  double uxsq_m_uysq_usq = (lambda_phi_ux_sq-lambda_phi_uy_sq)/uflow_lrf_sq;

  
  //cs_frame


  double unit_y[4] = {0.0, p_beam[2]*p_lab_Min_norm[3] - p_beam[3]*p_lab_Min_norm[2] , p_beam[3]*p_lab_Min_norm[1] - p_beam[1]*p_lab_Min_norm[3] , p_beam[1]*p_lab_Min_norm[2] - p_beam[2]*p_lab_Min_norm[1] };
  
  double unit_y_sq = sqrt(unit_y[1]* unit_y[1] + unit_y[2]* unit_y[2] + unit_y[3]* unit_y[3]);
  
  
  unit_y[1] =  unit_y[1]/unit_y_sq;
  unit_y[2] =  unit_y[2]/unit_y_sq;
  unit_y[3] =  unit_y[3]/unit_y_sq;


  double unit_x[4] = { 0.0, unit_y[2]*unit_z[3] - unit_y[3]*unit_z[2] , unit_y[3]*unit_z[1] - unit_y[1]*unit_z[3] , unit_y[1]*unit_z[2] - unit_y[2]*unit_z[1] };

  

  double unit_x_sq = sqrt(unit_x[1]* unit_x[1] + unit_x[2]* unit_x[2] + unit_x[3]* unit_x[3]);
  unit_x[1] =  unit_x[1]/unit_x_sq;
  unit_x[2] =  unit_x[2]/unit_x_sq;
  unit_x[3] =  unit_x[3]/unit_x_sq;



  double ux_sq_lrf = pow(uflow_dilepton_lrf[1]*unit_x[1] + uflow_dilepton_lrf[2]*unit_x[2] +uflow_dilepton_lrf[3]*unit_x[3],2);

  double uy_sq_lrf = pow(uflow_dilepton_lrf[1]*unit_y[1] + uflow_dilepton_lrf[2]*unit_y[2] +uflow_dilepton_lrf[3]*unit_y[3],2);
  
  if (CS_frame ==1){
      uxsq_m_uysq_usq = (ux_sq_lrf - uy_sq_lrf)/u_vec_sq;
  }

  //cs frame end






  double lambda_phi = (1.0 - 4.0*factor_pol_1)*rho_delta/lambda_theta2;

  if ( fabs(lambda_theta2) < 1e-18){
      lambda_phi = 0.0;
  }

  lambda_phi = uxsq_m_uysq_usq*lambda_phi;
  //lambda_phi = lambda_phi;

  double dNd2pTdphidy_cell_lambda_phi =
  lambda_phi*temp_eq_sum/(1.0+lambda_theta/3.);

  if(std::isnan(dNd2pTdphidy_cell_lambda_phi))
  {
      std::cout<< "ERROR: lambda_theta is NAN !!!!!! "<<std::endl;
      std::cout<< lambda_phi << " "<< lambda_theta2 <<" "<<temp_eq_sum<<" "<<factor_pol_1<<" "<<rho_V <<"  "<<lambda_theta1<<std::endl;
  }

  dNd2pTdphidy_pol_lambda_phi_all[n + CORES * i3] += dNd2pTdphidy_cell_lambda_phi;


  deleteA2DMatrix(lambda_munu, rows);


}

void ThermalPhoton::reduce_multile_core(){
  //std::cout << CORES <<" "<< nrapidity <<" "<< nphi <<" "<< np << " "<<nm<<std::endl;
  for (int k = 0; k < nrapidity; k++) {
      for (int m = 0; m < nphi; m++) {
          for (int l = 0; l < np; l++) {
              for (int j = 0; j < nm; j++) {

                  int i3 = (j+(l+(m+nphi * k) * np) * nm);

                  double dN_pTdpTdphidy_eq_tmp = 0.0; // reduction variable
                  double dN_pTdpTdphidy_eqT_tmp = 0.0;
                  double dN_pTdpTdphidy_eqL_tmp = 0.0;
                  double dN_pTdpTdphidy_visc_tmp = 0.0;
                  double dN_pTdpTdphidy_diff_tmp = 0.0;
                  double dN_pTdpTdphidy_tot_tmp = 0.0;
                  double dN_pTdpTdphidy_pol_lambda_theta_tmp = 0.0;
                  double dN_pTdpTdphidy_pol_lambda_norm_tmp = 0.0;
                  double dN_pTdpTdphidy_pol_lambda_phi_tmp = 0.0;

                  for(long n = 0; n < CORES; n++)
                  {   
                
                      dN_pTdpTdphidy_eq_tmp += dNd2pTdphidy_eq_all[n+CORES*i3]; 
                      dN_pTdpTdphidy_eqT_tmp += dNd2pTdphidy_eqT_all[n+CORES*i3];
                      dN_pTdpTdphidy_eqL_tmp += dNd2pTdphidy_eqL_all[n+CORES*i3];
                      dN_pTdpTdphidy_visc_tmp += dNd2pTdphidy_visc_all[n+CORES*i3];
                      dN_pTdpTdphidy_diff_tmp += dNd2pTdphidy_diff_all[n+CORES*i3];
                      dN_pTdpTdphidy_tot_tmp += dNd2pTdphidy_tot_all[n+CORES*i3];
                      dN_pTdpTdphidy_pol_lambda_theta_tmp += dNd2pTdphidy_pol_lambda_theta_all[n+CORES*i3];
                      dN_pTdpTdphidy_pol_lambda_norm_tmp += dNd2pTdphidy_pol_lambda_norm_all[n+CORES*i3];
                      dN_pTdpTdphidy_pol_lambda_phi_tmp += dNd2pTdphidy_pol_lambda_phi_all[n+CORES*i3];

                  } // sum over the cores

                  dNd2pTdphidy_eq[j][l][m][k] = dN_pTdpTdphidy_eq_tmp;
                  dNd2pTdphidy_eqT[j][l][m][k] = dN_pTdpTdphidy_eqT_tmp;
                  dNd2pTdphidy_eqL[j][l][m][k] = dN_pTdpTdphidy_eqL_tmp;
                  dNd2pTdphidy_visc[j][l][m][k] = dN_pTdpTdphidy_visc_tmp;
                  dNd2pTdphidy_diff[j][l][m][k] = dN_pTdpTdphidy_diff_tmp;
                  dNd2pTdphidy_tot[j][l][m][k] = dN_pTdpTdphidy_tot_tmp;

                  dNd2pTdphidy_pol_lambda_theta[j][l][m][k] = dN_pTdpTdphidy_pol_lambda_theta_tmp;
                  dNd2pTdphidy_pol_lambda_norm[j][l][m][k] = dN_pTdpTdphidy_pol_lambda_norm_tmp;
                  dNd2pTdphidy_pol_lambda_phi[j][l][m][k] = dN_pTdpTdphidy_pol_lambda_phi_tmp;

                  // distribution in (T, tau)
                //   if (differential_flag == 1) {
                //       for (int iT = 0; iT < nTcut; iT++) {
                //            for (int it = 0; it < n_tau_cut; it++) {

                //               double dN_pTdpTdphidydTdtau_eq_tmp = 0.0;
                //               double dN_pTdpTdphidydTdtau_visc_tmp = 0.0;
                //               double dN_pTdpTdphidydTdtau_diff_tmp = 0.0;
                //               double dN_pTdpTdphidydTdtau_tot_tmp = 0.0;

                //               #pragma omp simd reduction(+:dN_pTdpTdphidydTdtau_eq_tmp,dN_pTdpTdphidydTdtau_visc_tmp,\
                //                     dN_pTdpTdphidydTdtau_diff_tmp,dN_pTdpTdphidydTdtau_tot_tmp)
                //               for(long n = 0; n < CORES; n++) {
                //                   dN_pTdpTdphidydTdtau_eq_tmp +=
                //                   dNd2pTdphidydTdtau_eq_all[iT][it][n+CORES*i3];
                //                   dN_pTdpTdphidydTdtau_visc_tmp +=
                //                   dNd2pTdphidydTdtau_visc_all[iT][it][n+CORES*i3];
                //                   dN_pTdpTdphidydTdtau_diff_tmp +=
                //                   dNd2pTdphidydTdtau_diff_all[iT][it][n+CORES*i3];
                //                   dN_pTdpTdphidydTdtau_tot_tmp +=
                //                   dNd2pTdphidydTdtau_tot_all[iT][it][n+CORES*i3];
                //               }

                //               dNd2pTdphidydTdtau_eq[iT][it][j][l][m][k] =
                //               dN_pTdpTdphidydTdtau_eq_tmp;
                //               dNd2pTdphidydTdtau_visc[iT][it][j][l][m][k] =
                //               dN_pTdpTdphidydTdtau_visc_tmp;
                //               dNd2pTdphidydTdtau_diff[iT][it][j][l][m][k] =
                //               dN_pTdpTdphidydTdtau_diff_tmp;
                //               dNd2pTdphidydTdtau_tot[iT][it][j][l][m][k] =
                //               dN_pTdpTdphidydTdtau_tot_tmp;

                //           }
                //       }
                //   }

              } // M_ll
          } // p_T
      } // phi_p
  } // y

}


//  void ThermalPhoton::calThermalPhotonemission_3d() {

// }

// functions to get distributions in T and tau; contributions from all channels
// are included
void ThermalPhoton::calPhoton_SpMatrix_dTdtau(
    double ******dNd2pTdphidydTdtau_eq_temp,
    double ******dNd2pTdphidydTdtau_visc_temp,
    double ******dNd2pTdphidydTdtau_diff_temp,
    double ******dNd2pTdphidydTdtau_tot_temp) {
  for (int i = 0; i < nTcut; i++) {
    for (int j = 0; j < n_tau_cut; j++) {
      for (int k = 0; k < nm; k++) {
        for (int l = 0; l < np; l++) {
          for (int m = 0; m < nphi; m++) {
            for (int n = 0; n < nrapidity; n++) {
              dNd2pTdphidydTdtau_eq[i][j][k][l][m][n] =
                  dNd2pTdphidydTdtau_eq_temp[i][j][k][l][m][n];
              dNd2pTdphidydTdtau_visc[i][j][k][l][m][n] =
                  dNd2pTdphidydTdtau_visc_temp[i][j][k][l][m][n];
              dNd2pTdphidydTdtau_diff[i][j][k][l][m][n] =
                  dNd2pTdphidydTdtau_diff_temp[i][j][k][l][m][n];
              dNd2pTdphidydTdtau_tot[i][j][k][l][m][n] =
                  dNd2pTdphidydTdtau_tot_temp[i][j][k][l][m][n];
            }
          }
        }
      }
    }
  }
}



void ThermalPhoton::outputPhoton_Spectra_full_diff(string path,
                                                   double Tcut_high,
                                                   double Tcut_low,
                                                   double tau_cut_high,
                                                   double tau_cut_low) {

  double dT = (Tcut_high - Tcut_low) / (nTcut - 1);
  double dtau = (tau_cut_high - tau_cut_low) / (n_tau_cut - 1);

  ostringstream filename_Sp_full_diff_eq;
  ostringstream filename_Sp_full_diff_visc;
  ostringstream filename_Sp_full_diff_diff;
  ostringstream filename_Sp_full_diff_tot;

  filename_Sp_full_diff_eq << path << emissionProcess_name << "_Sp_full_eq.dat";
  filename_Sp_full_diff_visc << path << emissionProcess_name
                             << "_Sp_full_visc.dat";
  filename_Sp_full_diff_diff << path << emissionProcess_name
                             << "_Sp_full_diff.dat";
  filename_Sp_full_diff_tot << path << emissionProcess_name
                            << "_Sp_full_tot.dat";

  ofstream ofeq(filename_Sp_full_diff_eq.str().c_str());
  ofstream ofvisc(filename_Sp_full_diff_visc.str().c_str());
  ofstream ofdiff(filename_Sp_full_diff_diff.str().c_str());
  ofstream oftot(filename_Sp_full_diff_tot.str().c_str());

  for (int i = 0; i < nTcut; i++) {
    double T_local = Tcut_low + i * dT;
    for (int j = 0; j < n_tau_cut; j++) {
      double tau_local = tau_cut_low + j * dtau;
      for (int k = 0; k < nm; k++) {
        double Mll_local = getDileptonMass(k);
        for (int l = 0; l < np; l++) {
          double pT_local = getPhotonp(l);
          for (int m = 0; m < nphi; m++) {
            double phi_local = getPhotonphi(m);
            for (int n = 0; n < nrapidity; n++) {
              double y_local = getPhotonrapidity(n);

              ofeq << scientific << setw(18) << setprecision(8) << T_local
                   << " " << tau_local << " " << Mll_local << " " << pT_local
                   << " " << phi_local << " " << y_local << " "
                   << dNd2pTdphidydTdtau_eq[i][j][k][l][m][n] << std::endl;
              ofvisc << scientific << setw(18) << setprecision(8) << T_local
                     << " " << tau_local << " " << Mll_local << " " << pT_local
                     << " " << phi_local << " " << y_local << " "
                     << dNd2pTdphidydTdtau_visc[i][j][k][l][m][n] << std::endl;
              ofdiff << scientific << setw(18) << setprecision(8) << T_local
                     << " " << tau_local << " " << Mll_local << " " << pT_local
                     << " " << phi_local << " " << y_local << " "
                     << dNd2pTdphidydTdtau_diff[i][j][k][l][m][n] << std::endl;

              oftot << scientific << setw(18) << setprecision(8) << T_local
                    << " " << tau_local << " " << Mll_local << " " << pT_local
                    << " " << phi_local << " " << y_local << " "
                    << dNd2pTdphidydTdtau_tot[i][j][k][l][m][n] << std::endl;
            }
          }
        }
      }
    }
  }

  ofeq.close();
  ofvisc.close();
  ofdiff.close();
  oftot.close();
}

void ThermalPhoton::calPhoton_Spectra_dTdtau() {
  // calculate the photon spectra at T-tau interval
  // integrated out phi and rapidity

  for (int i = 0; i < nTcut; i++) {
    for (int j = 0; j < n_tau_cut; j++) {

      for (int m = 0; m < nm; m++) {
        for (int k = 0; k < np; k++) {
          for (int l = 0; l < nphi; l++) {
            for (int irap = 0; irap < nrapidity; irap++) {
              dNpTdpTdydTdtau_eq[i][j][m][k] +=
                  (dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap] * phi_weight[l] *
                   y_weight[irap] * dy / Dy);
              dNpTdpTdydTdtau_visc[i][j][m][k] +=
                  (dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap] *
                   phi_weight[l] * y_weight[irap] * dy / Dy);
              dNpTdpTdydTdtau_diff[i][j][m][k] +=
                  (dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap] *
                   phi_weight[l] * y_weight[irap] * dy / Dy);
              dNpTdpTdydTdtau_tot[i][j][m][k] +=
                  (dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap] * phi_weight[l] *
                   y_weight[irap] * dy / Dy);
            }
          }
          //
          dNpTdpTdydTdtau_eq[i][j][m][k] =
              dNpTdpTdydTdtau_eq[i][j][m][k] / (2 * M_PI);
          dNpTdpTdydTdtau_visc[i][j][m][k] =
              dNpTdpTdydTdtau_visc[i][j][m][k] / (2 * M_PI);
          dNpTdpTdydTdtau_diff[i][j][m][k] =
              dNpTdpTdydTdtau_diff[i][j][m][k] / (2 * M_PI);
          dNpTdpTdydTdtau_tot[i][j][m][k] =
              dNpTdpTdydTdtau_tot[i][j][m][k] / (2 * M_PI);
        }
      }
    }
  }

  for (int i = 0; i < nTcut; i++) {
    for (int j = 0; j < n_tau_cut; j++) {
      for (int m = 0; m < nm; m++) {
        for (int k = 0; k < np; k++) {

          dNpTdpTdydT_eq[i][m][k] += dNpTdpTdydTdtau_eq[i][j][m][k];
          dNpTdpTdydT_visc[i][m][k] += dNpTdpTdydTdtau_visc[i][j][m][k];
          dNpTdpTdydT_diff[i][m][k] += dNpTdpTdydTdtau_diff[i][j][m][k];
          dNpTdpTdydT_tot[i][m][k] += dNpTdpTdydTdtau_tot[i][j][m][k];

          dNpTdpTdydtau_eq[j][m][k] += dNpTdpTdydTdtau_eq[i][j][m][k];
          dNpTdpTdydtau_visc[j][m][k] += dNpTdpTdydTdtau_visc[i][j][m][k];
          dNpTdpTdydtau_diff[j][m][k] += dNpTdpTdydTdtau_diff[i][j][m][k];
          dNpTdpTdydtau_tot[j][m][k] += dNpTdpTdydTdtau_tot[i][j][m][k];
        }
      }
    }
  }
}

void ThermalPhoton::outputPhoton_Spectra_dTdtau(string path, double Tcut_high,
                                                double Tcut_low,
                                                double tau_cut_high,
                                                double tau_cut_low) {
  // calculate the inverse slope of the photon spectra at T-tau interval
  // integrated out phi and rapidity
  // pT differential spectra, dN/(2pi pTdpT MdM dy)

  double dT = (Tcut_high - Tcut_low) / (nTcut - 1);
  double dtau = (tau_cut_high - tau_cut_low) / (n_tau_cut - 1);

  ostringstream filename_SpdTdtau_eq;
  ostringstream filename_SpdTdtau_visc;
  ostringstream filename_SpdTdtau_diff;
  ostringstream filename_SpdTdtau_tot;

  filename_SpdTdtau_eq << path << emissionProcess_name << "_SpdTdtau_eq.dat";
  filename_SpdTdtau_visc << path << emissionProcess_name
                         << "_SpdTdtau_visc.dat";
  filename_SpdTdtau_diff << path << emissionProcess_name
                         << "_SpdTdtau_diff.dat";
  filename_SpdTdtau_tot << path << emissionProcess_name << "_SpdTdtau_tot.dat";

  ofstream ofeq(filename_SpdTdtau_eq.str().c_str());
  ofstream ofvisc(filename_SpdTdtau_visc.str().c_str());
  ofstream ofdiff(filename_SpdTdtau_diff.str().c_str());
  ofstream oftot(filename_SpdTdtau_tot.str().c_str());

  for (int i = 0; i < nTcut; i++) {
    double T_local = Tcut_low + i * dT;
    for (int j = 0; j < n_tau_cut; j++) {
      double tau_local = tau_cut_low + j * dtau;

      for (int m = 0; m < nm; m++) {
        ofeq << scientific << setw(18) << setprecision(8) << T_local << "   "
             << tau_local << "   ";
        ofdiff << scientific << setw(18) << setprecision(8) << T_local << "   "
               << tau_local << "   ";
        oftot << scientific << setw(18) << setprecision(8) << T_local << "   "
              << tau_local << "   ";
        for (int k = 0; k < np; k++) {
          ofeq << dNpTdpTdydTdtau_eq[i][j][m][k] << "   ";
          ofvisc << dNpTdpTdydTdtau_visc[i][j][m][k] << "   ";
          ofdiff << dNpTdpTdydTdtau_diff[i][j][m][k] << "   ";
          oftot << dNpTdpTdydTdtau_tot[i][j][m][k] << "   ";
        }
        ofeq << endl;
        ofvisc << endl;
        ofdiff << endl;
        oftot << endl;
      }
    }
  }

  ofeq.close();
  ofvisc.close();
  ofdiff.close();
  oftot.close();

  // in temperature
  ostringstream filename_SpdT_eq;
  ostringstream filename_SpdT_visc;
  ostringstream filename_SpdT_diff;
  ostringstream filename_SpdT_tot;
  filename_SpdT_eq << path << emissionProcess_name << "_SpdT_eq.dat";
  filename_SpdT_visc << path << emissionProcess_name << "_SpdT_visc.dat";
  filename_SpdT_diff << path << emissionProcess_name << "_SpdT_diff.dat";
  filename_SpdT_tot << path << emissionProcess_name << "_SpdT_tot.dat";

  ofstream ofdTeq(filename_SpdT_eq.str().c_str());
  ofstream ofdTvisc(filename_SpdT_visc.str().c_str());
  ofstream ofdTdiff(filename_SpdT_diff.str().c_str());
  ofstream ofdTtot(filename_SpdT_tot.str().c_str());

  for (int i = 0; i < nTcut; i++) {
    double T_local = Tcut_low + i * dT;

    for (int m = 0; m < nm; m++) {
      ofdTeq << scientific << setw(18) << setprecision(8) << T_local << "   ";
      ofdTvisc << scientific << setw(18) << setprecision(8) << T_local << "   ";
      ofdTdiff << scientific << setw(18) << setprecision(8) << T_local << "   ";
      ofdTtot << scientific << setw(18) << setprecision(8) << T_local << "   ";
      for (int k = 0; k < np; k++) {
        ofdTeq << dNpTdpTdydT_eq[i][m][k] << "   ";
        ofdTvisc << dNpTdpTdydT_visc[i][m][k] << "   ";
        ofdTdiff << dNpTdpTdydT_diff[i][m][k] << "   ";
        ofdTtot << dNpTdpTdydT_tot[i][m][k] << "   ";
      }
      ofdTeq << endl;
      ofdTvisc << endl;
      ofdTdiff << endl;
      ofdTtot << endl;
    }
  }

  ofdTeq.close();
  ofdTvisc.close();
  ofdTdiff.close();
  ofdTtot.close();

  // in proper time
  ostringstream filename_Spdtau_eq;
  ostringstream filename_Spdtau_visc;
  ostringstream filename_Spdtau_diff;
  ostringstream filename_Spdtau_tot;
  filename_Spdtau_eq << path << emissionProcess_name << "_Spdtau_eq.dat";
  filename_Spdtau_visc << path << emissionProcess_name << "_Spdtau_visc.dat";
  filename_Spdtau_diff << path << emissionProcess_name << "_Spdtau_diff.dat";
  filename_Spdtau_tot << path << emissionProcess_name << "_Spdtau_tot.dat";

  ofstream ofdtaueq(filename_Spdtau_eq.str().c_str());
  ofstream ofdtauvisc(filename_Spdtau_visc.str().c_str());
  ofstream ofdtaudiff(filename_Spdtau_diff.str().c_str());
  ofstream ofdtautot(filename_Spdtau_tot.str().c_str());

  for (int j = 0; j < n_tau_cut; j++) {
    double tau_local = tau_cut_low + j * dtau;

    for (int m = 0; m < nm; m++) {
      ofdtaueq << scientific << setw(18) << setprecision(8) << tau_local
               << "   ";
      ofdtauvisc << scientific << setw(18) << setprecision(8) << tau_local
                 << "   ";
      ofdtaudiff << scientific << setw(18) << setprecision(8) << tau_local
                 << "   ";
      ofdtautot << scientific << setw(18) << setprecision(8) << tau_local
                << "   ";
      for (int k = 0; k < np; k++) {
        ofdtaueq << dNpTdpTdydtau_eq[j][m][k] << "   ";
        ofdtauvisc << dNpTdpTdydtau_visc[j][m][k] << "   ";
        ofdtaudiff << dNpTdpTdydtau_diff[j][m][k] << "   ";
        ofdtautot << dNpTdpTdydtau_tot[j][m][k] << "   ";
      }
      ofdtaueq << endl;
      ofdtauvisc << endl;
      ofdtaudiff << endl;
      ofdtautot << endl;
    }
  }

  ofdtaueq.close();
  ofdtauvisc.close();
  ofdtaudiff.close();
  ofdtautot.close();
}

void ThermalPhoton::calPhoton_Spvn_dTdtau() {
  // calculate the dilepton, dN/dydM yields and vn
  // integrated out pT, phi and rapidity
  double eps = 1e-15;
  for (int i = 0; i < nTcut; i++) {
    for (int j = 0; j < n_tau_cut; j++) {

      for (int m = 0; m < nm; m++) {
        for (int k = 0; k < np; k++) {
          for (int l = 0; l < nphi; l++) {
            double weight = p[k] * p_weight[k] *
                            phi_weight[l]; // pT and phi_p integrated out
            for (int irap = 0; irap < nrapidity; irap++) {

              weight *= y_weight[irap] * dy / Dy;

              dNdydTdtau_eq[i][j][m] +=
                  (dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap] * weight);
              dNdydTdtau_visc[i][j][m] +=
                  (dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap] * weight);
              dNdydTdtau_diff[i][j][m] +=
                  (dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap] * weight);
              dNdydTdtau_tot[i][j][m] +=
                  (dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap] * weight);

              for (int order = 0; order < norder; order++) {
                vndTdtau_cos_eq[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap] * weight *
                     cos(order * phi[l]));
                vndTdtau_sin_eq[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap] * weight *
                     sin(order * phi[l]));
                vndTdtau_cos_visc[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap] * weight *
                     cos(order * phi[l]));
                vndTdtau_sin_visc[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap] * weight *
                     sin(order * phi[l]));
                vndTdtau_cos_diff[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap] * weight *
                     cos(order * phi[l]));
                vndTdtau_sin_diff[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap] * weight *
                     sin(order * phi[l]));
                vndTdtau_cos_tot[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap] * weight *
                     cos(order * phi[l]));
                vndTdtau_sin_tot[i][j][m][order] +=
                    (dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap] * weight *
                     sin(order * phi[l]));
              }
            }
          }
        }
      }

      for (int m = 0; m < nm; m++) {
        for (int order = 1; order < norder; order++) {
          vndTdtau_cos_eq[i][j][m][order] = (vndTdtau_cos_eq[i][j][m][order] /
                                             (dNdydTdtau_eq[i][j][m] + eps));
          vndTdtau_sin_eq[i][j][m][order] = (vndTdtau_sin_eq[i][j][m][order] /
                                             (dNdydTdtau_eq[i][j][m] + eps));
          vndTdtau_cos_visc[i][j][m][order] =
              (vndTdtau_cos_visc[i][j][m][order] /
               (dNdydTdtau_visc[i][j][m] + eps));
          vndTdtau_sin_visc[i][j][m][order] =
              (vndTdtau_sin_visc[i][j][m][order] /
               (dNdydTdtau_visc[i][j][m] + eps));
          vndTdtau_cos_diff[i][j][m][order] =
              (vndTdtau_cos_diff[i][j][m][order] /
               (dNdydTdtau_diff[i][j][m] + eps));
          vndTdtau_sin_diff[i][j][m][order] =
              (vndTdtau_sin_diff[i][j][m][order] /
               (dNdydTdtau_diff[i][j][m] + eps));
          vndTdtau_cos_tot[i][j][m][order] = (vndTdtau_cos_tot[i][j][m][order] /
                                              (dNdydTdtau_tot[i][j][m] + eps));
          vndTdtau_sin_tot[i][j][m][order] = (vndTdtau_sin_tot[i][j][m][order] /
                                              (dNdydTdtau_tot[i][j][m] + eps));
        }
      }
    }
  }

  // dNdydM on both sides
  for (int i = 0; i < nTcut; i++) {
    for (int j = 0; j < n_tau_cut; j++) {
      for (int m = 0; m < nm; m++) {

        dNdydT_eq[i][m] += dNdydTdtau_eq[i][j][m];
        dNdydT_visc[i][m] += dNdydTdtau_visc[i][j][m];
        dNdydT_diff[i][m] += dNdydTdtau_diff[i][j][m];
        dNdydT_tot[i][m] += dNdydTdtau_tot[i][j][m];

        dNdydtau_eq[j][m] += dNdydTdtau_eq[i][j][m];
        dNdydtau_visc[j][m] += dNdydTdtau_visc[i][j][m];
        dNdydtau_diff[j][m] += dNdydTdtau_diff[i][j][m];
        dNdydtau_tot[j][m] += dNdydTdtau_tot[i][j][m];
      }
    }
  }
}

void ThermalPhoton::outputPhoton_Spvn_dTdtau(string path, double Tcut_high,
                                             double Tcut_low,
                                             double tau_cut_high,
                                             double tau_cut_low) {

  double dT = (Tcut_high - Tcut_low) / (nTcut - 1);
  double dtau = (tau_cut_high - tau_cut_low) / (n_tau_cut - 1);

  // yields in (T, tau)
  ostringstream filename_stream_dNdydTdtau_eq;
  ostringstream filename_stream_dNdydTdtau_visc;
  ostringstream filename_stream_dNdydTdtau_diff;
  ostringstream filename_stream_dNdydTdtau_tot;

  filename_stream_dNdydTdtau_eq << path << emissionProcess_name
                                << "_dNdydTdtau_eq.dat";
  filename_stream_dNdydTdtau_visc << path << emissionProcess_name
                                  << "_dNdydTdtau_visc.dat";
  filename_stream_dNdydTdtau_diff << path << emissionProcess_name
                                  << "_dNdydTdtau_diff.dat";
  filename_stream_dNdydTdtau_tot << path << emissionProcess_name
                                 << "_dNdydTdtau_tot.dat";

  ofstream fphotondNdy_eq(filename_stream_dNdydTdtau_eq.str().c_str());
  ofstream fphotondNdy_visc(filename_stream_dNdydTdtau_visc.str().c_str());
  ofstream fphotondNdy_diff(filename_stream_dNdydTdtau_diff.str().c_str());
  ofstream fphotondNdy_tot(filename_stream_dNdydTdtau_tot.str().c_str());

  for (int i = 0; i < nTcut; i++) {
    double T_local = Tcut_low + i * dT;
    for (int j = 0; j < n_tau_cut; j++) {
      double tau_local = tau_cut_low + j * dtau;
      fphotondNdy_eq << scientific << setw(18) << setprecision(8) << T_local
                     << "   " << tau_local << "   ";
      fphotondNdy_visc << scientific << setw(18) << setprecision(8) << T_local
                       << "   " << tau_local << "   ";
      fphotondNdy_diff << scientific << setw(18) << setprecision(8) << T_local
                       << "   " << tau_local << "   ";
      fphotondNdy_tot << scientific << setw(18) << setprecision(8) << T_local
                      << "   " << tau_local << "   ";
      for (int m = 0; m < nm; m++) {
        double M_ll = M[m];
        fphotondNdy_eq << dNdydTdtau_eq[i][j][m] * M_ll << "    ";
        fphotondNdy_visc << dNdydTdtau_visc[i][j][m] * M_ll << "    ";
        fphotondNdy_diff << dNdydTdtau_diff[i][j][m] * M_ll << "    ";
        fphotondNdy_tot << dNdydTdtau_tot[i][j][m] * M_ll << "    ";
      }
      fphotondNdy_eq << endl;
      fphotondNdy_visc << endl;
      fphotondNdy_diff << endl;
      fphotondNdy_tot << endl;
    }
  }
  fphotondNdy_eq.close();
  fphotondNdy_visc.close();
  fphotondNdy_diff.close();
  fphotondNdy_tot.close();

  // yields in T
  ostringstream filename_stream_dNdydT_eq;
  ostringstream filename_stream_dNdydT_visc;
  ostringstream filename_stream_dNdydT_diff;
  ostringstream filename_stream_dNdydT_tot;

  filename_stream_dNdydT_eq << path << emissionProcess_name << "_dNdydT_eq.dat";
  filename_stream_dNdydT_visc << path << emissionProcess_name
                              << "_dNdydT_visc.dat";
  filename_stream_dNdydT_diff << path << emissionProcess_name
                              << "_dNdydT_diff.dat";
  filename_stream_dNdydT_tot << path << emissionProcess_name
                             << "_dNdydT_tot.dat";

  ofstream fphotondNdydT_eq(filename_stream_dNdydT_eq.str().c_str());
  ofstream fphotondNdydT_visc(filename_stream_dNdydT_visc.str().c_str());
  ofstream fphotondNdydT_diff(filename_stream_dNdydT_diff.str().c_str());
  ofstream fphotondNdydT_tot(filename_stream_dNdydT_tot.str().c_str());

  for (int i = 0; i < nTcut; i++) {
    double T_local = Tcut_low + i * dT;

    fphotondNdydT_eq << scientific << setw(18) << setprecision(8) << T_local
                     << "   ";
    fphotondNdydT_visc << scientific << setw(18) << setprecision(8) << T_local
                       << "   ";
    fphotondNdydT_diff << scientific << setw(18) << setprecision(8) << T_local
                       << "   ";
    fphotondNdydT_tot << scientific << setw(18) << setprecision(8) << T_local
                      << "   ";
    for (int m = 0; m < nm; m++) {
      double M_ll = M[m];
      fphotondNdydT_eq << dNdydT_eq[i][m] * M_ll << "    ";
      fphotondNdydT_visc << dNdydT_visc[i][m] * M_ll << "    ";
      fphotondNdydT_diff << dNdydT_diff[i][m] * M_ll << "    ";
      fphotondNdydT_tot << dNdydT_tot[i][m] * M_ll << "    ";
    }
    fphotondNdydT_eq << endl;
    fphotondNdydT_visc << endl;
    fphotondNdydT_diff << endl;
    fphotondNdydT_tot << endl;
  }

  fphotondNdydT_eq.close();
  fphotondNdydT_visc.close();
  fphotondNdydT_diff.close();
  fphotondNdydT_tot.close();

  // yields in tau
  ostringstream filename_stream_dNdydtau_eq;
  ostringstream filename_stream_dNdydtau_visc;
  ostringstream filename_stream_dNdydtau_diff;
  ostringstream filename_stream_dNdydtau_tot;

  filename_stream_dNdydtau_eq << path << emissionProcess_name
                              << "_dNdydtau_eq.dat";
  filename_stream_dNdydtau_visc << path << emissionProcess_name
                                << "_dNdydtau_visc.dat";
  filename_stream_dNdydtau_diff << path << emissionProcess_name
                                << "_dNdydtau_diff.dat";
  filename_stream_dNdydtau_tot << path << emissionProcess_name
                               << "_dNdydtau_tot.dat";

  ofstream fphotondNdydtau_eq(filename_stream_dNdydtau_eq.str().c_str());
  ofstream fphotondNdydtau_visc(filename_stream_dNdydtau_visc.str().c_str());
  ofstream fphotondNdydtau_diff(filename_stream_dNdydtau_diff.str().c_str());
  ofstream fphotondNdydtau_tot(filename_stream_dNdydtau_tot.str().c_str());

  for (int i = 0; i < n_tau_cut; i++) {
    double tau_local = tau_cut_low + i * dtau;

    fphotondNdydtau_eq << scientific << setw(18) << setprecision(8) << tau_local
                       << "   ";
    fphotondNdydtau_visc << scientific << setw(18) << setprecision(8)
                         << tau_local << "   ";
    fphotondNdydtau_diff << scientific << setw(18) << setprecision(8)
                         << tau_local << "   ";
    fphotondNdydtau_tot << scientific << setw(18) << setprecision(8)
                        << tau_local << "   ";
    for (int m = 0; m < nm; m++) {
      double M_ll = M[m];
      fphotondNdydtau_eq << dNdydtau_eq[i][m] * M_ll << "    ";
      fphotondNdydtau_visc << dNdydtau_visc[i][m] * M_ll << "    ";
      fphotondNdydtau_diff << dNdydtau_diff[i][m] * M_ll << "    ";
      fphotondNdydtau_tot << dNdydtau_tot[i][m] * M_ll << "    ";
    }
    fphotondNdydtau_eq << endl;
    fphotondNdydtau_visc << endl;
    fphotondNdydtau_diff << endl;
    fphotondNdydtau_tot << endl;
  }

  fphotondNdydtau_eq.close();
  fphotondNdydtau_visc.close();
  fphotondNdydtau_diff.close();
  fphotondNdydtau_tot.close();

  // flow coefficients

  ostringstream filename_stream_vndTdtau_eq;
  ostringstream filename_stream_vndTdtau_visc;
  ostringstream filename_stream_vndTdtau_diff;
  ostringstream filename_stream_vndTdtau_tot;

  filename_stream_vndTdtau_eq << path << emissionProcess_name
                              << "_vn_dTdtau_eq.dat";
  filename_stream_vndTdtau_visc << path << emissionProcess_name
                                << "_vn_dTdtau_visc.dat";
  filename_stream_vndTdtau_diff << path << emissionProcess_name
                                << "_vn_dTdtau_diff.dat";
  filename_stream_vndTdtau_tot << path << emissionProcess_name
                               << "_vn_dTdtau_tot.dat";

  ofstream fphotonvn_eq(filename_stream_vndTdtau_eq.str().c_str());
  ofstream fphotonvn_visc(filename_stream_vndTdtau_visc.str().c_str());
  ofstream fphotonvn_diff(filename_stream_vndTdtau_diff.str().c_str());
  ofstream fphotonvn_tot(filename_stream_vndTdtau_tot.str().c_str());

  for (int i = 0; i < nTcut; i++) {
    double T_local = Tcut_low + i * dT;
    for (int j = 0; j < n_tau_cut; j++) {
      double tau_local = tau_cut_low + j * dtau;
      fphotonvn_eq << scientific << setw(18) << setprecision(8) << T_local
                   << "   " << tau_local << "   ";
      fphotonvn_visc << scientific << setw(18) << setprecision(8) << T_local
                     << "   " << tau_local << "   ";
      fphotonvn_diff << scientific << setw(18) << setprecision(8) << T_local
                     << "   " << tau_local << "   ";
      fphotonvn_tot << scientific << setw(18) << setprecision(8) << T_local
                    << "   " << tau_local << "   ";

      for (int m = 0; m < nm; m++) {
        for (int order = 1; order < norder; order++) {
          fphotonvn_eq << order << "   " << vndTdtau_cos_eq[i][j][m][order]
                       << "    " << vndTdtau_sin_eq[i][j][m][order] << "    "
                       << sqrt(pow(vndTdtau_cos_eq[i][j][m][order], 2) +
                               pow(vndTdtau_sin_eq[i][j][m][order], 2))
                       << "  ";
          fphotonvn_visc << order << "   " << vndTdtau_cos_visc[i][j][m][order]
                         << "    " << vndTdtau_sin_visc[i][j][m][order]
                         << "    "
                         << sqrt(pow(vndTdtau_cos_visc[i][j][m][order], 2) +
                                 pow(vndTdtau_sin_visc[i][j][m][order], 2))
                         << "  ";
          fphotonvn_diff << order << "   " << vndTdtau_cos_diff[i][j][m][order]
                         << "    " << vndTdtau_sin_diff[i][j][m][order]
                         << "    "
                         << sqrt(pow(vndTdtau_cos_diff[i][j][m][order], 2) +
                                 pow(vndTdtau_sin_diff[i][j][m][order], 2))
                         << "  ";
          fphotonvn_tot << order << "   " << vndTdtau_cos_tot[i][j][m][order]
                        << "    " << vndTdtau_sin_tot[i][j][m][order] << "    "
                        << sqrt(pow(vndTdtau_cos_tot[i][j][m][order], 2) +
                                pow(vndTdtau_sin_tot[i][j][m][order], 2))
                        << "  ";
        }
      }
      fphotonvn_eq << endl;
      fphotonvn_visc << endl;
      fphotonvn_diff << endl;
      fphotonvn_tot << endl;
    }
  }
  fphotonvn_eq.close();
  fphotonvn_visc.close();
  fphotonvn_diff.close();
  fphotonvn_tot.close();
}


void ThermalPhoton::calPhoton_SpvnpT_shell() {

    calPhoton_SpvnpT(dNd2pTdphidy_eq, vnpT_cos_eq,vnpT_sin_eq,
                    vn_cos_eq, vn_sin_eq,
                    dNd2Mdy_eq, dNd2pTd2M_eq,
                    dNd2pTd2Mdy_eq,vnMpTy_cos_eq,vnMpTy_sin_eq);
    
    
    calPhoton_SpvnpT(dNd2pTdphidy_eqT,  vnpT_cos_eqT, vnpT_sin_eqT,
                    vn_cos_eqT,  vn_sin_eqT,
                    dNd2Mdy_eqT,  dNd2pTd2M_eqT,
                    dNd2pTd2Mdy_eqT, vnMpTy_cos_eqT, vnMpTy_sin_eqT);

    calPhoton_SpvnpT( dNd2pTdphidy_eqL, vnpT_cos_eqL, vnpT_sin_eqL,
                         vn_cos_eqL,  vn_sin_eqL,
                         dNd2Mdy_eqL,dNd2pTd2M_eqL,
                         dNd2pTd2Mdy_eqL, vnMpTy_cos_eqL, vnMpTy_sin_eqL);
                        
    calPhoton_SpvnpT(dNd2pTdphidy_visc, vnpT_cos_visc, vnpT_sin_visc,
                         vn_cos_visc,  vn_sin_visc,
                          dNd2Mdy_visc,  dNd2pTd2M_visc,
                         dNd2pTd2Mdy_visc, vnMpTy_cos_visc, vnMpTy_sin_visc);
    
    calPhoton_SpvnpT(dNd2pTdphidy_diff,  vnpT_cos_diff, vnpT_sin_diff,
                             vn_cos_diff,  vn_sin_diff,
                              dNd2Mdy_diff,  dNd2pTd2M_diff,
                             dNd2pTd2Mdy_diff, vnMpTy_cos_diff, vnMpTy_sin_diff);
    

    calPhoton_SpvnpT( dNd2pTdphidy_tot,  vnpT_cos_tot, vnpT_sin_tot,
                                 vn_cos_tot, vn_sin_tot,
                                  dNd2Mdy_tot, dNd2pTd2M_tot,
                                 dNd2pTd2Mdy_tot,vnMpTy_cos_tot,vnMpTy_sin_tot);
    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_theta,dNd2pTd2M_pol_lambda_theta,dNd2Mdy_pol_lambda_theta,dNd2pTd2Mdy_pol_lambda_theta);
    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_norm,dNd2pTd2M_pol_lambda_norm,dNd2Mdy_pol_lambda_norm,dNd2pTd2Mdy_pol_lambda_norm);
    calPhoton_SpvnpT_pol(dNd2pTdphidy_pol_lambda_phi,dNd2pTd2M_pol_lambda_phi,dNd2Mdy_pol_lambda_phi,dNd2pTd2Mdy_pol_lambda_phi);

}

// pT-integrated and pT-differential spectra and flows for individual channels
void ThermalPhoton::calPhoton_SpvnpT_pol(double ****dNd2pTdphidy_pol_lambda_theta, double **dNd2pTd2M_pol_lambda_theta, double*dNd2Mdy_pol_lambda_theta, double*** dNd2pTd2Mdy_pol_lambda_theta ){



  for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
          double p_local = getPhotonp(i); // p_T
          double pweight_local = getPhoton_pweight(i);
          for (int j = 0; j < nphi; j++) {
              double phi_local = getPhotonphi(j); // phi
              double phi_weight_local = getPhoton_phiweight(j);
              
              for (int k = 0; k < nrapidity; k++) {
                double y_weight_local = getPhoton_yweight(k); // y
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
      double y_weight_local = getPhoton_yweight(k); // y
      for (int m = 0; m < nm; m++) {
          for (int i = 0; i < np; i++) {
              double p_local = getPhotonp(i); // p_T
              double pweight_local = getPhoton_pweight(i);
              
              for (int j = 0; j < nphi; j++) {
                  double phi_local = getPhotonphi(j); // phi
                  double phi_weight_local = getPhoton_phiweight(j);
                  dNd2pTd2Mdy_pol_lambda_theta[m][i][k] += dNd2pTdphidy_pol_lambda_theta[m][i][j][k]*phi_weight_local;
              
              }
          }
      }
  }

  
}

void ThermalPhoton::calPhoton_SpvnpT(
        double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
        double **vn_cos_eq, double **vn_sin_eq,
        double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
        double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq) {
    

       
  for (int m = 0; m < nm; m++) {
    for (int i = 0; i < np; i++) {
      double p_local = getPhotonp(i); // p_T
      double pweight_local = getPhoton_pweight(i);
      for (int j = 0; j < nphi; j++) {
        double phi_local = getPhotonphi(j); // phi
        double phi_weight_local = getPhoton_phiweight(j);
        for (int k = 0; k < nrapidity; k++) {
          double y_weight_local = getPhoton_yweight(k); // y
          // below dy/Dy to make everything rapidity density, Dy is the rapidity
          // range
          double weight = phi_weight_local * y_weight_local * dy / Dy;
          // integrate over rapidity, azimuthal angle of momentum
          dNd2pTd2M_eq[m][i] += dNd2pTdphidy_eq[m][i][j][k] * weight;
          // //if( std::isnan(dNd2pTdphidy_eq[m][i][j][k] ))
          // {
          //   std::cout<< dNd2pTdphidy_eq[m][i][j][k]  << " err "<<dNd2pTd2M_eq[m][i] <<std::endl;
          //   //exit(1);
          // }

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
    double y_weight_local = getPhoton_yweight(k); // y
    for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
        double p_local = getPhotonp(i); // p_T
        double pweight_local = getPhoton_pweight(i);

        for (int j = 0; j < nphi; j++) {
          double phi_local = getPhotonphi(j); // phi
          double phi_weight_local = getPhoton_phiweight(j);

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

void ThermalPhoton::outputPhoton_SpvnpT(string path, string type_str,string type_str2,
        double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
        double **vn_cos_eq, double **vn_sin_eq,
        double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
        double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq) {

          ostringstream filename_stream_eq_SpMatrix;
          ostringstream filename_stream_eq_Spvn;
          ostringstream filename_stream_eq_inte_Spvn;
          ostringstream filename_stream_eq_SpMatrix_dy;
          string file_end = ".dat";
          filename_stream_eq_SpMatrix << path << emissionProcess_name << "_"<< type_str<<"_"<<type_str2 <<"_SpMatrix"
          << file_end;
          filename_stream_eq_Spvn << path << emissionProcess_name << "_"<< type_str<<"_"<<type_str2 << "_Spvn" << file_end;
          filename_stream_eq_inte_Spvn << path << emissionProcess_name << "_"<< type_str<<"_"<<type_str2 << "_Spvn_inte"
          << file_end;
          filename_stream_eq_SpMatrix_dy << path << emissionProcess_name << "_"<< type_str<<"_"<<type_str2 << "_Spvn_MpTy"
           << file_end;

          ofstream fphoton_eq_SpMatrix_dy(filename_stream_eq_SpMatrix_dy.str().c_str());
          ofstream fphoton_eq_SpMatrix(filename_stream_eq_SpMatrix.str().c_str());
          ofstream fphoton_eq_Spvn(filename_stream_eq_Spvn.str().c_str());
          ofstream fphoton_eq_inte_Spvn(filename_stream_eq_inte_Spvn.str().c_str());

          
          for (int m = 0; m < nm; m++) {
            for (int i = 0; i < nphi; i++) {
              double phi_local = getPhotonphi(i);
              fphoton_eq_SpMatrix << phi_local << "  ";
              for (int j = 0; j < np; j++) {
                double temp_eq = 0.0; 
                for (int k = 0; k < nrapidity; k++) {
                  double y_weight_local = getPhoton_yweight(k);
                  
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
                double mll_local = getDileptonMass(m); // Mll
                double p_local = getPhotonp(i);        // p_T
                double y_local = getPhotonrapidity(k);
        
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
          double M_ll_local = getDileptonMass(m);
          double pT_local = getPhotonp(i);
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
    double M_ll_local = getDileptonMass(m);
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


  fphoton_eq_SpMatrix_dy.close();
  fphoton_eq_SpMatrix.close();
  fphoton_eq_Spvn.close();
  fphoton_eq_inte_Spvn.close();



}


void ThermalPhoton::outputPhoton_SpvnpT_pol(std::string path, std::string type_str, std::string type_str2, double **** dNd2pTdphidy_pol_lambda_theta,double *** dNd2pTd2Mdy_pol_lambda_theta, double** dNd2pTd2M_pol_lambda_theta,double* dNd2Mdy_pol_lambda_theta, double* dNd2Mdy_pol_lambda_norm) {

  ostringstream filename_stream_pol_lambda_theta_SpMatrix;
  ostringstream filename_stream_pol_lambda_theta_Spvn;
  ostringstream filename_stream_pol_lambda_theta_inte_Spvn;
  
  // ostringstream filename_stream_pol_lambda_norm_SpMatrix;
  // ostringstream filename_stream_pol_lambda_norm_Spvn;
  // ostringstream filename_stream_pol_lambda_norm_inte_Spvn;

  // ostringstream filename_stream_pol_lambda_phi_SpMatrix;
  // ostringstream filename_stream_pol_lambda_phi_Spvn;
  // ostringstream filename_stream_pol_lambda_phi_inte_Spvn;

  //rapdity
  ostringstream filename_stream_pol_lambda_theta_SpMatrix_dy;
  // ostringstream filename_stream_pol_lambda_norm_SpMatrix_dy;
  // ostringstream filename_stream_pol_lambda_phi_SpMatrix_dy;


  string file_end = ".dat";
  filename_stream_pol_lambda_theta_SpMatrix << path << emissionProcess_name << "_"<< type_str <<"_"<<type_str2 <<"_SpMatrix" << file_end;
  filename_stream_pol_lambda_theta_Spvn << path << emissionProcess_name << "_"<< type_str<<"_"<<type_str2 << "_Spvn" << file_end;
  filename_stream_pol_lambda_theta_inte_Spvn << path << emissionProcess_name << "_"<< type_str<<"_"<<type_str2 << "_Spvn_inte"<< file_end;
  filename_stream_pol_lambda_theta_SpMatrix_dy << path << emissionProcess_name << "_"<< type_str<<"_"<<type_str2 << "_Spvn_MpTy"<< file_end;



  ofstream fphoton_pol_lambda_theta_SpMatrix_dy(filename_stream_pol_lambda_theta_SpMatrix_dy.str().c_str());
  ofstream fphoton_pol_lambda_theta_SpMatrix(filename_stream_pol_lambda_theta_SpMatrix.str().c_str());
  ofstream fphoton_pol_lambda_theta_Spvn(filename_stream_pol_lambda_theta_Spvn.str().c_str());
  ofstream fphoton_pol_lambda_theta_inte_Spvn(filename_stream_pol_lambda_theta_inte_Spvn.str().c_str());

  
  for (int m = 0; m < nm; m++) {
      for (int i=0; i < nphi; i++) {
          double phi_local = getPhotonphi(i);
          fphoton_pol_lambda_theta_SpMatrix << phi_local << "  ";
         
          for (int j = 0; j < np; j++) {
              double temp_pol_lambda_theta = 0.0;

              for (int k = 0; k < nrapidity; k++) {
                  double y_weight_local = getPhoton_yweight(k);
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
                  double mll_local = getDileptonMass(m); //Mll
                  double p_local = getPhotonp(i); //p_T
                  double y_local = getPhotonrapidity(k);
                  
                  
                  fphoton_pol_lambda_theta_SpMatrix_dy<< scientific << setprecision(6) << setw(16) 
                              << mll_local << "  "<<p_local<<" "<<y_local<<" "<< dNd2pTd2Mdy_pol_lambda_theta[m][i][k] << "  ";
                  fphoton_pol_lambda_theta_SpMatrix_dy<< endl;
              
              }
      }
  }


  // pT differential, dN/(2pi pTdpT MdM dy) and vn(M, pT)
  for (int m = 0; m < nm; m++) {
      for (int i = 0; i < np; i++) {
          double M_ll_local = getDileptonMass(m);
          double pT_local = getPhotonp(i);
         
          fphoton_pol_lambda_theta_Spvn << scientific << setprecision(6) << setw(16)
                          << M_ll_local << "  "<< pT_local << "  " << dNd2pTd2M_pol_lambda_theta[m][i] << "  ";
          fphoton_pol_lambda_theta_Spvn << endl;
          
      }
  }
  

  if( type_str != "pol_lambda_norm")
  {
  // pT integrated
  for (int m = 0; m < nm; m++) {
      double M_ll_local = getDileptonMass(m);
      // get dN/dMdy from dN/MdMdy
      fphoton_pol_lambda_theta_inte_Spvn << scientific << setprecision(6) << setw(16)
                      << M_ll_local << "  " << dNd2Mdy_pol_lambda_theta[m]/dNd2Mdy_pol_lambda_norm[m] << "  ";
      fphoton_pol_lambda_theta_inte_Spvn << endl;

  }
  }


  fphoton_pol_lambda_theta_SpMatrix_dy.close();
  fphoton_pol_lambda_theta_SpMatrix.close();
  fphoton_pol_lambda_theta_Spvn.close();
  fphoton_pol_lambda_theta_inte_Spvn.close();



}



void ThermalPhoton::outputPhoton_SpvnpT_shell(string path,std::string type_str) {


    outputPhoton_SpvnpT(path, "eq", type_str, dNd2pTdphidy_eq, vnpT_cos_eq,vnpT_sin_eq,
      vn_cos_eq, vn_sin_eq,
      dNd2Mdy_eq, dNd2pTd2M_eq,
      dNd2pTd2Mdy_eq,vnMpTy_cos_eq,vnMpTy_sin_eq);
          
          
      outputPhoton_SpvnpT(path, "eqT", type_str,dNd2pTdphidy_eqT,  vnpT_cos_eqT, vnpT_sin_eqT,
                    vn_cos_eqT,  vn_sin_eqT,
                    dNd2Mdy_eqT,  dNd2pTd2M_eqT,
                    dNd2pTd2Mdy_eqT, vnMpTy_cos_eqT, vnMpTy_sin_eqT);

      outputPhoton_SpvnpT(path, "eqL", type_str, dNd2pTdphidy_eqL, vnpT_cos_eqL, vnpT_sin_eqL,
            vn_cos_eqL,  vn_sin_eqL,
            dNd2Mdy_eqL,dNd2pTd2M_eqL,
            dNd2pTd2Mdy_eqL, vnMpTy_cos_eqL, vnMpTy_sin_eqL);
          
      outputPhoton_SpvnpT(path, "visc", type_str,dNd2pTdphidy_visc, vnpT_cos_visc, vnpT_sin_visc,
      vn_cos_visc,  vn_sin_visc,
      dNd2Mdy_visc,  dNd2pTd2M_visc,
      dNd2pTd2Mdy_visc, vnMpTy_cos_visc, vnMpTy_sin_visc);

      outputPhoton_SpvnpT(path, "diff", type_str,dNd2pTdphidy_diff,  vnpT_cos_diff, vnpT_sin_diff,
          vn_cos_diff,  vn_sin_diff,
          dNd2Mdy_diff,  dNd2pTd2M_diff,
          dNd2pTd2Mdy_diff, vnMpTy_cos_diff, vnMpTy_sin_diff);


      outputPhoton_SpvnpT(path, "tot", type_str, dNd2pTdphidy_tot,  vnpT_cos_tot, vnpT_sin_tot,
              vn_cos_tot, vn_sin_tot,
              dNd2Mdy_tot, dNd2pTd2M_tot,
              dNd2pTd2Mdy_tot,vnMpTy_cos_tot,vnMpTy_sin_tot);
      
      
      outputPhoton_SpvnpT_pol(path, "pol_lambda_theta", type_str, dNd2pTdphidy_pol_lambda_theta, dNd2pTd2Mdy_pol_lambda_theta, dNd2pTd2M_pol_lambda_theta, dNd2Mdy_pol_lambda_theta,  dNd2Mdy_pol_lambda_norm);
      outputPhoton_SpvnpT_pol(path, "pol_lambda_phi", type_str, dNd2pTdphidy_pol_lambda_phi, dNd2pTd2Mdy_pol_lambda_phi, dNd2pTd2M_pol_lambda_phi, dNd2Mdy_pol_lambda_phi,  dNd2Mdy_pol_lambda_norm);
      outputPhoton_SpvnpT_pol(path, "pol_lambda_norm", type_str, dNd2pTdphidy_pol_lambda_norm, dNd2pTd2Mdy_pol_lambda_norm, dNd2pTd2M_pol_lambda_norm, dNd2Mdy_pol_lambda_norm,  dNd2Mdy_pol_lambda_norm);


}
