#ifndef SRC_PHOTONEMISSION_H_
#define SRC_PHOTONEMISSION_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Hydroinfo_MUSIC.h"
#include "Hydroinfo_h5.h"
#include "ParameterReader.h"
#include "ThermalPhoton.h"

class PhotonEmission {
private:
  std::shared_ptr<ParameterReader> paraRdr;
  std::string output_path;
  // photon production processes
  std::unique_ptr<ThermalPhoton> dilepton_QGP_thermal_LO;
  std::unique_ptr<ThermalPhoton> dilepton_QGP_thermal;
  std::unique_ptr<ThermalPhoton> HadronGas_rho_meson;
  std::unique_ptr<ThermalPhoton> HadronGas_4piV;

  int CORES;

  int neta;
  int nm;
  int np, nphi, nrapidity;
  int norder;
  double Dy;

  double gridDx, gridDy, gridDtau;
  double gridX0, gridY0, gridTau0, gridTauf;
  double ETAmax;
  int gridNx, gridNy;

  double tau0, tau_max;

  double T_dec, T_sw_high, T_sw_low;
  double T_cuthigh, T_cutlow;
  double tau_cut_high, tau_cut_low;
  int nTcut, n_tau_cut;

  double T_test;
  double muB_test;
  double rhoB_eplusp_test;
  double inv_eplusp_test;

  int hydro_flag;
  int differential_flag;
  int turn_off_transverse_flow;
  int turn_on_muB_;
  int calHGIdFlag;
  int emission_rate_flag;
  int test_code_flag;
  int rho_rate_flag;

  double **lambda; // Lorentz boost



  double ******dNd2pTdphidydTdtau_eq, ******dNd2pTdphidydTdtau_tot;
  double ******dNd2pTdphidydTdtau_visc, ******dNd2pTdphidydTdtau_diff;

  double ***dNd2pTdphidydTdtau_eq_all, ***dNd2pTdphidydTdtau_tot_all;
  double ***dNd2pTdphidydTdtau_visc_all, ***dNd2pTdphidydTdtau_diff_all;

  //LO
  double ****dNd2pTdphidy_eq_lo;
  double ****dNd2pTdphidy_tot_lo;
  double ****dNd2pTdphidy_pol_lambda_theta_lo;
  double ****dNd2pTdphidy_pol_lambda_norm_lo;
  double ****dNd2pTdphidy_pol_lambda_phi_lo;


  double ***dNd2pTd2Mdy_eq_lo;
  double ***dNd2pTd2Mdy_tot_lo;


  double ***vnpT_cos_eq_lo, ***vnpT_sin_eq_lo;
  double ***vnpT_cos_tot_lo, ***vnpT_sin_tot_lo;
 
  double ****vnMpTy_cos_eq_lo, ****vnMpTy_sin_eq_lo;
  double ****vnMpTy_cos_tot_lo, ****vnMpTy_sin_tot_lo;

  double **vn_sin_eq_lo, **vn_cos_eq_lo;
  double **vn_cos_tot_lo, **vn_sin_tot_lo;

  double **dNd2pTd2M_eq_lo;
  double **dNd2pTd2M_tot_lo;

  double ***dNd2pTd2Mdy_pol_lambda_theta_lo;
  double ***dNd2pTd2Mdy_pol_lambda_norm_lo;
  double ***dNd2pTd2Mdy_pol_lambda_phi_lo;
  double **dNd2pTd2M_pol_lambda_theta_lo;
  double **dNd2pTd2M_pol_lambda_norm_lo;
  double **dNd2pTd2M_pol_lambda_phi_lo;

  
  
  double * dNd2Mdy_eq_lo;
  double *  dNd2Mdy_tot_lo;
  double * dNd2Mdy_pol_lambda_norm_lo;
  double *  dNd2Mdy_pol_lambda_theta_lo;
  double *  dNd2Mdy_pol_lambda_phi_lo;
 
  //NLO
  double ****dNd2pTdphidy_eq_nlo;
  double ****dNd2pTdphidy_tot_nlo;
  double ****dNd2pTdphidy_pol_lambda_theta_nlo;
  double ****dNd2pTdphidy_pol_lambda_norm_nlo;
  double ****dNd2pTdphidy_pol_lambda_phi_nlo;


  double ***dNd2pTd2Mdy_eq_nlo;
  double ***dNd2pTd2Mdy_tot_nlo;


  double ***vnpT_cos_eq_nlo, ***vnpT_sin_eq_nlo;
  double ***vnpT_cos_tot_nlo, ***vnpT_sin_tot_nlo;
 
  double ****vnMpTy_cos_eq_nlo, ****vnMpTy_sin_eq_nlo;
  double ****vnMpTy_cos_tot_nlo, ****vnMpTy_sin_tot_nlo;

  double **vn_sin_eq_nlo, **vn_cos_eq_nlo;
  double **vn_cos_tot_nlo, **vn_sin_tot_nlo;

  double **dNd2pTd2M_eq_nlo;
  double **dNd2pTd2M_tot_nlo;

  double ***dNd2pTd2Mdy_pol_lambda_theta_nlo;
  double ***dNd2pTd2Mdy_pol_lambda_norm_nlo;
  double ***dNd2pTd2Mdy_pol_lambda_phi_nlo;
  double **dNd2pTd2M_pol_lambda_theta_nlo;
  double **dNd2pTd2M_pol_lambda_norm_nlo;
  double **dNd2pTd2M_pol_lambda_phi_nlo;

  
  
  double * dNd2Mdy_eq_nlo;
  double *  dNd2Mdy_tot_nlo;
  double * dNd2Mdy_pol_lambda_norm_nlo;
  double *  dNd2Mdy_pol_lambda_theta_nlo;
  double *  dNd2Mdy_pol_lambda_phi_nlo;


public:
  PhotonEmission(std::shared_ptr<ParameterReader> paraRdr_in);
  ~PhotonEmission();

  void set_hydroGridinfo();
  void print_hydroGridinfo();
  void InitializePhotonEmissionRateTables();
  void calPhotonemission_3d(void *hydroinfo_ptr_in, int hydro_mode = -1);
  void calPhotonemission_2d(void *hydroinfo_ptr_in, int hydro_mode = -1);
  void calPhoton_total_Spvn();
  void calPhoton_total_Spvn_sum(const PhotonEmission &spvn_tem);
  void calPhoton_SpvnpT_individualchannel();
  void outputPhoton_total_SpMatrix_and_SpvnpT(std::string type_str);
  void outputPhotonSpvn_individualchannel(std::string type_str);
  double suppression_factor(double tau, double T);

  void calPhoton_SpvnpT( double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
    double **vn_cos_eq, double **vn_sin_eq,
    double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
    double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq);
    void calPhoton_total_Spvn_shell();

  void calPhoton_SpvnpT_pol(double ****dNd2pTdphidy_pol_lambda_theta, double **dNd2pTd2M_pol_lambda_theta, double*dNd2Mdy_pol_lambda_theta, double*** dNd2pTd2Mdy_pol_lambda_theta );

  void outputPhoton_SpvnpT_pol(std::string path, std::string type_str, double **** dNd2pTdphidy_pol_lambda_theta,double *** dNd2pTd2Mdy_pol_lambda_theta, double** dNd2pTd2M_pol_lambda_theta,double* dNd2Mdy_pol_lambda_theta, double* dNd2Mdy_pol_lambda_norm);

  void outputPhoton_SpvnpT(string path, string type_str,
    double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
    double **vn_cos_eq, double **vn_sin_eq,
    double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
    double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq);
  

 
};

#endif // SRC_PHOTONEMISSION_H_
