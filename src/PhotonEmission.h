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

  double **lambda; // Lorentz boost



  double ******dNd2pTdphidydTdtau_eq, ******dNd2pTdphidydTdtau_tot;
  double ******dNd2pTdphidydTdtau_visc, ******dNd2pTdphidydTdtau_diff;

  double ***dNd2pTdphidydTdtau_eq_all, ***dNd2pTdphidydTdtau_tot_all;
  double ***dNd2pTdphidydTdtau_visc_all, ***dNd2pTdphidydTdtau_diff_all;



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
  void outputPhoton_total_SpMatrix_and_SpvnpT(int hydro_mode = -1);
  void outputPhotonSpvn_individualchannel();
  double suppression_factor(double tau, double T);
};

#endif // SRC_PHOTONEMISSION_H_
