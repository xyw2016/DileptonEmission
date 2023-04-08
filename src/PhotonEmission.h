#ifndef SRC_PHOTONEMISSION_H_
#define SRC_PHOTONEMISSION_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <memory>
#include <vector>

#include "Hydroinfo_h5.h"
#include "Hydroinfo_MUSIC.h"
#include "ThermalPhoton.h"
#include "ParameterReader.h"


class PhotonEmission {
 private:
    std::shared_ptr<ParameterReader> paraRdr;
    std::string output_path;
    //photon production processes
    std::unique_ptr<ThermalPhoton> dilepton_QGP_thermal;

    int CORES;

    int neta;
    int nm;
    int np, nphi, nrapidity;
    int norder;

    double gridDx, gridDy, gridDtau;
    double gridX0, gridY0;

    double T_dec, T_sw_high, T_sw_low;
    double T_cuthigh, T_cutlow;
    double tau_cut_high, tau_cut_low;
    int nTcut, n_tau_cut;

    double T_test;
    double muB_test;
    double rhoB_eplusp_test;

    int hydro_flag;
    int differential_flag;
    int turn_off_transverse_flow;
    int turn_on_muB_;
    int calHGIdFlag;
    int diff_flag;
    int emission_rate_flag;
    int test_code_flag;

    double **lambda;  // Lorentz boost transverse only

    double **dNd2pTd2M_eq;
    double **dNd2pTd2M_tot;
    double ****dNd2pTdphidy_eq;
    double ***vnpT_cos_eq, ***vnpT_sin_eq;
    double ****dNd2pTdphidy_tot;
    double ***vnpT_cos_tot, ***vnpT_sin_tot;

    std::vector<double> dNd2Mdy_eq, dNd2Mdy_tot;
    double **vn_sin_eq;
    double **vn_cos_eq;
    double **vn_cos_tot;
    double **vn_sin_tot;

    double ******dNd2pTdphidydTdtau_eq, ******dNd2pTdphidydTdtau_tot;
    double ******dNd2pTdphidydTdtau_diff;

    double ***dNd2pTdphidydTdtau_eq_all, ***dNd2pTdphidydTdtau_tot_all;
    double ***dNd2pTdphidydTdtau_diff_all;

 public:
    PhotonEmission(std::shared_ptr<ParameterReader> paraRdr_in);
    ~PhotonEmission();

    void set_hydroGridinfo();
    void print_hydroGridinfo();
    void InitializePhotonEmissionRateTables();
    void calPhotonemission_3d(void *hydroinfo_ptr_in);
    void calPhoton_total_Spvn();
    void calPhoton_SpvnpT_individualchannel();
    void outputPhoton_total_SpMatrix_and_SpvnpT();
    void outputPhotonSpvn_individualchannel();
};

#endif   // SRC_PHOTONEMISSION_H_

