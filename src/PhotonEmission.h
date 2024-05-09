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
    double Dy;

    double gridDx, gridDy, gridDtau;
    double gridX0, gridY0;
    double ETAmax;

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

    double  sigmael_over_T;
    double  eB_over_mpi2;
    int include_EM_deltaf;
    

    double **lambda;  // Lorentz boost

    double **dNd2pTd2M_eq;
    double **dNd2pTd2M_eqT;
    double **dNd2pTd2M_eqL;
    double **dNd2pTd2M_visc;
    double **dNd2pTd2M_diff;
    double **dNd2pTd2M_em;
    double **dNd2pTd2M_em_E;
    
    double **dNd2pTd2M_tot;
    double ****dNd2pTdphidy_eq;
    double ****dNd2pTdphidy_eqT;
    double ****dNd2pTdphidy_eqL;
    double ***vnpT_cos_eq, ***vnpT_sin_eq;
    double ****dNd2pTdphidy_visc;
    double ***vnpT_cos_visc, ***vnpT_sin_visc;
    double ****dNd2pTdphidy_diff;
    double ***vnpT_cos_diff, ***vnpT_sin_diff;
    
    double ****dNd2pTdphidy_em;
    double ***vnpT_cos_em, ***vnpT_sin_em;
    
    double ****dNd2pTdphidy_em_E;
    double ***vnpT_cos_em_E, ***vnpT_sin_em_E;
    

    double ****dNd2pTdphidy_tot;
    double ***vnpT_cos_tot, ***vnpT_sin_tot;

    std::vector<double> dNd2Mdy_eq, dNd2Mdy_visc, dNd2Mdy_diff, dNd2Mdy_em, dNd2Mdy_em_E, dNd2Mdy_tot;
    std::vector<double> dNd2Mdy_eqT, dNd2Mdy_eqL;
    double **vn_sin_eq;
    double **vn_cos_eq;
    double **vn_sin_visc;
    double **vn_cos_visc;
    double **vn_sin_diff;
    double **vn_cos_diff;

    double **vn_sin_em;
    double **vn_cos_em;
    
    double **vn_sin_em_E;
    double **vn_cos_em_E;

    double **vn_cos_tot;
    double **vn_sin_tot;

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
    void calPhotonemission_3d(void *hydroinfo_ptr_in,int hydro_mode=-1);
    void calPhoton_total_Spvn();
    void calPhoton_total_Spvn_sum(const PhotonEmission& spvn_tem); 
    void calPhoton_SpvnpT_individualchannel();
    void outputPhoton_total_SpMatrix_and_SpvnpT(int hydro_mode=-1);
    void outputPhotonSpvn_individualchannel();
    double suppression_factor(double tau,double T);
    void EM_profile_0(double sigma_el,double local_t, double local_x, double local_y, double local_z, double& eB, double& eEx, double& eEz);
};

#endif   // SRC_PHOTONEMISSION_H_

