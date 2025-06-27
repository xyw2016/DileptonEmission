#ifndef SRC_THERMALPHOTON_H_
#define SRC_THERMALPHOTON_H_

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Arsenal.h"
#include "ParameterReader.h"
#include "Table2D.h"

class ThermalPhoton {
private:
  std::shared_ptr<ParameterReader> paraRdr;

  int np, nphi, nrapidity;
  int norder;
  int neta;
  int nm;
  std::string rate_path_;

  int use_logarithmic_mass_grid;

  int CORES;

  int CS_frame;

  double dy;
  double dM;
  double Dy;

  bool bRateTable_;
  bool bShearVisCorr_;
  bool bBulkVisCorr_;
  bool bDiffusionCorr_;

  int include_diff_deltaf;
  int include_visc_deltaf;
  int turn_on_muB_;


  // double em_eqrate;
  // double em_eqrateT;
  // double em_eqrateL;
  // double em_visrate;
  // double em_bulkvis;
  // double em_diffrate;




  double alpha_s;

  // photon emission rate
  std::unique_ptr<Table2D> Photonemission_eqrateTable_ptr;
  std::unique_ptr<Table2D> Photonemission_viscous_rateTable_ptr;
  std::unique_ptr<Table2D> Photonemission_bulkvis_rateTable_ptr;

  double **Emission_eqrateTb_ptr;
  double **Emission_viscous_rateTb_ptr;
  double **Emission_bulkvis_rateTb_ptr;
  std::vector<double> EmissionrateTb_Yidxptr;
  double EmissionrateTb_Xmin;
  double EmissionrateTb_Ymin;
  int EmissionrateTb_sizeX;
  int EmissionrateTb_sizeY;
  double EmissionrateTb_dX;
  double EmissionrateTb_dY;

  std::vector<double> a_list; // alpha_s values
  std::vector<double> B_list; // chemical potential (muB/T)
  std::vector<double> M_list; // invariant mass (units of T)
  std::vector<double>
      k_list; // 3-momentum (k/T), defined in the local rest frame
  std::vector<double> rhoT_list, rhoL_list; // rho list

  // photon spectra parameters
  std::string emissionProcess_name;
  double *p, *p_weight;
  double *phi, *phi_weight;
  double *y_weight;
  std::vector<double> M;
  std::vector<double> y;



  double ***vnypT_cos_eq, ***vnypT_sin_eq;

  double **vnpT_cos_vis, **vnpT_sin_vis;
  double ***vnypT_cos_vis, ***vnypT_sin_vis;

  double **vnpT_cos_bulkvis, **vnpT_sin_bulkvis;
  double ***vnypT_cos_bulkvis, ***vnypT_sin_bulkvis;





  // matrix for cuts on temperature and proper time
  int nTcut, n_tau_cut;
  double ******dNd2pTdphidydTdtau_eq, ******dNd2pTdphidydTdtau_tot;
  double ******dNd2pTdphidydTdtau_visc, *****dNd2pTdphidydTdtau_bulkvis;
  double ******dNd2pTdphidydTdtau_diff;
  double ****dNpTdpTdydTdtau_eq, ****dNpTdpTdydTdtau_visc;
  double ****dNpTdpTdydTdtau_diff, ****dNpTdpTdydTdtau_tot;
  double ***dNdydTdtau_eq, ***dNdydTdtau_tot;
  double ***dNdydTdtau_visc; //, **dNdydTdtau_bulkvis;
  double ***dNdydTdtau_diff;
  double ****vndTdtau_cos_eq, ****vndTdtau_sin_eq;
  double ****vndTdtau_cos_visc, ****vndTdtau_sin_visc;
  // double ***vndTdtau_cos_bulkvis, ***vndTdtau_sin_bulkvis;
  double ****vndTdtau_cos_diff, ****vndTdtau_sin_diff;
  double ****vndTdtau_cos_tot, ****vndTdtau_sin_tot;

  // yields and spectra in temperature or proper time
  double ***dNpTdpTdydT_eq, ***dNpTdpTdydT_visc;
  double ***dNpTdpTdydT_diff, ***dNpTdpTdydT_tot;

  double ***dNpTdpTdydtau_eq, ***dNpTdpTdydtau_visc;
  double ***dNpTdpTdydtau_diff, ***dNpTdpTdydtau_tot;

  double **dNdydT_eq;
  double **dNdydT_visc;
  double **dNdydT_diff;
  double **dNdydT_tot;

  double **dNdydtau_eq;
  double **dNdydtau_visc;
  double **dNdydtau_diff;
  double **dNdydtau_tot;

//xyw
    double *dNd2pTdphidy_eq_all;
    double *dNd2pTdphidy_eqT_all;   
    double *dNd2pTdphidy_eqL_all;
    double *dNd2pTdphidy_visc_all;
    double *dNd2pTdphidy_diff_all; 
    double *dNd2pTdphidy_tot_all;
    double *dNd2pTdphidy_pol_lambda_theta_all;
    double *dNd2pTdphidy_pol_lambda_norm_all;
    double *dNd2pTdphidy_pol_lambda_phi_all;

    double ****dNd2pTdphidy_eq;
    double ****dNd2pTdphidy_eqT;
    double ****dNd2pTdphidy_eqL;
    double ****dNd2pTdphidy_visc;
    double ****dNd2pTdphidy_diff;
    double ****dNd2pTdphidy_tot;
    double ****dNd2pTdphidy_pol_lambda_theta;
    double ****dNd2pTdphidy_pol_lambda_norm;
    double ****dNd2pTdphidy_pol_lambda_phi;


    double ***dNd2pTd2Mdy_eq;
    double ***dNd2pTd2Mdy_eqT;
    double ***dNd2pTd2Mdy_eqL;
    double ***dNd2pTd2Mdy_visc;
    double ***dNd2pTd2Mdy_diff;
    double ***dNd2pTd2Mdy_tot;

    



    double ***vnpT_cos_eq, ***vnpT_sin_eq;
    double ***vnpT_cos_eqT, ***vnpT_sin_eqT;
    double ***vnpT_cos_eqL, ***vnpT_sin_eqL;
    double ***vnpT_cos_visc, ***vnpT_sin_visc;
    double ***vnpT_cos_diff, ***vnpT_sin_diff;
    double ***vnpT_cos_tot, ***vnpT_sin_tot;
   
    double ****vnMpTy_cos_eq, ****vnMpTy_sin_eq;
    double ****vnMpTy_cos_eqT, ****vnMpTy_sin_eqT;
    double ****vnMpTy_cos_eqL, ****vnMpTy_sin_eqL;
    double ****vnMpTy_cos_visc, ****vnMpTy_sin_visc;
    double ****vnMpTy_cos_diff, ****vnMpTy_sin_diff;
    double ****vnMpTy_cos_tot, ****vnMpTy_sin_tot;

    double **vn_sin_eq, **vn_cos_eq;
    double **vn_sin_eqT, **vn_cos_eqT;
    double **vn_sin_eqL, **vn_cos_eqL;
    double **vn_sin_visc, **vn_cos_visc;
    double **vn_sin_diff, **vn_cos_diff;
    double **vn_cos_tot, **vn_sin_tot;

    double **dNd2pTd2M_eq;
    double **dNd2pTd2M_eqT;
    double **dNd2pTd2M_eqL;
    double **dNd2pTd2M_visc;
    double **dNd2pTd2M_diff;
    double **dNd2pTd2M_tot;

    double ***dNd2pTd2Mdy_pol_lambda_theta;
    double ***dNd2pTd2Mdy_pol_lambda_norm;
    double ***dNd2pTd2Mdy_pol_lambda_phi;
    double **dNd2pTd2M_pol_lambda_theta;
    double **dNd2pTd2M_pol_lambda_norm;
    double **dNd2pTd2M_pol_lambda_phi;

    
    
    double * dNd2Mdy_eq;
    double *  dNd2Mdy_visc;
    double *  dNd2Mdy_diff;
    double *  dNd2Mdy_tot;
    double * dNd2Mdy_pol_lambda_norm;
    double *  dNd2Mdy_pol_lambda_theta;
    double *  dNd2Mdy_pol_lambda_phi;
    double * dNd2Mdy_eqT;
    double * dNd2Mdy_eqL;




public:
  ThermalPhoton(std::shared_ptr<ParameterReader> paraRdr_in,
                std::string emissionProcess);

  virtual ~ThermalPhoton();

  void readEmissionrateFromFile(bool bRateTable);
  void initialize(std::string fname, std::vector<double> &a_list,
                  std::vector<double> &B_list, std::vector<double> &M_list,
                  std::vector<double> &k_list, std::vector<double> &rhoT_list,
                  std::vector<double> &rhoL_list);

  double get_dy() { return (dy); }
  double get_Dy() { return (Dy); }

  double getPhotonp(int i) { return (p[i]); }
  double getPhoton_pweight(int i) { return (p_weight[i]); }
  double getPhotonphi(int i) { return (phi[i]); }
  double getPhoton_phiweight(int i) { return (phi_weight[i]); }
  double getDileptonMass(int i) { return (M[i]); }
  // double getPhotontheta(int i) {return(theta[i]);}
  double getPhotonrapidity(int i) { return (y[i]); }
  double getPhoton_yweight(int i) { return (y_weight[i]); }

  inline double get_dNd2pTdphidy_eq(int m, int i, int j, int k) const {
    return dNd2pTdphidy_eq[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_eqT(int m, int i, int j, int k) const {
    return dNd2pTdphidy_eqT[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_eqL(int m, int i, int j, int k) const {
    return dNd2pTdphidy_eqL[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_visc(int m, int i, int j, int k) const {
    return dNd2pTdphidy_visc[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_diff(int m, int i, int j, int k) const {
    return dNd2pTdphidy_diff[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_tot(int m, int i, int j, int k) const {
    return dNd2pTdphidy_tot[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_pol_lambda_theta(int m, int i, int j, int k) const {
    return dNd2pTdphidy_pol_lambda_theta[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_pol_lambda_norm(int m, int i, int j, int k) const {
    return dNd2pTdphidy_pol_lambda_norm[m][i][j][k];
  }

  inline double get_dNd2pTdphidy_pol_lambda_phi(int m, int i, int j, int k) const {
    return dNd2pTdphidy_pol_lambda_phi[m][i][j][k];
  }
  

  virtual void analyticRates(double T, double muB, std::vector<double> &Eq,
                             double *M_ll, std::vector<double> &eqrate_ptr,
                             int nm, int np, int nphi, int nrapidity);
  virtual void analyticRatesShearVis(double T, std::vector<double> &Eq,
                                     double *M_ll,
                                     std::vector<double> &eqrate_ptr);
  virtual void analyticRatesBulkVis(double T, std::vector<double> &Eq,
                                    double *M_ll,
                                    std::vector<double> &eqrate_ptr);
  virtual void FiniteBaryonRates(double T, double muB, double inv_eplusp,
                                 double rhoB_over_eplusp, double Eq,
                                 double M_ll, double &eqrate_ptr,
                                 double &eqrateT_ptr, double &eqrateL_ptr,
                                 double &viscrate_ptr, double &diffrate_ptr,
                                 int include_visc_deltaf,
                                 int include_diff_deltaf);

  void getPhotonemissionRate(double &Eq, double &M_ll, double &pi_zz,
                             double &bulkPi, double &diff_factor, double &T,
                             double &muB, double &inv_eplusp,
                             double &rhoB_over_eplusp, double& em_eqrate ,double& em_eqrateT,
                             double &em_eqrateL, double& em_visrate,double& em_bulkvis, double& em_diffrate);
                             
  void calThermalPhotonemission_3d(
      double (&p_lab_Min)[4], double (&flow_u_mu_Min)[4], double Eq,
      double M_ll, double pi_zz, double bulkPi, double diff_factor, double T,
      double muB, double inv_eplusp, double rhoB_over_eplusp, double volume,
      double fraction, double &dNd2pTdphidy_cell_eq,
      double &dNd2pTdphidy_cell_eqT, double &dNd2pTdphidy_cell_eqL,
      double &dNd2pTdphidy_cell_visc, double &dNd2pTdphidy_cell_diff,
      double &dNd2pTdphidy_cell_tot, double &dNd2pTdphidy_cell_lambda_norm,
      double &dNd2pTdphidy_cell_lambda_theta,
      double &dNd2pTdphidy_cell_lambda_phi);
  void calThermalPhotonemission_3d(double (&p_lab_Min)[4],
                                   double (&flow_u_mu_Min)[4], double &Eq,
                                   double &M_ll, double &pi_zz, double &bulkPi,
                                   double &diff_factor, double &T, double &muB,
                                   double &inv_eplusp, double &rhoB_over_eplusp,
                                   double &volume, double &fraction,long& n, int& i3);

  void calThermalPhotonemission_3d();

  void calPhoton_SpvnpT( double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
    double **vn_cos_eq, double **vn_sin_eq,
    double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
    double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq);

  void calPhoton_SpvnpT_pol(double ****dNd2pTdphidy_pol_lambda_theta, double **dNd2pTd2M_pol_lambda_theta, double*dNd2Mdy_pol_lambda_theta, double*** dNd2pTd2Mdy_pol_lambda_theta);

  void calPhoton_SpvnpT_shell();
  void calPhoton_Spvn_dTdtau();
  void calPhoton_Spectra_dTdtau();
  void calPhoton_SpMatrix_dTdtau(double ******dNd2pTdphidydTdtau_eq,
                                 double ******dNd2pTdphidydTdtau_visc,
                                 double ******dNd2pTdphidydTdtau_diff,
                                 double ******dNd2pTdphidydTdtau_tot);
  void outputPhoton_SpvnpT(std::string path, std::string type_str,std::string type_str2,
                            double ****dNd2pTdphidy_eq, double ***vnpT_cos_eq,double ***vnpT_sin_eq,
                            double **vn_cos_eq, double **vn_sin_eq,
                            double * dNd2Mdy_eq, double **dNd2pTd2M_eq,
                            double ***dNd2pTd2Mdy_eq,double ****vnMpTy_cos_eq,double ****vnMpTy_sin_eq);
  void outputPhoton_SpvnpT_pol(std::string path,std::string type_str,std::string type_str2, double **** dNd2pTdphidy_pol_lambda_theta,double *** dNd2pTd2Mdy_pol_lambda_theta, double** dNd2pTd2M_pol_lambda_theta,double* dNd2Mdy_pol_lambda_theta, double* dNd2Mdy_pol_lambda_norm);
  void outputPhoton_SpvnpT_shell(std::string path,std::string type_str);
  void outputPhoton_Spvn_dTdtau(std::string path, double Tcut_high,
                                double Tcut_low, double tau_cut_high,
                                double tau_cut_low);
  void outputPhoton_Spectra_dTdtau(std::string path, double Tcut_high,
                                   double Tcut_low, double tau_cut_high,
                                   double tau_cut_low);
  void outputPhoton_Spectra_full_diff(std::string path, double Tcut_high,
                                      double Tcut_low, double tau_cut_high,
                                      double tau_cut_low);

  struct Table {
    int nx, ny, nz, nw;
    double *x, *y, *z, *w, ****F;
    double x_min, x_max, y_min, y_max, z_min, z_max, w_min, w_max;

    Table() = default;
    Table(std::vector<double> &_x, std::vector<double> &_y,
          std::vector<double> &_z, std::vector<double> &_w,
          std::vector<double> &_F);
    double interp(double _x, double _y, double _z, double _w);
    // ~Table();
  };

  void NLO_rate(struct Table grid_T, struct Table grid_L, double o, double k,
                double alpha_s, double muB, double T, double m_l,
                double &rateTot, double &rateT, double &rateL);
   void reduce_multile_core();
   virtual void getRateFromTable(const double E,
    const double T_local, const double k_local, const double M_local, double &rateTot,
    double &rateT, double &rateL){
      rateTot = 0.0;
      rateT = 0.0;
      rateL = 0.0;
  };

private:
  Table grid_T;
  Table grid_L;
};
#endif // SRC_THERMALPHOTON_H_
