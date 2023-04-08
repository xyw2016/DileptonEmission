#ifndef SRC_THERMALPHOTON_H_
#define SRC_THERMALPHOTON_H_

#include <string>
#include <fstream>
#include <memory>
#include <vector>

#include "Arsenal.h"
#include "Table2D.h"
#include "ParameterReader.h"

class ThermalPhoton {
 private:
    std::shared_ptr<ParameterReader> paraRdr;

    int np, nphi, nrapidity;
    int norder;
    int neta;
    int nm;
    std::string rate_path_;

    double dy;
    double dM;

    bool bRateTable_;
    bool bShearVisCorr_;
    bool bBulkVisCorr_;
    bool bDiffusionCorr_;

    int include_diff_deltaf;
    int turn_on_muB_;

    double alpha_s;

    // photon emission rate
    std::unique_ptr<Table2D> Photonemission_eqrateTable_ptr;
    std::unique_ptr<Table2D> Photonemission_viscous_rateTable_ptr;
    std::unique_ptr<Table2D> Photonemission_bulkvis_rateTable_ptr;

    double** Emission_eqrateTb_ptr;
    double** Emission_viscous_rateTb_ptr;
    double** Emission_bulkvis_rateTb_ptr;
    std::vector<double> EmissionrateTb_Yidxptr;
    double EmissionrateTb_Xmin;
    double EmissionrateTb_Ymin;
    int EmissionrateTb_sizeX;
    int EmissionrateTb_sizeY;
    double EmissionrateTb_dX;
    double EmissionrateTb_dY;

    std::vector<double>  a_list; // alpha_s values
    std::vector<double>  B_list; // chemical potential (muB/T)
    std::vector<double>  M_list; // invariant mass (units of T)
    std::vector<double>  k_list; // 3-momentum (k/T), defined in the local rest frame
    std::vector<double>  rhoT_list, rhoL_list; // rho list

    // photon spectra parameters
    std::string emissionProcess_name;
    double *p, *p_weight;
    double *phi, *phi_weight;
    std::vector<double> M;
    std::vector<double> y;
    std::vector<double> theta;

    double ****dNd2pTdphidy_eq, ***dNd2pTdphidy_vis, ****dNd2pTdphidy_tot;
    double ***dNd2pTdphidy_bulkvis, ****dNd2pTdphidy_diff;

    double **vnpT_cos_eq, **vnpT_sin_eq;
    double ***vnypT_cos_eq, ***vnypT_sin_eq;

    double **vnpT_cos_vis, **vnpT_sin_vis;
    double ***vnypT_cos_vis, ***vnypT_sin_vis;

    double **vnpT_cos_bulkvis, **vnpT_sin_bulkvis;
    double ***vnypT_cos_bulkvis, ***vnypT_sin_bulkvis;

    double **vnpT_cos_tot, **vnpT_sin_tot;
    double ***vnypT_cos_tot, ***vnypT_sin_tot;

    std::vector<double> vn_cos_eq;
    std::vector<double> vn_sin_eq;
    std::vector<double> vn_cos_vis;
    std::vector<double> vn_sin_vis;
    std::vector<double> vn_cos_bulkvis;
    std::vector<double> vn_sin_bulkvis;
    std::vector<double> vn_cos_tot;
    std::vector<double> vn_sin_tot;

    // matrix for cuts on temperature and proper time
    int nTcut, n_tau_cut;
    double Tcut_high, Tcut_low;
    double tau_cut_high, tau_cut_low;
    double ******dNd2pTdphidydTdtau_eq, ******dNd2pTdphidydTdtau_tot;
    // double *****dNd2pTdphidydTdtau_vis, *****dNd2pTdphidydTdtau_bulkvis;
    double ******dNd2pTdphidydTdtau_diff;
    double ***dNdydTdtau_eq, ***dNdydTdtau_tot;
    // double **dNdydTdtau_vis, **dNdydTdtau_bulkvis;
    double ***dNdydTdtau_diff;
    double ****vndTdtau_cos_eq, ****vndTdtau_sin_eq;
    // double ***vndTdtau_cos_vis, ***vndTdtau_sin_vis;
    // double ***vndTdtau_cos_bulkvis, ***vndTdtau_sin_bulkvis;
    double ****vndTdtau_cos_diff, ****vndTdtau_sin_diff;
    double ****vndTdtau_cos_tot, ****vndTdtau_sin_tot;

 public:
    ThermalPhoton(std::shared_ptr<ParameterReader> paraRdr_in,
                  std::string emissionProcess);

    virtual ~ThermalPhoton();

    void readEmissionrateFromFile(bool bRateTable);
    void initialize(std::string fname, std::vector<double> &a_list, std::vector<double> &B_list, std::vector<double> &M_list, 
        std::vector<double> &k_list, std::vector<double> &rhoT_list, std::vector<double> &rhoL_list);

    double get_dy() {return(dy);}

    double getPhotonp(int i) {return(p[i]);}
    double getPhoton_pweight(int i) {return(p_weight[i]);}
    double getPhotonphi(int i) {return(phi[i]);}
    double getPhoton_phiweight(int i) {return(phi_weight[i]);}
    double getDileptonMass(int i) {return(M[i]);}
    double getPhotontheta(int i) {return(theta[i]);}
    double getPhotonrapidity(int i) {return(y[i]);}

    virtual void analyticRates(double T, double muB, std::vector<double> &Eq,
        double *M_ll, std::vector<double> &eqrate_ptr, int nm,
        int np, int nphi, int nrapidity);
    virtual void analyticRatesShearVis(double T, std::vector<double> &Eq,
        double *M_ll, std::vector<double> &eqrate_ptr);
    virtual void analyticRatesBulkVis(double T, std::vector<double> &Eq,
        double *M_ll, std::vector<double> &eqrate_ptr);
    virtual void FiniteBaryonRates(double T, double muB, double rhoB_over_eplusp, double Eq, 
    double M_ll, double &eqrate_ptr, double &diffrate_ptr, int include_diff_deltaf);

    void getPhotonemissionRate(double Eq, double M_ll, double pi_zz, double bulkPi,
        double diff_factor, double T, double muB, double rhoB_over_eplusp, double &eqrate_ptr, 
        double &visrate_ptr, double &bulkvis_ptr, double &diffrate_ptr);
    void calThermalPhotonemission_3d(double Eq, double M_ll, double pi_zz, double bulkPi, 
        double diff_factor, double T, double muB, double rhoB_over_eplusp, double volume, double fraction,
        double &dNd2pTdphidy_cell_eq, double &dNd2pTdphidy_cell_diff, double &dNd2pTdphidy_cell_tot);

    void calPhoton_SpvnpT(double ***dNd2pTdphipy,
                          double ***vnypT_cos, double *** vnypT_sin,
                          double **vnpT_cos, double **vnpT_sin,
                          std::vector<double> &vn_cos,
                          std::vector<double> &vn_sin);
    void calPhoton_SpvnpT_shell();
    void calPhoton_SpvnpT_dTdtau();
    void calPhoton_SpMatrix_dTdtau(double ******dNd2pTdphidydTdtau_eq, 
            double ******dNd2pTdphidydTdtau_tot, double ******dNd2pTdphidydTdtau_diff);
    void outputPhoton_SpvnpT(std::string path, std::string type_str,
                             double ***dNd2pTdphidy,
                             double ***vnypT_cos, double ***vnypT_sin,
                             double **vnpT_cos, double **vnpT_sin,
                             std::vector<double> &vn_cos,
                             std::vector<double> &vn_sin);
    void outputPhoton_SpvnpT_shell(std::string path);
    void outputPhoton_SpvnpT_dTdtau(std::string path);
    void outputPhoton_spectra_dTdtau(std::string path);
    void interpolation2D_bilinear(double varX, std::vector<double> &varY,
                                  double** Table2D_ptr,
                                  std::vector<double> &results);

    struct Table {
        int nx, ny, nz, nw;
        double *x, *y, *z, *w, ****F;
        double x_min, x_max, y_min, y_max, z_min, z_max, w_min, w_max;

        Table() = default;
        Table(std::vector<double> &_x, std::vector<double> &_y, std::vector<double> &_z, std::vector<double> &_w, std::vector<double> &_F);
        double interp(double _x, double _y, double _z, double _w);
        // ~Table();
    };

    double rate(struct Table grid_T, struct Table grid_L,
        double o, double k, double alpha_s, double muB, double T, double m_l);

private:
    Table grid_T;
    Table grid_L;
};
#endif  // SRC_THERMALPHOTON_H_
