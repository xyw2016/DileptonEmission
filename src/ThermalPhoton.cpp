/////////////////////////////////////////////////////////////////////////
//  To do in the future:
//      change the integration routines into gaussian integration
/////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
//#include <omp.h>
#ifndef _OPENMP
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif

#include "Table2D.h"
#include "ThermalPhoton.h"
#include "ParameterReader.h"
#include "gauss_quadrature.h"
#include "Arsenal.h"
#include "data_struct.h"
#include "QGP_NLO.h"

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
using ARSENAL::logarithmic_mass_grid;

using PhysConsts::me;

ThermalPhoton::ThermalPhoton(std::shared_ptr<ParameterReader> paraRdr_in,
                             std::string emissionProcess): grid_T(), grid_L(){

    paraRdr = paraRdr_in;
    emissionProcess_name = emissionProcess;

    neta = paraRdr->getVal("neta");
    nm = paraRdr->getVal("nm");
    np = paraRdr->getVal("np");
    nphi = paraRdr->getVal("nphi");
    nrapidity = paraRdr->getVal("nrapidity");
    norder = paraRdr->getVal("norder");
    rate_path_ = "ph_rates/";

    bRateTable_    = false;
    bShearVisCorr_ = false;
    bBulkVisCorr_  = false;
    bDiffusionCorr_  = false;

    turn_on_muB_ = static_cast<int>(paraRdr->getVal("turn_on_muB", 1));
    include_diff_deltaf = paraRdr->getVal("include_baryondiff_deltaf");
    include_visc_deltaf = paraRdr->getVal("include_shearvisc_deltaf");

    alpha_s = paraRdr->getVal("alpha_s");

    // if muB is off, no need to include diffusion correction
    if(turn_on_muB_==0)
    	include_diff_deltaf = 0;

    //initial variables for photon spectra 
    double p_i = paraRdr->getVal("photon_q_i"); 
    double p_f = paraRdr->getVal("photon_q_f");
    double phi_i = paraRdr->getVal("photon_phi_q_i");
    double phi_f = paraRdr->getVal("photon_phi_q_f");
    double y_i = paraRdr->getVal("photon_y_i");
    double y_f = paraRdr->getVal("photon_y_f");

    double m_i = paraRdr->getVal("dilepton_mass_i");
    double m_f = paraRdr->getVal("dilepton_mass_f");

    // Choose between equal steps and logarithmically spaced mass grid
    use_logarithmic_mass_grid =  paraRdr->getVal("use_logarithmic_mass_grid");

    p = new double [np];
    p_weight = new double [np];
    phi = new double [nphi];
    phi_weight = new double [nphi];

    gauss_quadrature(np, 1, 0.0, 0.0, p_i, p_f, p, p_weight);
    gauss_quadrature(nphi, 1, 0.0, 0.0, phi_i, phi_f, phi, phi_weight);

    // dilepton rapidity
    y_weight = new double [nrapidity];
    
    if (nrapidity > 1) {
        Dy = y_f - y_i;
        dy = Dy/(nrapidity - 1);
        trapezoidal_weights(nrapidity, y_weight);
    } else {
        Dy = 1.0;
        dy = 1.0;
        y_weight[0] = 1.0;
    }

    y.resize(nrapidity, 0);
    
    for (int i=0;i<nrapidity;i++) {
        y[i] = y_i + i*dy;
    }

    // dilepton invariant mass
    if (use_logarithmic_mass_grid) {
        M = logarithmic_mass_grid(m_i, m_f, nm);
    } else {
        M.resize(nm, 0);
        dM = (m_f - m_i)/(nm - 1);
        for (int i=0; i<nm; i++) {
            M[i] = m_i + i*dM;
        }
    }


    // dNd2pTdphidy_eq = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    // dNd2pTdphidy_visc = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    // dNd2pTdphidy_bulkvis = createA3DMatrix(np, nphi, nrapidity, 0.);
    // dNd2pTdphidy_diff = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    // dNd2pTdphidy_tot = createA4DMatrix(nm, np, nphi, nrapidity, 0.);

    // vnypT_cos_eq = createA3DMatrix(norder, np, nrapidity, 0.);
    // vnypT_sin_eq = createA3DMatrix(norder, np, nrapidity, 0.);
    // vnypT_cos_vis = createA3DMatrix(norder, np, nrapidity, 0.);
    // vnypT_sin_vis = createA3DMatrix(norder, np, nrapidity, 0.);
    // vnypT_cos_bulkvis = createA3DMatrix(norder, np, nrapidity, 0.);
    // vnypT_sin_bulkvis = createA3DMatrix(norder, np, nrapidity, 0.);
    // vnypT_cos_tot = createA3DMatrix(norder, np, nrapidity, 0.);
    // vnypT_sin_tot = createA3DMatrix(norder, np, nrapidity, 0.);

    // vnpT_cos_eq = createA2DMatrix(norder, np, 0.);
    // vnpT_sin_eq = createA2DMatrix(norder, np, 0.);
    // vnpT_cos_vis = createA2DMatrix(norder, np, 0.);
    // vnpT_sin_vis = createA2DMatrix(norder, np, 0.);
    // vnpT_cos_bulkvis = createA2DMatrix(norder, np, 0.);
    // vnpT_sin_bulkvis = createA2DMatrix(norder, np, 0.);
    // vnpT_cos_tot = createA2DMatrix(norder, np, 0.);
    // vnpT_sin_tot = createA2DMatrix(norder, np, 0.);

    // vn_cos_eq.resize(norder,0.);
    // vn_sin_eq.resize(norder,0.);
    // vn_cos_vis.resize(norder,0.);
    // vn_sin_vis.resize(norder,0.);
    // vn_cos_bulkvis.resize(norder,0.);
    // vn_sin_bulkvis.resize(norder,0.);
    // vn_cos_tot.resize(norder,0.);
    // vn_sin_tot.resize(norder,0.);

    int differential_flag = paraRdr->getVal("differential_flag");

    if (differential_flag == 1 or differential_flag > 10) {
        nTcut = paraRdr->getVal("nTcut");
        n_tau_cut = paraRdr->getVal("n_tau_cut");

        dNd2pTdphidydTdtau_eq = createA6DMatrix(
                    nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_visc = createA6DMatrix(
                    nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_diff = createA6DMatrix(
                    nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_tot = createA6DMatrix(
                    nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);

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

    delete [] p;
    delete [] p_weight;
    delete [] phi;
    delete [] phi_weight;
    delete [] y_weight;

    // deleteA4DMatrix(dNd2pTdphidy_eq, nm, np, nphi);
    // deleteA4DMatrix(dNd2pTdphidy_visc, nm, np, nphi);
    // deleteA3DMatrix(dNd2pTdphidy_bulkvis, np, nphi);
    // deleteA4DMatrix(dNd2pTdphidy_diff, nm, np, nphi);
    // deleteA4DMatrix(dNd2pTdphidy_tot, nm, np, nphi);

    // deleteA2DMatrix(vnpT_cos_eq, norder);
    // deleteA2DMatrix(vnpT_sin_eq, norder);
    // deleteA2DMatrix(vnpT_cos_vis, norder);
    // deleteA2DMatrix(vnpT_sin_vis, norder);
    // deleteA2DMatrix(vnpT_cos_bulkvis, norder);
    // deleteA2DMatrix(vnpT_sin_bulkvis, norder);
    // deleteA2DMatrix(vnpT_cos_tot, norder);
    // deleteA2DMatrix(vnpT_sin_tot, norder);

    // deleteA3DMatrix(vnypT_cos_eq, norder, np);
    // deleteA3DMatrix(vnypT_sin_eq, norder, np);
    // deleteA3DMatrix(vnypT_cos_vis, norder, np);
    // deleteA3DMatrix(vnypT_sin_vis, norder, np);
    // deleteA3DMatrix(vnypT_cos_bulkvis, norder, np);
    // deleteA3DMatrix(vnypT_sin_bulkvis, norder, np);
    // deleteA3DMatrix(vnypT_cos_tot, norder, np);
    // deleteA3DMatrix(vnypT_sin_tot, norder, np);

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

    //read in equilibrium rate
    ostringstream eqrate_filename_stream;
    eqrate_filename_stream << rate_path_ << "rate_"
                           << emissionProcess_name << "_eqrate.dat";

    string fname = eqrate_filename_stream.str();

    initialize(fname, a_list, B_list, M_list, k_list, rhoT_list, rhoL_list);

    if (a_list.empty() ||B_list.empty() || M_list.empty() || k_list.empty() || rhoT_list.empty()) {
	  // Handle the error here

	    printf("Error: the vectors are still empty. Emission rate table was not read properly.\n");
	  }

  	// Update the tables
    grid_T = Table(a_list,B_list,M_list,k_list,rhoT_list);
    grid_L = Table(a_list,B_list,M_list,k_list,rhoL_list);

	// print some details:
	cout << "-> boundaries of the table ..."  << endl;
	cout <<   " min alpha: " << grid_T.x_min
	     << " , max alpha: " << grid_T.x_max << endl;
	cout <<   " min muB: "   << grid_T.y_min
	     << " , max muB: "   << grid_T.y_max << endl;
	cout <<   " min M: "     << grid_T.z_min
	     << " , max M: "     << grid_T.z_max << endl;
	cout <<   " min k: "     << grid_T.w_min
	     << " , max k: "     << grid_T.w_max << endl;
	cout << " (muB, M and k are given in units of T!)" << endl << endl;

}


void ThermalPhoton::analyticRates(double T, double muB, vector<double> &Eq, 
    double *M_ll, std::vector<double> &eqrate_ptr, int nm,
    int np, int nphi, int nrapidity) {
    for (unsigned int i = 0; i < eqrate_ptr.size(); i++) {
        eqrate_ptr[i] = 1e-16;
    }
}


void ThermalPhoton::analyticRatesShearVis(double T, vector<double> &Eq, 
    double *M_ll, std::vector<double> &visrate_ptr) {
    for (unsigned int i = 0; i < visrate_ptr.size(); i++) {
        visrate_ptr[i] = 0.;
    }
}


void ThermalPhoton::analyticRatesBulkVis(double T, vector<double> &Eq, 
    double *M_ll, std::vector<double> &bulkvis_ptr) {
    for (unsigned int i = 0; i < bulkvis_ptr.size(); i++) {
        bulkvis_ptr[i] = 0.;
    }
}


void ThermalPhoton::FiniteBaryonRates(double T, double muB, double inv_eplusp, double rhoB_over_eplusp, double Eq, 
    double M_ll, double &eqrate_ptr, double &eqrateT_ptr, double &eqrateL_ptr, double &viscrate_ptr, double &diffrate_ptr, 
    int include_visc_deltaf, int include_diff_deltaf) {
    eqrate_ptr = 1.e-16;
    eqrateT_ptr = 1.e-16;
    eqrateL_ptr = 1.e-16;
    viscrate_ptr = 0.;
    diffrate_ptr = 0.;
}


void ThermalPhoton::getPhotonemissionRate(double Eq, double M_ll, double pi_factor, double bulkPi_factor, double diff_factor, 
    double T, double muB, double inv_eplusp, double rhoB_over_eplusp, double &eqrate_ptr, double &eqrateT_ptr, double &eqrateL_ptr, 
	double &viscrate_ptr, double &bulkvis_ptr, double &diffrate_ptr) 
{

	if (bRateTable_) {
        // interpolate NLO equilibrium rate
        double k = sqrt(Eq*Eq-M_ll*M_ll);

  		NLO_rate(grid_T,grid_L,Eq,k,alpha_s,muB,T,me, eqrate_ptr, eqrateT_ptr, eqrateL_ptr);
  		diffrate_ptr = 0.;
    } else {
    	// use LO analytical form
    	FiniteBaryonRates(T, muB, inv_eplusp, rhoB_over_eplusp, Eq, M_ll, eqrate_ptr, eqrateT_ptr, eqrateL_ptr, 
            viscrate_ptr, diffrate_ptr, include_visc_deltaf, include_diff_deltaf);
    }

    viscrate_ptr = pi_factor * viscrate_ptr;
    diffrate_ptr = diff_factor * diffrate_ptr;
}


// this function calculates the spectra for a specific Eq and M at a fluid cell
void ThermalPhoton::calThermalPhotonemission_3d(double (&p_lab_Min)[4], double (&flow_u_mu_Min)[4], double Eq, double M_ll, double pi_zz, double bulkPi, 
	double diff_factor, double T, double muB, double inv_eplusp, double rhoB_over_eplusp, double volume, double fraction,
	double &dNd2pTdphidy_cell_eq, double &dNd2pTdphidy_cell_eqT, double &dNd2pTdphidy_cell_eqL, double &dNd2pTdphidy_cell_visc, 
    double &dNd2pTdphidy_cell_diff, double &dNd2pTdphidy_cell_tot, double &dNd2pTdphidy_cell_lambda_theta) {

    const double volfrac = volume*fraction;

    // Given (T, mu), the emission rate depends on momentum and invariant mass of dilepton
    // photon emission equilibrium rate at local rest cell
    double em_eqrate = 0.;
    // transverse
    double em_eqrateT = 0.;
    // longitudinal
    double em_eqrateL = 0.;
    // photon emission viscous correction at local rest cell
    double em_visrate = 0.;
    // photon emission bulk viscous correction at local rest cell
    double em_bulkvis = 0.;
    // dilepton emission diffusion correction at local rest cell
    double em_diffrate = 0.;
    getPhotonemissionRate(Eq, M_ll, pi_zz, bulkPi, diff_factor, T, muB, inv_eplusp, rhoB_over_eplusp,
                          em_eqrate, em_eqrateT, em_eqrateL, em_visrate, em_bulkvis, em_diffrate);

    double temp_eq_sum   = em_eqrate*volfrac;
    double temp_eqT_sum  = em_eqrateT*volfrac;
    double temp_eqL_sum  = em_eqrateL*volfrac;
    double temp_visc_sum = em_visrate*volfrac;
    // double temp_bulkvis_sum = em_bulkvis*volfrac;
    double temp_diff_sum = em_diffrate*volfrac;

    // spectra
    dNd2pTdphidy_cell_eq   = temp_eq_sum;
    
    //dNd2pTdphidy_cell_eq   = 1.0;

    dNd2pTdphidy_cell_eqT  = temp_eqT_sum;
    dNd2pTdphidy_cell_eqL  = temp_eqL_sum;
    dNd2pTdphidy_cell_visc = temp_eq_sum  + temp_visc_sum;
    dNd2pTdphidy_cell_diff = temp_eq_sum  + temp_diff_sum;
    dNd2pTdphidy_cell_tot  = (temp_eq_sum + temp_visc_sum + temp_diff_sum);


    // dilepton polarization
    double p_lab_Min_norm[3]; 
    double p_lab_min_vec_square = sqrt( p_lab_Min[1]*p_lab_Min[1]
                                        + p_lab_Min[2]*p_lab_Min[2]
                                        + p_lab_Min[3]*p_lab_Min[3]);

    p_lab_Min_norm[0]=p_lab_Min[1]/p_lab_min_vec_square;
    p_lab_Min_norm[1]=p_lab_Min[2]/p_lab_min_vec_square;
    p_lab_Min_norm[2]=p_lab_Min[3]/p_lab_min_vec_square;

   
    double u_dot_kvec_norm = flow_u_mu_Min[1]*p_lab_Min_norm[0]
                            +flow_u_mu_Min[2]*p_lab_Min_norm[1]
                            +flow_u_mu_Min[3]*p_lab_Min_norm[2];

    

    double u_dot_kvec = flow_u_mu_Min[1]*p_lab_Min[1]
                        +flow_u_mu_Min[2]*p_lab_Min[2]
                        +flow_u_mu_Min[3]*p_lab_Min[3];
    
    double factor_pol_0 = pow((p_lab_Min[0]*u_dot_kvec_norm - p_lab_min_vec_square*flow_u_mu_Min[0]),2)/
                          (pow((p_lab_Min[0]*flow_u_mu_Min[0]-u_dot_kvec),2)-M_ll*M_ll);
    factor_pol_0 = factor_pol_0 - 1./3.0;

    double factor_pol_1 = me*me/(p_lab_min_vec_square*p_lab_min_vec_square);
    factor_pol_1 = 0.0;
    double rho_delta= temp_eqT_sum - temp_eqL_sum;
    double rho_V = temp_eqT_sum*2 + temp_eqL_sum;


    double lambda_theta1 = factor_pol_0*(1.0 - 4.0*factor_pol_1)*rho_delta;
    double lambda_theta2 = ( 4.0*(1+2.0*factor_pol_1)*rho_V/3.0 - lambda_theta1 );
    

    double lambda_theta = 3.0*lambda_theta1/lambda_theta2;

    if ( fabs(lambda_theta2) < 1e-18){
        lambda_theta = 0.0;
    }


    dNd2pTdphidy_cell_lambda_theta = lambda_theta*dNd2pTdphidy_cell_eq;
    //dNd2pTdphidy_cell_lambda_theta = lambda_theta;
    
    if(std::isnan(dNd2pTdphidy_cell_lambda_theta))
    { 
        std::cout<< rho_delta << " "<< rho_V <<" "<<lambda_theta1<<" "<<factor_pol_0<<" "<<lambda_theta2<< " "<<dNd2pTdphidy_cell_lambda_theta<<std::endl;
    }


    

}


// functions to get distributions in T and tau; contributions from all channels are included
void ThermalPhoton::calPhoton_SpMatrix_dTdtau(double ******dNd2pTdphidydTdtau_eq_temp, double ******dNd2pTdphidydTdtau_visc_temp, 
    double ******dNd2pTdphidydTdtau_diff_temp, double ******dNd2pTdphidydTdtau_tot_temp) {
    for (int i = 0; i < nTcut; i++) {
        for (int j = 0; j < n_tau_cut; j++) {
		    for (int k = 0; k < nm; k++) {
		        for (int l = 0; l < np; l++) {
		            for (int m = 0; m < nphi; m++) {
		                for (int n = 0; n < nrapidity; n++) {
		                    dNd2pTdphidydTdtau_eq[i][j][k][l][m][n] = dNd2pTdphidydTdtau_eq_temp[i][j][k][l][m][n];
                            dNd2pTdphidydTdtau_visc[i][j][k][l][m][n] = dNd2pTdphidydTdtau_visc_temp[i][j][k][l][m][n];
		                    dNd2pTdphidydTdtau_diff[i][j][k][l][m][n] = dNd2pTdphidydTdtau_diff_temp[i][j][k][l][m][n];
		                    dNd2pTdphidydTdtau_tot[i][j][k][l][m][n] = dNd2pTdphidydTdtau_tot_temp[i][j][k][l][m][n];
		                }
		            }
		        }
		    }
		}
	}
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
                            dNpTdpTdydTdtau_eq[i][j][m][k] += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]
                                    *phi_weight[l]*y_weight[irap]*dy/Dy);
                            dNpTdpTdydTdtau_visc[i][j][m][k] += (
                                    dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap]
                                    *phi_weight[l]*y_weight[irap]*dy/Dy);
                            dNpTdpTdydTdtau_diff[i][j][m][k] += (
                                    dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap]
                                    *phi_weight[l]*y_weight[irap]*dy/Dy);
                            dNpTdpTdydTdtau_tot[i][j][m][k] += (
                                    dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap]
                                    *phi_weight[l]*y_weight[irap]*dy/Dy);
                        }
                    }
                    // 
                    dNpTdpTdydTdtau_eq[i][j][m][k] = dNpTdpTdydTdtau_eq[i][j][m][k]/(2*M_PI);
                    dNpTdpTdydTdtau_visc[i][j][m][k] = dNpTdpTdydTdtau_visc[i][j][m][k]/(2*M_PI);
                    dNpTdpTdydTdtau_diff[i][j][m][k] = dNpTdpTdydTdtau_diff[i][j][m][k]/(2*M_PI);
                    dNpTdpTdydTdtau_tot[i][j][m][k] = dNpTdpTdydTdtau_tot[i][j][m][k]/(2*M_PI);
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


void ThermalPhoton::outputPhoton_Spectra_dTdtau(string path, double Tcut_high, double Tcut_low, double tau_cut_high, double tau_cut_low) {
    // calculate the inverse slope of the photon spectra at T-tau interval
    // integrated out phi and rapidity
    // pT differential spectra, dN/(2pi pTdpT MdM dy)

    double dT = (Tcut_high - Tcut_low)/(nTcut - 1);
    double dtau = (tau_cut_high - tau_cut_low)/(n_tau_cut - 1);

    ostringstream filename_SpdTdtau_eq;
    ostringstream filename_SpdTdtau_visc;
    ostringstream filename_SpdTdtau_diff;
    ostringstream filename_SpdTdtau_tot;

    filename_SpdTdtau_eq 	<< path << emissionProcess_name
                         	<< "_SpdTdtau_eq.dat";
    filename_SpdTdtau_visc  << path << emissionProcess_name
                            << "_SpdTdtau_visc.dat";
    filename_SpdTdtau_diff 	<< path << emissionProcess_name
                          	<< "_SpdTdtau_diff.dat";
    filename_SpdTdtau_tot 	<< path << emissionProcess_name
                          	<< "_SpdTdtau_tot.dat";

    ofstream ofeq(filename_SpdTdtau_eq.str().c_str());
    ofstream ofvisc(filename_SpdTdtau_visc.str().c_str());
    ofstream ofdiff(filename_SpdTdtau_diff.str().c_str());
    ofstream oftot(filename_SpdTdtau_tot.str().c_str());

    for (int i = 0; i < nTcut; i++) {
        double T_local = Tcut_low + i*dT;
        for (int j = 0; j < n_tau_cut; j++) {
            double tau_local = tau_cut_low + j*dtau;

            for (int m = 0; m < nm; m++) {
            	ofeq 	<< scientific << setw(18) << setprecision(8) 
                 		<< T_local << "   "  << tau_local << "   ";
	            ofdiff 	<< scientific << setw(18) << setprecision(8) 
	                   	<< T_local << "   "  << tau_local << "   ";
	            oftot 	<< scientific << setw(18) << setprecision(8) 
	                  	<< T_local << "   "  << tau_local << "   ";
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
    filename_SpdT_eq    << path << emissionProcess_name
                            << "_SpdT_eq.dat";
    filename_SpdT_visc  << path << emissionProcess_name
                            << "_SpdT_visc.dat";
    filename_SpdT_diff  << path << emissionProcess_name
                            << "_SpdT_diff.dat";
    filename_SpdT_tot   << path << emissionProcess_name
                            << "_SpdT_tot.dat";

    ofstream ofdTeq(filename_SpdT_eq.str().c_str());
    ofstream ofdTvisc(filename_SpdT_visc.str().c_str());
    ofstream ofdTdiff(filename_SpdT_diff.str().c_str());
    ofstream ofdTtot(filename_SpdT_tot.str().c_str());

    for (int i = 0; i < nTcut; i++) {
        double T_local = Tcut_low + i*dT;

        for (int m = 0; m < nm; m++) {
            ofdTeq    << scientific << setw(18) << setprecision(8) 
                    << T_local << "   ";
            ofdTvisc  << scientific << setw(18) << setprecision(8) 
                    << T_local << "   ";
            ofdTdiff  << scientific << setw(18) << setprecision(8) 
                    << T_local << "   ";
            ofdTtot   << scientific << setw(18) << setprecision(8) 
                    << T_local << "   ";
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
    filename_Spdtau_eq    << path << emissionProcess_name
                            << "_Spdtau_eq.dat";
    filename_Spdtau_visc  << path << emissionProcess_name
                            << "_Spdtau_visc.dat";
    filename_Spdtau_diff  << path << emissionProcess_name
                            << "_Spdtau_diff.dat";
    filename_Spdtau_tot   << path << emissionProcess_name
                            << "_Spdtau_tot.dat";

    ofstream ofdtaueq(filename_Spdtau_eq.str().c_str());
    ofstream ofdtauvisc(filename_Spdtau_visc.str().c_str());
    ofstream ofdtaudiff(filename_Spdtau_diff.str().c_str());
    ofstream ofdtautot(filename_Spdtau_tot.str().c_str());

    for (int j = 0; j < n_tau_cut; j++) {
        double tau_local = tau_cut_low + j*dtau;

        for (int m = 0; m < nm; m++) {
            ofdtaueq    << scientific << setw(18) << setprecision(8) 
                    << tau_local << "   ";
            ofdtauvisc  << scientific << setw(18) << setprecision(8) 
                    << tau_local << "   ";
            ofdtaudiff  << scientific << setw(18) << setprecision(8) 
                    << tau_local << "   ";
            ofdtautot   << scientific << setw(18) << setprecision(8) 
                    << tau_local << "   ";
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
                        double weight = p[k]*p_weight[k]*phi_weight[l]; // pT and phi_p integrated out
                        for (int irap = 0; irap < nrapidity; irap++) {

                            weight *=y_weight[irap]*dy/Dy;

                            dNdydTdtau_eq[i][j][m] += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]*weight);
                            dNdydTdtau_visc[i][j][m] += (
                                    dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap]*weight);
                            dNdydTdtau_diff[i][j][m] += (
                                    dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap]*weight);
                            dNdydTdtau_tot[i][j][m] += (
                                    dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap]*weight);

                            for (int order = 0; order < norder; order++) {
                                vndTdtau_cos_eq[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]
                                    *weight*cos(order*phi[l]));
                                vndTdtau_sin_eq[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]
                                    *weight*sin(order*phi[l]));
                                vndTdtau_cos_visc[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap]
                                    *weight*cos(order*phi[l]));
                                vndTdtau_sin_visc[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_visc[i][j][m][k][l][irap]
                                    *weight*sin(order*phi[l]));
                                vndTdtau_cos_diff[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap]
                                    *weight*cos(order*phi[l]));
                                vndTdtau_sin_diff[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap]
                                    *weight*sin(order*phi[l]));
                                vndTdtau_cos_tot[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap]
                                    *weight*cos(order*phi[l]));
                                vndTdtau_sin_tot[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap]
                                    *weight*sin(order*phi[l]));
                            }
                        }
                    }
                }
            }

            for (int m = 0; m < nm; m++) {
                for (int order = 1; order < norder ; order++) {
                    vndTdtau_cos_eq[i][j][m][order] = (
                            vndTdtau_cos_eq[i][j][m][order]
                            /(dNdydTdtau_eq[i][j][m] + eps));
                    vndTdtau_sin_eq[i][j][m][order] = (
                            vndTdtau_sin_eq[i][j][m][order]
                            /(dNdydTdtau_eq[i][j][m] + eps));
                    vndTdtau_cos_visc[i][j][m][order] = (
                            vndTdtau_cos_visc[i][j][m][order]
                            /(dNdydTdtau_visc[i][j][m] + eps));
                    vndTdtau_sin_visc[i][j][m][order] = (
                            vndTdtau_sin_visc[i][j][m][order]
                            /(dNdydTdtau_visc[i][j][m] + eps));
                    vndTdtau_cos_diff[i][j][m][order] = (
                            vndTdtau_cos_diff[i][j][m][order]
                            /(dNdydTdtau_diff[i][j][m] + eps));
                    vndTdtau_sin_diff[i][j][m][order] = (
                            vndTdtau_sin_diff[i][j][m][order]
                            /(dNdydTdtau_diff[i][j][m] + eps));
                    vndTdtau_cos_tot[i][j][m][order] = (
                            vndTdtau_cos_tot[i][j][m][order]
                            /(dNdydTdtau_tot[i][j][m] + eps));
                    vndTdtau_sin_tot[i][j][m][order] = (
                            vndTdtau_sin_tot[i][j][m][order]
                            /(dNdydTdtau_tot[i][j][m] + eps));
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


void ThermalPhoton::outputPhoton_Spvn_dTdtau(string path, double Tcut_high, double Tcut_low, double tau_cut_high, double tau_cut_low) {

    double dT = (Tcut_high - Tcut_low)/(nTcut - 1);
    double dtau = (tau_cut_high - tau_cut_low)/(n_tau_cut - 1);

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

    for(int i = 0; i < nTcut; i++)
    {
       double T_local = Tcut_low + i*dT;
       for(int j = 0; j < n_tau_cut; j++)
       {
          	double tau_local = tau_cut_low + j*dtau;
          	fphotondNdy_eq << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
            fphotondNdy_visc << scientific << setw(18) << setprecision(8) 
                    << T_local << "   "  << tau_local << "   ";
          	fphotondNdy_diff << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
          	fphotondNdy_tot << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
          	for (int m = 0; m < nm; m++) {
          		double M_ll = M[m];
              	fphotondNdy_eq << dNdydTdtau_eq[i][j][m]*M_ll << "    ";
              	fphotondNdy_visc << dNdydTdtau_visc[i][j][m]*M_ll << "    ";
              	fphotondNdy_diff << dNdydTdtau_diff[i][j][m]*M_ll << "    ";
              	fphotondNdy_tot << dNdydTdtau_tot[i][j][m]*M_ll << "    ";
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

    filename_stream_dNdydT_eq << path << emissionProcess_name
                                   << "_dNdydT_eq.dat";
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

    for(int i = 0; i < nTcut; i++)
    {
        double T_local = Tcut_low + i*dT;

        fphotondNdydT_eq << scientific << setw(18) << setprecision(8) 
                << T_local << "   ";
        fphotondNdydT_visc << scientific << setw(18) << setprecision(8) 
                << T_local << "   ";
        fphotondNdydT_diff << scientific << setw(18) << setprecision(8) 
                << T_local << "   ";
        fphotondNdydT_tot << scientific << setw(18) << setprecision(8) 
                << T_local << "   ";
        for (int m = 0; m < nm; m++) {
            double M_ll = M[m];
            fphotondNdydT_eq << dNdydT_eq[i][m]*M_ll << "    ";
            fphotondNdydT_visc << dNdydT_visc[i][m]*M_ll << "    ";
            fphotondNdydT_diff << dNdydT_diff[i][m]*M_ll << "    ";
            fphotondNdydT_tot << dNdydT_tot[i][m]*M_ll << "    ";
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

    for(int i = 0; i < n_tau_cut; i++)
    {
        double tau_local = tau_cut_low + i*dtau;

        fphotondNdydtau_eq << scientific << setw(18) << setprecision(8) 
                << tau_local << "   ";
        fphotondNdydtau_visc << scientific << setw(18) << setprecision(8) 
                << tau_local << "   ";
        fphotondNdydtau_diff << scientific << setw(18) << setprecision(8) 
                << tau_local << "   ";
        fphotondNdydtau_tot << scientific << setw(18) << setprecision(8) 
                << tau_local << "   ";
        for (int m = 0; m < nm; m++) {
            double M_ll = M[m];
            fphotondNdydtau_eq << dNdydtau_eq[i][m]*M_ll << "    ";
            fphotondNdydtau_visc << dNdydtau_visc[i][m]*M_ll << "    ";
            fphotondNdydtau_diff << dNdydtau_diff[i][m]*M_ll << "    ";
            fphotondNdydtau_tot << dNdydtau_tot[i][m]*M_ll << "    ";
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

    filename_stream_vndTdtau_eq  << path << emissionProcess_name
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

    for(int i = 0; i < nTcut; i++)
    {
        double T_local = Tcut_low + i*dT;
        for(int j = 0; j < n_tau_cut; j++)
        {
            double tau_local = tau_cut_low + j*dtau;
            fphotonvn_eq << scientific << setw(18) << setprecision(8) 
                << T_local << "   "  << tau_local << "   ";
            fphotonvn_visc << scientific << setw(18) << setprecision(8) 
                << T_local << "   "  << tau_local << "   ";
            fphotonvn_diff << scientific << setw(18) << setprecision(8) 
                << T_local << "   "  << tau_local << "   ";
            fphotonvn_tot << scientific << setw(18) << setprecision(8) 
                << T_local << "   "  << tau_local << "   ";

            for (int m = 0; m < nm; m++) {
                for(int order = 1; order < norder; order++)
                    {
                        fphotonvn_eq << order << "   " 
                            << vndTdtau_cos_eq[i][j][m][order] << "    " << vndTdtau_sin_eq[i][j][m][order] << "    "
                            << sqrt(pow(vndTdtau_cos_eq[i][j][m][order], 2) + pow(vndTdtau_sin_eq[i][j][m][order], 2)) << "  ";
                        fphotonvn_visc << order << "   " 
                            << vndTdtau_cos_visc[i][j][m][order] << "    " << vndTdtau_sin_visc[i][j][m][order] << "    "
                            << sqrt(pow(vndTdtau_cos_visc[i][j][m][order], 2) + pow(vndTdtau_sin_visc[i][j][m][order], 2)) << "  ";
                        fphotonvn_diff << order << "   " 
                            << vndTdtau_cos_diff[i][j][m][order] << "    " << vndTdtau_sin_diff[i][j][m][order] << "    "
                            << sqrt(pow(vndTdtau_cos_diff[i][j][m][order], 2) + pow(vndTdtau_sin_diff[i][j][m][order], 2)) << "  ";
                        fphotonvn_tot << order << "   " 
                            << vndTdtau_cos_tot[i][j][m][order] << "    " << vndTdtau_sin_tot[i][j][m][order] << "    "
                            << sqrt(pow(vndTdtau_cos_tot[i][j][m][order], 2) + pow(vndTdtau_sin_tot[i][j][m][order], 2)) << "  ";
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


// pT-integrated and pT-differential spectra and flows for individual channels

// void ThermalPhoton::calPhoton_SpvnpT(
//         double ***dNd2pTdphidy, double ***vnypT_cos, double ***vnypT_sin,
//         double **vnpT_cos, double **vnpT_sin,
//         vector<double> &vn_cos, vector<double> &vn_sin) {
//     // calculate the photon spectra and differential vn
//     const double eps = 1e-15;
//     for (int i = 0; i < np; i++) {
//         for (int k = 0; k < nrapidity; k++) {
//             for (int j = 0; j < nphi; j++) {
//                 double weight = phi_weight[j];
//                 for (int order = 0; order < norder; order++) {
//                     double cos_tmp = (dNd2pTdphidy[i][j][k]
//                                       *cos(order*phi[j])*weight);
//                     double sin_tmp = (dNd2pTdphidy[i][j][k]
//                                       *sin(order*phi[j])*weight);
//                     vnypT_cos[order][i][k] += cos_tmp;
//                     vnypT_sin[order][i][k] += sin_tmp;
//                     if (std::abs(y[k]) < 0.5) {
//                         // only integrate mid-rapidity
//                         vnpT_cos[order][i] += cos_tmp*dy;
//                         vnpT_sin[order][i] += sin_tmp*dy;
//                     }
//                 }
//             }
//         }
//         double p_weight_factor = p[i]*p_weight[i];
//         for (int order = 0; order < norder ; order++) {
//             vn_cos[order] += vnpT_cos[order][i]*p_weight_factor;
//             vn_sin[order] += vnpT_sin[order][i]*p_weight_factor;

//             // vn(pT)
//             if (order > 0) {
//                 vnpT_cos[order][i] /= (vnpT_cos[0][i] + eps);
//                 vnpT_sin[order][i] /= (vnpT_cos[0][i] + eps);
//                 for (int k = 0; k < nrapidity; k++) {
//                     vnypT_cos[order][i][k] /= (vnypT_cos[0][i][k] + eps);
//                     vnypT_sin[order][i][k] /= (vnypT_cos[0][i][k] + eps);
//                 }
//             }
//         }
//         vnpT_cos[0][i] /= (2*M_PI);  // dN/(2pi dy pT dpT)
//     }
//     for (int order = 1; order < norder ; order++) {
//         // vn
//         vn_cos[order] /= (vn_cos[0] + eps);
//         vn_sin[order] /= (vn_cos[0] + eps);
//     }
// }


// void ThermalPhoton::outputPhoton_SpvnpT(string path, string type_str,
//         double ***dNd2pTdphidy, double ***vnypT_cos, double ***vnypT_sin,
//         double **vnpT_cos, double **vnpT_sin,
//         vector<double> &vn_cos, vector<double> &vn_sin) {
//     ostringstream filename_stream_SpMatrix;
//     ostringstream filename_stream_Spvn;
//     ostringstream filename_stream_inte_Spvn;

//     filename_stream_SpMatrix << path << emissionProcess_name
//                              << "_Spvn_" << type_str << "_ypTdiff.dat";
//     filename_stream_Spvn << path << emissionProcess_name << "_Spvn_"
//                          << type_str << ".dat";
//     filename_stream_inte_Spvn << path << emissionProcess_name
//                               << "_Spvn_" << type_str << "_inte.dat";

//     ofstream fphotonSpMatrix(filename_stream_SpMatrix.str().c_str());
//     ofstream fphotonSpvn(filename_stream_Spvn.str().c_str());
//     ofstream fphotoninteSpvn(filename_stream_inte_Spvn.str().c_str());

//     for (int k = 0; k < nrapidity; k++) {
//         for (int j = 0; j < np; j++) {
//             fphotonSpMatrix << scientific << setprecision(6) << setw(16)
//                             << y[k] << "  " << p[j] << "  "
//                             << vnypT_cos[0][j][k] << "  ";
//             for (int order = 1; order < norder; order++) {
//                 fphotonSpMatrix << scientific << setprecision(6) << setw(16)
//                                 << vnypT_cos[order][j][k] << "  "
//                                 << vnypT_sin[order][j][k] << "  ";
//             }
//             fphotonSpMatrix << endl;
//         }
//     }

//     for (int i = 0; i < np; i++) {
//         fphotonSpvn << scientific << setprecision(6) << setw(16)
//                     << p[i] << "  " << vnpT_cos[0][i] << "  " ;
//         for (int order = 1; order < norder; order++) {
//             fphotonSpvn << scientific << setprecision(6) << setw(16)
//                         << vnpT_cos[order][i] << "  "
//                         << vnpT_sin[order][i] << "  ";
//         }
//         fphotonSpvn << endl;
//     }

//     for (int order = 0; order < norder; order++) {
//         fphotoninteSpvn << scientific << setprecision(6) << setw(16)
//                         << order << "   " << vn_cos[order] << "   "
//                         << vn_sin[order] << endl;
//     }

//     fphotonSpMatrix.close();
//     fphotonSpvn.close();
//     fphotoninteSpvn.close();
// }


// void ThermalPhoton::calPhoton_SpvnpT_shell() {
//     calPhoton_SpvnpT(dNd2pTdphidy_eq, vnypT_cos_eq, vnypT_sin_eq,
//                      vnpT_cos_eq, vnpT_sin_eq,
//                      vn_cos_eq, vn_sin_eq);
//     calPhoton_SpvnpT(dNd2pTdphidy_vis, vnypT_cos_vis, vnypT_sin_vis,
//                      vnpT_cos_vis, vnpT_sin_vis,
//                      vn_cos_vis, vn_sin_vis);
//     calPhoton_SpvnpT(dNd2pTdphidy_bulkvis,
//                      vnypT_cos_bulkvis, vnypT_sin_bulkvis,
//                      vnpT_cos_bulkvis, vnpT_sin_bulkvis,
//                      vn_cos_bulkvis, vn_sin_bulkvis);
//     calPhoton_SpvnpT(dNd2pTdphidy_tot, vnypT_cos_tot, vnypT_sin_tot,
//                      vnpT_cos_tot, vnpT_sin_tot,
//                      vn_cos_tot, vn_sin_tot);
// }


// void ThermalPhoton::outputPhoton_SpvnpT_shell(string path) {
//     outputPhoton_SpvnpT(path, "eq", dNd2pTdphidy_eq,
//                         vnypT_cos_eq, vnypT_sin_eq,
//                         vnpT_cos_eq, vnpT_sin_eq,
//                         vn_cos_eq, vn_sin_eq);
//     outputPhoton_SpvnpT(path, "vis", dNd2pTdphidy_vis,
//                         vnypT_cos_vis, vnypT_sin_vis,
//                         vnpT_cos_vis, vnpT_sin_vis,
//                         vn_cos_vis, vn_sin_vis);
//     outputPhoton_SpvnpT(path, "bulkvis", dNd2pTdphidy_bulkvis,
//                         vnypT_cos_bulkvis, vnypT_sin_bulkvis,
//                         vnpT_cos_bulkvis, vnpT_sin_bulkvis,
//                         vn_cos_bulkvis, vn_sin_bulkvis);
//     outputPhoton_SpvnpT(path, "tot", dNd2pTdphidy_tot,
//                         vnypT_cos_tot, vnypT_sin_tot,
//                         vnpT_cos_tot, vnpT_sin_tot,
//                         vn_cos_tot, vn_sin_tot);
// }


