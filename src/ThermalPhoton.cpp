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
#include <omp.h>

#include "Table2D.h"
#include "ThermalPhoton.h"
#include "ParameterReader.h"
#include "gauss_quadrature.h"
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

ThermalPhoton::ThermalPhoton(std::shared_ptr<ParameterReader> paraRdr_in,
                             std::string emissionProcess) {
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

    p = new double [np];
    p_weight = new double [np];
    phi = new double [nphi];
    phi_weight = new double [nphi];

    gauss_quadrature(np, 1, 0.0, 0.0, p_i, p_f, p, p_weight);
    gauss_quadrature(nphi, 1, 0.0, 0.0, phi_i, phi_f, phi, phi_weight);

    if (nrapidity > 1) {
        dy = (y_f - y_i)/(nrapidity - 1 + 1e-100);
    } else {
        dy = 1.0;
    }

    y.resize(nrapidity, 0);
    theta.resize(nrapidity, 0);
    for (int i=0;i<nrapidity;i++) {
        y[i] = y_i + i*dy;
        theta[i] = acos(tanh(y[i]));  //rapidity's corresponding polar angle
    }

    M.resize(nm, 0);
    dM = (m_f - m_i)/(nm + 1e-100);
    for (int i=0;i<nm;i++) {
        M[i] = m_i + i*dM;
    }


    dNd2pTdphidy_eq = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    // dNd2pTdphidy_vis = createA3DMatrix(np, nphi, nrapidity, 0.);
    // dNd2pTdphidy_bulkvis = createA3DMatrix(np, nphi, nrapidity, 0.);
    dNd2pTdphidy_diff = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidy_tot = createA4DMatrix(nm, np, nphi, nrapidity, 0.);

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

    int diff_flag = paraRdr->getVal("differential_flag");

    if (diff_flag == 1 or diff_flag > 10) {
        nTcut = paraRdr->getVal("nTcut");
        n_tau_cut = paraRdr->getVal("n_tau_cut");

        Tcut_high = paraRdr->getVal("T_cuthigh");
        Tcut_low = paraRdr->getVal("T_cutlow");
        tau_cut_high = paraRdr->getVal("tau_end");
        tau_cut_low = paraRdr->getVal("tau_start");

        dNd2pTdphidydTdtau_eq = createA6DMatrix(
                    nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        // dNd2pTdphidydTdtau_vis = createA5DMatrix(
        //             nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_diff = createA6DMatrix(
                    nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_tot = createA6DMatrix(
                    nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);

        dNdydTdtau_eq = createA3DMatrix(nTcut, n_tau_cut, nm, 0.);
        // dNdydTdtau_vis = createA2DMatrix(nTcut, n_tau_cut, nm, 0.);
        dNdydTdtau_diff = createA3DMatrix(nTcut, n_tau_cut, nm, 0.);
        dNdydTdtau_tot = createA3DMatrix(nTcut, n_tau_cut, nm, 0.);

        vndTdtau_cos_eq = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
        vndTdtau_sin_eq = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
        // vndTdtau_cos_vis = createA3DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
        // vndTdtau_sin_vis = createA3DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
        vndTdtau_cos_diff = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
        vndTdtau_sin_diff = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
        vndTdtau_cos_tot = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
        vndTdtau_sin_tot = createA4DMatrix(nTcut, n_tau_cut, nm, norder, 0.);
    }
}


ThermalPhoton::~ThermalPhoton() {
    if (bRateTable_) {
        int TbsizeX = EmissionrateTb_sizeX;
        deleteA2DMatrix(Emission_eqrateTb_ptr, TbsizeX);
        if (bShearVisCorr_) {
            deleteA2DMatrix(Emission_viscous_rateTb_ptr, TbsizeX);
        }
        if (bBulkVisCorr_) {
            deleteA2DMatrix(Emission_bulkvis_rateTb_ptr, TbsizeX);
        }
    }

    delete [] p;
    delete [] p_weight;
    delete [] phi;
    delete [] phi_weight;

    deleteA4DMatrix(dNd2pTdphidy_eq, nm, np, nphi);
    // deleteA3DMatrix(dNd2pTdphidy_vis, np, nphi);
    // deleteA3DMatrix(dNd2pTdphidy_bulkvis, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_diff, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_tot, nm, np, nphi);

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

    int diff_flag = paraRdr->getVal("differential_flag");
    if (diff_flag == 1 or diff_flag > 10) {
        deleteA3DMatrix(dNdydTdtau_eq, nTcut, nm);
        // deleteA3DMatrix(dNdydTdtau_vis, nTcut, nm);
        deleteA3DMatrix(dNdydTdtau_diff, nTcut, nm);
        deleteA3DMatrix(dNdydTdtau_tot, nTcut, nm);

        deleteA4DMatrix(vndTdtau_cos_eq, nTcut, n_tau_cut, nm);
        deleteA4DMatrix(vndTdtau_sin_eq, nTcut, n_tau_cut, nm);
        // deleteA4DMatrix(vndTdtau_cos_vis, nTcut, n_tau_cut, nm);
        // deleteA4DMatrix(vndTdtau_sin_vis, nTcut, n_tau_cut, nm);
        deleteA4DMatrix(vndTdtau_cos_diff, nTcut, n_tau_cut, nm);
        deleteA4DMatrix(vndTdtau_sin_diff, nTcut, n_tau_cut, nm);
        deleteA4DMatrix(vndTdtau_cos_tot, nTcut, n_tau_cut, nm);
        deleteA4DMatrix(vndTdtau_sin_tot, nTcut, n_tau_cut, nm);

        deleteA6DMatrix(dNd2pTdphidydTdtau_eq, nTcut, n_tau_cut, nm, np, nphi);
        // deleteA6DMatrix(dNd2pTdphidydTdtau_vis, nTcut, n_tau_cut, nm, np, nphi);
        deleteA6DMatrix(dNd2pTdphidydTdtau_diff, nTcut, n_tau_cut, nm, np, nphi);
        deleteA6DMatrix(dNd2pTdphidydTdtau_tot, nTcut, n_tau_cut, nm, np, nphi);
    }
}


void ThermalPhoton::setupEmissionrateFromFile(
        double Xmin, double dX, double Ymin, double dY,
        bool bShearVisCorr, bool bBulkVisCorr, bool bDiffusionCorr) {
    bRateTable_    = true;
    bShearVisCorr_ = bShearVisCorr;
    bBulkVisCorr_  = bBulkVisCorr;
    bDiffusionCorr_  = bDiffusionCorr;
    EmissionrateTb_Xmin = Xmin;
    EmissionrateTb_Ymin = Ymin;
    EmissionrateTb_dX = dX;
    EmissionrateTb_dY = dY;
    readEmissionrate(emissionProcess_name);
}


void ThermalPhoton::readEmissionrate(string emissionProcess) {

    //equilibrium rate
    ostringstream eqrate_filename_stream;
    eqrate_filename_stream << rate_path_ << "rate_"
                           << emissionProcess << "_eqrate.dat";
    Photonemission_eqrateTable_ptr = std::unique_ptr<Table2D>(
                        new Table2D(eqrate_filename_stream.str().c_str()));
    EmissionrateTb_sizeX = Photonemission_eqrateTable_ptr->getTbsizeX();
    EmissionrateTb_sizeY = Photonemission_eqrateTable_ptr->getTbsizeY();
    Emission_eqrateTb_ptr = createA2DMatrix(EmissionrateTb_sizeX,
                                            EmissionrateTb_sizeY, 0);

    // shear correction
    if (bShearVisCorr_) {
        ostringstream visrate_filename_stream;
        visrate_filename_stream << rate_path_ << "rate_"
                                << emissionProcess << "_viscous.dat";
        Photonemission_viscous_rateTable_ptr = std::unique_ptr<Table2D>(
                    new Table2D(visrate_filename_stream.str().c_str()));
        Emission_viscous_rateTb_ptr = createA2DMatrix(EmissionrateTb_sizeX,
                                                      EmissionrateTb_sizeY, 0);
    }

    // bulk correction
    if (bBulkVisCorr_) {
        ostringstream bulkvisrate_filename_stream;
        bulkvisrate_filename_stream << rate_path_ << "rate_"
                                    << emissionProcess << "_bulkvis.dat";
        Photonemission_bulkvis_rateTable_ptr = std::unique_ptr<Table2D>(
                    new Table2D(bulkvisrate_filename_stream.str().c_str()));
        Emission_bulkvis_rateTb_ptr = createA2DMatrix(EmissionrateTb_sizeX,
                                                      EmissionrateTb_sizeY, 0);
    }

    EmissionrateTb_Yidxptr.resize(EmissionrateTb_sizeY, 0);
    for (int i = 0; i < EmissionrateTb_sizeY; i++)
        EmissionrateTb_Yidxptr[i] = EmissionrateTb_Ymin + i*EmissionrateTb_dY;

    // take log for the equilibrium emission rates for better
    // interpolation precisions;
    // take ratio of viscous rates w.r.t eq. rates for faster performance
    for (int i = 0; i < EmissionrateTb_sizeX; i++) {
        for (int j = 0; j < EmissionrateTb_sizeY; j++) {
            Emission_eqrateTb_ptr[i][j] =
                log(Photonemission_eqrateTable_ptr->getTbdata(i,j) + 1e-30);
            if (bShearVisCorr_) {
                Emission_viscous_rateTb_ptr[i][j] = (
                    Photonemission_viscous_rateTable_ptr->getTbdata(i,j)
                    /(Photonemission_eqrateTable_ptr->getTbdata(i,j) + 1e-30));
            }
            if (bBulkVisCorr_) {
                Emission_bulkvis_rateTb_ptr[i][j] = (
                    Photonemission_bulkvis_rateTable_ptr->getTbdata(i,j)
                    /(Photonemission_eqrateTable_ptr->getTbdata(i,j) + 1e-30));
            }
        }
    }
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


void ThermalPhoton::FiniteBaryonRates(double T, double muB, double rhoB_over_eplusp, double Eq, 
    double M_ll, double &eqrate_ptr, double &diffrate_ptr, int include_diff_deltaf) {
    eqrate_ptr = 1.e-16;
    diffrate_ptr = 0.;
}


void ThermalPhoton::getPhotonemissionRate(double Eq, double M_ll, double pi_zz, double bulkPi,
	double diff_factor, double T, double muB, double rhoB_over_eplusp, double &eqrate_ptr, 
	double &visrate_ptr, double &bulkvis_ptr, double &diffrate_ptr) 
{

    FiniteBaryonRates(T, muB, rhoB_over_eplusp, Eq, M_ll, eqrate_ptr, diffrate_ptr, include_diff_deltaf);

    diffrate_ptr = diff_factor * diffrate_ptr;
}


// this function calculates the spectra for a specific Eq and M at a fluid cell
void ThermalPhoton::calThermalPhotonemission_3d(double Eq, double M_ll, double pi_zz, double bulkPi, 
	double diff_factor, double T, double muB, double rhoB_over_eplusp, double volume, double fraction,
	double &dNd2pTdphidy_cell_eq, double &dNd2pTdphidy_cell_diff, double &dNd2pTdphidy_cell_tot) {

    const double volfrac = volume*fraction;

    // Given (T, mu), the emission rate depends on momentum and invariant mass of dilepton
    // photon emission equilibrium rate at local rest cell
    double em_eqrate = 0.;
    // photon emission viscous correction at local rest cell
    double em_visrate = 0.;
    // photon emission bulk viscous correction at local rest cell
    double em_bulkvis = 0.;
    // dilepton emission diffusion correction at local rest cell
    double em_diffrate = 0.;

    getPhotonemissionRate(Eq, M_ll, pi_zz, bulkPi, diff_factor, T, muB, rhoB_over_eplusp,
                          em_eqrate, em_visrate, em_bulkvis, em_diffrate);

    double temp_eq_sum = em_eqrate*volfrac;
    // double temp_vis_sum = em_visrate*volfrac;
    // double temp_bulkvis_sum = em_bulkvis*volfrac;
    double temp_diff_sum = em_diffrate*volfrac;

    // spectra
    dNd2pTdphidy_cell_eq = temp_eq_sum;
    dNd2pTdphidy_cell_diff = temp_eq_sum + temp_diff_sum;
    dNd2pTdphidy_cell_tot = (temp_eq_sum + temp_diff_sum);

}


// functions to get distributions in T and tau; contributions from all channels are included

void ThermalPhoton::calPhoton_SpMatrix_dTdtau(double ******dNd2pTdphidydTdtau_eq_temp, 
            double ******dNd2pTdphidydTdtau_tot_temp, double ******dNd2pTdphidydTdtau_diff_temp) {
    //#pragma omp parallel for collapse(4)
    for (int i = 0; i < nTcut; i++) {
        for (int j = 0; j < n_tau_cut; j++) {
		    for (int k = 0; k < nm; k++) {
		        for (int l = 0; l < np; l++) {
		            for (int m = 0; m < nphi; m++) {
		                for (int n = 0; n < nrapidity; n++) {
		                    dNd2pTdphidydTdtau_eq[i][j][k][l][m][n] = dNd2pTdphidydTdtau_eq_temp[i][j][k][l][m][n];
		                    dNd2pTdphidydTdtau_diff[i][j][k][l][m][n] = dNd2pTdphidydTdtau_diff_temp[i][j][k][l][m][n];
		                    dNd2pTdphidydTdtau_tot[i][j][k][l][m][n] = dNd2pTdphidydTdtau_tot_temp[i][j][k][l][m][n];
		                }
		            }
		        }
		    }
		}
	}
}


void ThermalPhoton::outputPhoton_spectra_dTdtau(string path) {
    // calculate the inverse slope of the photon spectra at T-tau interval
    ostringstream filename_SpdTdtau_eq;
    ostringstream filename_SpdTdtau_diff;
    ostringstream filename_SpdTdtau_tot;
    filename_SpdTdtau_eq 	<< path << emissionProcess_name
                         	<< "_SpdTdtau_eq.dat";
    filename_SpdTdtau_diff 	<< path << emissionProcess_name
                          	<< "_SpdTdtau_diff.dat";
    filename_SpdTdtau_tot 	<< path << emissionProcess_name
                          	<< "_SpdTdtau_tot.dat";

    ofstream ofeq(filename_SpdTdtau_eq.str().c_str());
    ofstream ofdiff(filename_SpdTdtau_diff.str().c_str());
    ofstream oftot(filename_SpdTdtau_tot.str().c_str());

    double dT = (Tcut_high - Tcut_low)/(nTcut - 1);
    double dtau = (tau_cut_high - tau_cut_low)/(n_tau_cut - 1);

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
                    double temp_dNdypTdpT_eq = 0.0;
                    double temp_dNdypTdpT_diff = 0.0;
                    double temp_dNdypTdpT_tot = 0.0;
                    for (int l = 0; l < nphi; l++) {
                        for (int irap = 0; irap < nrapidity; irap++) {
                            temp_dNdypTdpT_eq += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]
                                    *phi_weight[l]*dy);
                            temp_dNdypTdpT_diff += (
                                    dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap]
                                    *phi_weight[l]*dy);
                            temp_dNdypTdpT_tot += (
                                    dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap]
                                    *phi_weight[l]*dy);
                        }
                    }
                    ofeq << temp_dNdypTdpT_eq << "   ";
                    ofdiff << temp_dNdypTdpT_diff << "   ";
                    oftot << temp_dNdypTdpT_tot << "   ";
                }
                ofeq << endl;
                ofdiff << endl;
                oftot << endl;
            }
        }
    }

    ofeq.close();
    ofdiff.close();
    oftot.close();
}


void ThermalPhoton::calPhoton_SpvnpT_dTdtau() {
    // calculate the photon spectra and differential vn at mid-rapidity
    double eps = 1e-15;
    for (int i = 0; i < nTcut; i++) {
        for (int j = 0; j < n_tau_cut; j++) {

            for (int m = 0; m < nm; m++) {
                for (int k = 0; k < np; k++) {
                    for (int l = 0; l < nphi; l++) {
                        double weight = p[k]*p_weight[k]*phi_weight[l]*dy; // pT and phi_p integrated out
                        for (int irap = 0; irap < nrapidity; irap++) {
                            dNdydTdtau_eq[i][j][m] += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]
                                    *weight);
                            // dNdydTdtau_vis[i][j] += (
                            //         dNd2pTdphidydTdtau_vis[i][j][k][l][irap]
                            //         *weight);
                            dNdydTdtau_diff[i][j][m] += (
                                    dNd2pTdphidydTdtau_diff[i][j][m][k][l][irap]
                                    *weight);
                            dNdydTdtau_tot[i][j][m] += (
                                    dNd2pTdphidydTdtau_tot[i][j][m][k][l][irap]
                                    *weight);
                            for (int order = 0; order < norder; order++) {
                                vndTdtau_cos_eq[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]
                                    *cos(order*phi[l])*weight);
                                vndTdtau_sin_eq[i][j][m][order] += (
                                    dNd2pTdphidydTdtau_eq[i][j][m][k][l][irap]
                                    *weight*sin(order*phi[l]));
                                // vndTdtau_cos_vis[i][j][order] += (
                                //     dNd2pTdphidydTdtau_vis[i][j][k][l][irap]
                                //     *weight*cos(order*phi[l]));
                                // vndTdtau_sin_vis[i][j][order] += (
                                //     dNd2pTdphidydTdtau_vis[i][j][k][l][irap]
                                //     *weight*sin(order*phi[l]));
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
                    vndTdtau_cos_eq[i][j][m][order] = (vndTdtau_cos_eq[i][j][m][order]
                                                    /(dNdydTdtau_eq[i][j][m] + eps));
                    vndTdtau_sin_eq[i][j][m][order] = (vndTdtau_sin_eq[i][j][m][order]
                                                    /(dNdydTdtau_eq[i][j][m] + eps));
                    // vndTdtau_cos_vis[i][j][order] = (
                    //         vndTdtau_cos_vis[i][j][order]
                    //         /(dNdydTdtau_vis[i][j] + eps));
                    // vndTdtau_sin_vis[i][j][order] = (
                    //         vndTdtau_sin_vis[i][j][order]
                    //         /(dNdydTdtau_vis[i][j] + eps));
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
}


void ThermalPhoton::outputPhoton_SpvnpT_dTdtau(string path) {
    double dT = (Tcut_high - Tcut_low)/(nTcut - 1);
    double dtau = (tau_cut_high - tau_cut_low)/(n_tau_cut - 1);
    ostringstream filename_stream_dNdydTdtau_eq;
    // ostringstream filename_stream_dNdydTdtau_vis;
    ostringstream filename_stream_dNdydTdtau_diff;
    ostringstream filename_stream_dNdydTdtau_tot;

    filename_stream_dNdydTdtau_eq << path << emissionProcess_name
                                  << "_dNdydTdtau_eq.dat";
    // filename_stream_dNdydTdtau_vis << path << emissionProcess_name
    //                                << "_dNdydTdtau_vis.dat";
    filename_stream_dNdydTdtau_diff << path << emissionProcess_name
                                       << "_dNdydTdtau_diff.dat";
    filename_stream_dNdydTdtau_tot << path << emissionProcess_name
                                   << "_dNdydTdtau_tot.dat";

    ofstream fphotondNdy_eq(filename_stream_dNdydTdtau_eq.str().c_str());
    // ofstream fphotondNdy_vis(filename_stream_dNdydTdtau_vis.str().c_str());
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
          	fphotondNdy_diff << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
          	fphotondNdy_tot << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
          	for (int m = 0; m < nm; m++) {
          		double M_ll = M[m];
              	fphotondNdy_eq << dNdydTdtau_eq[i][j][m]*M_ll << "    ";
              	// fphotondNdy_vis << dNdydTdtau_vis[i][j]*M_ll << "    ";
              	fphotondNdy_diff << dNdydTdtau_diff[i][j][m]*M_ll << "    ";
              	fphotondNdy_tot << dNdydTdtau_tot[i][j][m]*M_ll << "    ";
          	}
          	fphotondNdy_eq << endl;
           	// fphotondNdy_vis << endl;
           	fphotondNdy_diff << endl;
           	fphotondNdy_tot << endl;
       	}
	}
    fphotondNdy_eq.close();
    // fphotondNdy_vis.close();
    fphotondNdy_diff.close();
    fphotondNdy_tot.close();

    for(int order = 1; order < norder; order++)
    {
       	ostringstream filename_stream_vncosdTdtau_eq;
       	// ostringstream filename_stream_vncosdTdtau_vis;
       	ostringstream filename_stream_vncosdTdtau_diff;
       	ostringstream filename_stream_vncosdTdtau_tot;
       	ostringstream filename_stream_vnsindTdtau_eq;
       	// ostringstream filename_stream_vnsindTdtau_vis;
       	ostringstream filename_stream_vnsindTdtau_diff;
       	ostringstream filename_stream_vnsindTdtau_tot;
       	filename_stream_vncosdTdtau_eq 	<< path << emissionProcess_name
                                      	<< "_v_" << order
                                      	<< "_cos_dTdtau_eq.dat";
       	// filename_stream_vncosdTdtau_vis << path << emissionProcess_name
       	//                                 << "_v_" << order
       	//                                 << "_cos_dTdtau_vis.dat";
       	filename_stream_vncosdTdtau_diff 	<< path << emissionProcess_name
                                           	<< "_v_" << order
                                           	<< "_cos_dTdtau_diff.dat";
       	filename_stream_vncosdTdtau_tot << path << emissionProcess_name
                                       	<< "_v_" << order
                                       	<< "_cos_dTdtau_tot.dat";
       	filename_stream_vnsindTdtau_eq 	<< path << emissionProcess_name
                                      	<< "_v_" << order
                                      	<< "_sin_dTdtau_eq.dat";
       	// filename_stream_vnsindTdtau_vis << path << emissionProcess_name
       	//                                 << "_v_" << order
       	//                                 << "_sin_dTdtau_vis.dat";
       	filename_stream_vnsindTdtau_diff << path << emissionProcess_name
                                      	<< "_v_" << order
                                        << "_sin_dTdtau_diff.dat";
       	filename_stream_vnsindTdtau_tot 	<< path << emissionProcess_name
                                       	<< "_v_" << order
                                       	<< "_sin_dTdtau_tot.dat";

       	ofstream fphotonvncos_eq(filename_stream_vncosdTdtau_eq.str().c_str());
       	// ofstream fphotonvncos_vis(filename_stream_vncosdTdtau_vis.str().c_str());
       	ofstream fphotonvncos_diff(filename_stream_vncosdTdtau_diff.str().c_str());
       	ofstream fphotonvncos_tot(filename_stream_vncosdTdtau_tot.str().c_str());
       	ofstream fphotonvnsin_eq(filename_stream_vnsindTdtau_eq.str().c_str());
       	// ofstream fphotonvnsin_vis(filename_stream_vnsindTdtau_vis.str().c_str());
       	ofstream fphotonvnsin_diff(filename_stream_vnsindTdtau_diff.str().c_str());
       	ofstream fphotonvnsin_tot(filename_stream_vnsindTdtau_tot.str().c_str());
       	for(int i = 0; i < nTcut; i++)
       	{
          	double T_local = Tcut_low + i*dT;
          	for(int j = 0; j < n_tau_cut; j++)
      		{
             	double tau_local = tau_cut_low + j*dtau;
             	fphotonvncos_eq << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
             	fphotonvncos_diff << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
             	fphotonvncos_tot << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
             	fphotonvnsin_eq << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
             	fphotonvnsin_diff << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
             	fphotonvnsin_tot << scientific << setw(18) << setprecision(8) 
                   	<< T_local << "   "  << tau_local << "   ";
             	for (int m = 0; m < nm; m++) {
                 	fphotonvncos_eq << vndTdtau_cos_eq[i][j][m][order] << "    ";
                 	// fphotonvncos_vis << vndTdtau_cos_vis[i][j][order] << "    ";
                 	fphotonvncos_diff << vndTdtau_cos_diff[i][j][m][order] << "    ";
                 	fphotonvncos_tot << vndTdtau_cos_tot[i][j][m][order] << "    ";
                 	fphotonvnsin_eq << vndTdtau_sin_eq[i][j][m][order] << "    ";
                 	// fphotonvnsin_vis << vndTdtau_sin_vis[i][j][order] << "    ";
                 	fphotonvnsin_diff << vndTdtau_sin_diff[i][j][m][order] << "    ";
                 	fphotonvnsin_tot << vndTdtau_sin_tot[i][j][m][order] << "    ";
             	}
             	fphotonvncos_eq << endl;
             	// fphotonvncos_vis << endl;
             	fphotonvncos_diff << endl;
             	fphotonvncos_tot << endl;
             	fphotonvnsin_eq << endl;
             	// fphotonvnsin_vis << endl;
             	fphotonvnsin_diff << endl;
             	fphotonvnsin_tot << endl;
          	}
       	}
       	fphotonvncos_eq.close();
       	fphotonvnsin_eq.close();
       	// fphotonvncos_vis.close();
       	// fphotonvnsin_vis.close();
       	fphotonvncos_diff.close();
       	fphotonvnsin_diff.close();
       	fphotonvncos_tot.close();
       	fphotonvnsin_tot.close();
    }
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


// //! this function is used the most frequent one,
// //! it needs to be as fast as possible
// void ThermalPhoton::interpolation2D_bilinear(
//         double varX, vector<double> &varY,
//         double** Table2D_ptr, vector<double> &results) {
//     const int Y_length = varY.size();
//     double varX_min = EmissionrateTb_Xmin;
//     double varY_min = EmissionrateTb_Ymin;
//     double dX = EmissionrateTb_dX;
//     double dY = EmissionrateTb_dY;
//     //int tb_sizeX = EmissionrateTb_sizeX;
//     int tb_sizeY = EmissionrateTb_sizeY;

//     double dXdYinverse = 1/(dX*dY);
//     double dYinverse = 1/dY;

//     int idx_Y;
//     double f00, f10, f01, f11;
//     double ddy;
//     double dYhigh;

//     int idx_X;
//     double ddx;
//     double dXhigh;
//     idx_X  = static_cast<int>((varX-varX_min)/dX);
//     ddx    = (varX - (varX_min + idx_X*dX));
//     dXhigh = (dX - ddx)*dXdYinverse;
//     ddx    = ddx*dXdYinverse;

//     for (int ii = 0; ii < Y_length; ii++) {
//         idx_Y = static_cast<int>((varY[ii] - varY_min)*dYinverse);

//         //if it is out the table, doing linear extrapolation
//         //if (idx_X<0 || idx_X>=tb_sizeX-1) {
//         //    if (Iwarning == 1) {
//         //        cout << "interpolation2D_bilinear: varX out of bounds!"
//         //             << endl;
//         //        cout << "varX= " << varX << "idx_X= " << idx_X <<endl;
//         //    }
//         //    if (idx_X < 0)
//         //        idx_X = 0;
//         //    if (idx_X >= tb_sizeX-1)
//         //        idx_X = tb_sizeX - 1;
//         //}

//         //if (Iwarning == 1) {
//         //    cout << "interpolation2D_bilinear: varY out of bounds!" << endl;
//         //    cout << "varY= " << varY << "idx_Y= "<< idx_Y << endl;
//         //}
//         if (idx_Y < 0)
//             idx_Y = 0;
//         if (idx_Y >= tb_sizeY-1)
//             idx_Y = tb_sizeY - 2;

//         ddy = varY[ii] - EmissionrateTb_Yidxptr[idx_Y];
//         dYhigh = dY - ddy;

//         f00 = Table2D_ptr[idx_X][idx_Y];
//         f10 = Table2D_ptr[idx_X+1][idx_Y];
//         f01 = Table2D_ptr[idx_X][idx_Y+1];
//         f11 = Table2D_ptr[idx_X+1][idx_Y+1];

//         results[ii] = (
//             ((f00*dYhigh + f01*ddy)*dXhigh + (f10*dYhigh + f11*ddy)*ddx));
//      }
// }

