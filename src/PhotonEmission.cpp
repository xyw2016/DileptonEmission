// Copyright 2016 Chun Shen
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
//#include <omp.h>
#include <vector>
#include <unordered_set>
#ifndef _OPENMP
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif

#include "Hydroinfo_h5.h"
#include "ThermalPhoton.h"
#include "QGP_LO_analytic.h"
#include "QGP_LO.h"
#include "QGP_NLO.h"
#include "tensor_trans.h"
#include "PhotonEmission.h"
#include "ParameterReader.h"
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
using TENSORTRANSFORM::getTransverseflow_u_mu_low;
using TENSORTRANSFORM::lorentz_boost_matrix;

using PhysConsts::hbarC;
using PhysConsts::eps;


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
    InitializePhotonEmissionRateTables();

    dNd2pTdphidy_eq = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidy_eqT = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidy_eqL = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidy_visc = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidy_diff = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTdphidy_tot = createA4DMatrix(nm, np, nphi, nrapidity, 0.);
    dNd2pTd2M_eq = createA2DMatrix(nm, np, 0);
    dNd2pTd2M_eqT = createA2DMatrix(nm, np, 0);
    dNd2pTd2M_eqL = createA2DMatrix(nm, np, 0);
    dNd2pTd2M_visc = createA2DMatrix(nm, np, 0);
    dNd2pTd2M_diff = createA2DMatrix(nm, np, 0);
    dNd2pTd2M_tot = createA2DMatrix(nm, np, 0);
    dNd2Mdy_eq.resize(nm, 0);
    dNd2Mdy_eqT.resize(nm, 0);
    dNd2Mdy_eqL.resize(nm, 0);
    dNd2Mdy_visc.resize(nm, 0);
    dNd2Mdy_diff.resize(nm, 0);
    dNd2Mdy_tot.resize(nm, 0);

    vnpT_cos_eq = createA3DMatrix(norder, nm, np, 0.);
    vnpT_sin_eq = createA3DMatrix(norder, nm, np, 0.);
    vnpT_cos_visc = createA3DMatrix(norder, nm, np, 0.);
    vnpT_sin_visc = createA3DMatrix(norder, nm, np, 0.);
    vnpT_cos_diff = createA3DMatrix(norder, nm, np, 0.);
    vnpT_sin_diff = createA3DMatrix(norder, nm, np, 0.);
    vnpT_cos_tot = createA3DMatrix(norder, nm, np, 0.);
    vnpT_sin_tot = createA3DMatrix(norder, nm, np, 0.);

    vn_cos_eq = createA2DMatrix(norder, nm, 0.);
    vn_sin_eq = createA2DMatrix(norder, nm, 0.);
    vn_cos_visc = createA2DMatrix(norder, nm, 0.);
    vn_sin_visc = createA2DMatrix(norder, nm, 0.);
    vn_cos_diff = createA2DMatrix(norder, nm, 0.);
    vn_sin_diff = createA2DMatrix(norder, nm, 0.);
    vn_cos_tot = createA2DMatrix(norder, nm, 0.);
    vn_sin_tot = createA2DMatrix(norder, nm, 0.);
     
    if (differential_flag == 1) {
        dNd2pTdphidydTdtau_eq = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_visc = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_diff = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);
        dNd2pTdphidydTdtau_tot = createA6DMatrix(nTcut, n_tau_cut, nm, np, nphi, nrapidity, 0.);

        dNd2pTdphidydTdtau_eq_all = createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
        dNd2pTdphidydTdtau_visc_all = createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
        dNd2pTdphidydTdtau_diff_all = createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
        dNd2pTdphidydTdtau_tot_all = createA3DMatrix(nTcut, n_tau_cut, CORES*nrapidity*np*nphi*nm, 0.);
    }

}


PhotonEmission::~PhotonEmission() {

    deleteA4DMatrix(dNd2pTdphidy_eq, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_eqT, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_eqL, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_visc, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_diff, nm, np, nphi);
    deleteA4DMatrix(dNd2pTdphidy_tot, nm, np, nphi);
    deleteA2DMatrix(dNd2pTd2M_eq, nm);
    deleteA2DMatrix(dNd2pTd2M_eqT, nm);
    deleteA2DMatrix(dNd2pTd2M_eqL, nm);
    deleteA2DMatrix(dNd2pTd2M_visc, nm);
    deleteA2DMatrix(dNd2pTd2M_diff, nm);
    deleteA2DMatrix(dNd2pTd2M_tot, nm);

    deleteA3DMatrix(vnpT_cos_eq, norder, nm);
    deleteA3DMatrix(vnpT_sin_eq, norder, nm);
    deleteA3DMatrix(vnpT_cos_visc, norder, nm);
    deleteA3DMatrix(vnpT_sin_visc, norder, nm);
    deleteA3DMatrix(vnpT_cos_diff, norder, nm);
    deleteA3DMatrix(vnpT_sin_diff, norder, nm);
    deleteA3DMatrix(vnpT_cos_tot, norder, nm);
    deleteA3DMatrix(vnpT_sin_tot, norder, nm);

    deleteA2DMatrix(vn_cos_eq, norder);
    deleteA2DMatrix(vn_sin_eq, norder);
    deleteA2DMatrix(vn_cos_visc, norder);
    deleteA2DMatrix(vn_sin_visc, norder);
    deleteA2DMatrix(vn_cos_diff, norder);
    deleteA2DMatrix(vn_sin_diff, norder);
    deleteA2DMatrix(vn_cos_tot, norder);
    deleteA2DMatrix(vn_sin_tot, norder);

    if (differential_flag == 1) {
        nTcut = paraRdr->getVal("nTcut");
        n_tau_cut = paraRdr->getVal("n_tau_cut");
        deleteA6DMatrix(dNd2pTdphidydTdtau_eq, nTcut, n_tau_cut, nm, np, nphi);
        deleteA6DMatrix(dNd2pTdphidydTdtau_visc, nTcut, n_tau_cut, nm, np, nphi);
        deleteA6DMatrix(dNd2pTdphidydTdtau_diff, nTcut, n_tau_cut, nm, np, nphi);
        deleteA6DMatrix(dNd2pTdphidydTdtau_tot, nTcut, n_tau_cut, nm, np, nphi);

        deleteA3DMatrix(dNd2pTdphidydTdtau_eq_all, nTcut, n_tau_cut);
        deleteA3DMatrix(dNd2pTdphidydTdtau_visc_all, nTcut, n_tau_cut);
        deleteA3DMatrix(dNd2pTdphidydTdtau_diff_all, nTcut, n_tau_cut);
        deleteA3DMatrix(dNd2pTdphidydTdtau_tot_all, nTcut, n_tau_cut);
    }
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

    T_test =  paraRdr->getVal("T_test");
    muB_test =  paraRdr->getVal("muB_test");
    rhoB_eplusp_test =  paraRdr->getVal("rhoB_eplusp_test");
    inv_eplusp_test =  paraRdr->getVal("inv_eplusp_test");

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
    cout << "T_sw = " << T_sw_low <<  " to " << T_sw_high << " GeV."<< endl;
    cout << endl;

    cout << "dilepton_emission_rate = " << emission_rate_flag << endl;

    cout << "Photon momentum: " << paraRdr->getVal("photon_q_i")
         << " to " << paraRdr->getVal("photon_q_f") << " GeV, "
         << "n_q =" << np << endl;
    cout << "Photon momentum angles: " << paraRdr->getVal("photon_phi_q_i")
         << " to " << paraRdr->getVal("photon_phi_q_f")
         << ", n_phi=" << nphi << endl;
    cout << "Photon momentum rapidity: " << paraRdr->getVal("photon_y_i")
         << " to " << paraRdr->getVal("photon_y_f")
         << ", n_y =" << nrapidity << endl;
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
        cout << "rhoB_eplusp_test = " << rhoB_eplusp_test  << " fm^-1." << endl;
    }

}


void PhotonEmission::InitializePhotonEmissionRateTables() {

    if(emission_rate_flag==0){
        dilepton_QGP_thermal = std::unique_ptr<ThermalPhoton>(
            new QGP_LO_analytic(paraRdr, "QGP_LO_analytic_total"));
    } else if(emission_rate_flag==1) { // To use LO analytical form
        dilepton_QGP_thermal = std::unique_ptr<ThermalPhoton>(
            new QGP_LO(paraRdr, "QGP_LO_total"));
    } else if(emission_rate_flag==2) { // To use NLO table from GJ
        dilepton_QGP_thermal = std::unique_ptr<ThermalPhoton>(
            new QGP_NLO(paraRdr, "QGP_NLO_total"));
        dilepton_QGP_thermal->readEmissionrateFromFile(true);// true to read in the NLO emission table
    }
}

double PhotonEmission::suppression_factor(double tau){
    //double lambda= 0.651; //tau_chem = 1.5
    //double lambda= 0.434; //tau_chem = 1.0
    double sig_lambda = paraRdr->getVal("sig_lambda");
    

    return 1-exp(- tau/sig_lambda);
}


void PhotonEmission::calPhotonemission_3d(void *hydroinfo_ptr_in) {

    // hydro data read in main.cpp
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;
    hydroinfo_MUSIC_ptr =
                reinterpret_cast<Hydroinfo_MUSIC*>(hydroinfo_ptr_in);
    
    int hydro_flag = paraRdr->getVal("hydro_flag");

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
    
    double volume_base = Nskip_tau*dtau*Nskip_x*dx*Nskip_x*dx*Nskip_eta*deta;


    tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
    tau_max = hydroinfo_MUSIC_ptr->get_hydro_tau_max();
    tau_cut_low = tau0 - dtau;
    tau_cut_high = tau_max + dtau;
    T_cutlow = hydroinfo_MUSIC_ptr->get_hydro_T_min();
    T_cuthigh =hydroinfo_MUSIC_ptr->get_hydro_T_max();

    // output results in (T, tau)
    double dT_cut = (T_cuthigh - T_cutlow)/(nTcut - 1);
    double dtau_cut = (tau_cut_high - tau_cut_low)/(n_tau_cut - 1);

    if (differential_flag == 1){
        if (tau_max>tau_cut_high||tau0<tau_cut_low){
            printf("Warning: proper time out of [tau_cut_low, tau_cut_high].\n");
            printf("The boundaries should be enlarged to encolse all hydro cells...\n");
        }
    }
    
    // number of fluid cells
    long int number_of_cells =(
                hydroinfo_MUSIC_ptr->get_number_of_fluid_cells_3d());
    cout << "number of cells:" << number_of_cells << endl;

    // multi-threads setup
    long FO_chunk = number_of_cells / CORES;
    long remainder = number_of_cells  -  CORES * FO_chunk;

    cout << "Number of cores : " << CORES << endl;
    cout << "Chunk size = " << FO_chunk << endl;
    cout << "Remainder cells = " << remainder << endl;
    cout << "Volume base = " << volume_base << endl;

    if(remainder != 0) FO_chunk++;

    cout << "----------------------------------------" << endl;

    // arrays to store values across cores
    double *dNd2pTdphidy_eq_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
    double *dNd2pTdphidy_eqT_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
    double *dNd2pTdphidy_eqL_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
    double *dNd2pTdphidy_visc_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
    double *dNd2pTdphidy_diff_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));
    double *dNd2pTdphidy_tot_all = (double*)calloc(CORES * nrapidity*np*nphi*nm, sizeof(double));

    // main loop begins ...
    // loop over all fluid cells
    // subdivide bite size chunks of freezeout surface across cores
    int ncells = 0;
    double tau_now = 0.0;

    #pragma omp parallel for reduction(+:ncells)
    for(long n = 0; n < CORES; n++)
    {
        long endFO = FO_chunk;

        int processedCells = 0;  // Counter for processed cells
        int printInterval = 0.25 * endFO;  // Print at 25% intervals

        for(long int icell = 0; icell < endFO; icell++)  // cell index inside each chunk
        {
            if((icell == endFO - 1) && (remainder != 0) && (n > remainder - 1)) continue;
            
	    long cell_id = n  +  icell * CORES;

            // fluid velocity in Minkowski
            double flow_u_mu_low[4];
            double flow_u_mu_Min[4]; 
            // dilepton 4-momentum in Minkowski coordinates
            double p_lab_local[4];
            double p_lab_Min[4]; 

            fluidCell_3D_new fluidCellptr;

            hydroinfo_MUSIC_ptr->get_hydro_cell_info_3d(cell_id, fluidCellptr);

            double tau_local = tau0 + fluidCellptr.itau*dtau;
            double eta_local = -eta_max + fluidCellptr.ieta*deta;

            #ifdef _OPENMP
            processedCells++;  // Increment the counter for processed cells

            if (processedCells % printInterval == 0)
            {
                double progress = (processedCells * 100) / FO_chunk;
                int threadNum = omp_get_thread_num();
                std::cout << "Thread " << threadNum << ": Processing " << progress << "% of cells" << std::endl;
            }
            #else

            if (fabs(tau_now - tau_local) > 1e-10) {
                tau_now = tau_local;
                cout << "Calculating tau = " << setw(4) << setprecision(3)
                     << tau_now << " fm/c..." << endl;
             }
            #endif         

            // wxy
            // volume element: tau*dtau*dx*dy*deta,
            double volume = tau_local*volume_base;
            if(hydro_flag == 22){
                volume = volume_base;
            }
            

            double ed_local = fluidCellptr.ed;
            double pd_local = fluidCellptr.pressure;
            double temp_local = fluidCellptr.temperature;
            double temp_inv = 1/temp_local;


            // note that when turn_on_muB_=0, muB and rhoB are set to zero as follows
            double muB_local = fluidCellptr.muB;
            double rhoB_local = fluidCellptr.rhoB;
            double inv_eplusp = 1./(ed_local+pd_local);
            double rhoB_over_eplusp = rhoB_local*inv_eplusp;

            double T_sw = 0.166 - 0.139 * pow(muB_local, 2) - 0.053 * pow(muB_local, 4); // Cleymans et al

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

            // fluid cell is out of interest, temperature below T_dec or eta larger than ETAmax
            if (hydro_flag==2 && (temp_local < T_dec || eta_local > ETAmax))
                continue;
            if (hydro_flag==22 && (temp_local < T_dec || eta_local > ETAmax))
                continue;

            if (differential_flag == 1){
                if (temp_local > T_cuthigh||temp_local < T_cutlow||tau_local>tau_cut_high||tau_local<tau_cut_low){
                    printf("Warning: local temperature or proper time out of [T_cutlow, T_cuthigh] or [tau_cut_low, tau_cut_high].\n");
                    printf("The boundaries should be enlarged to encolse all hydro cells...\n");
                    continue;
                }
            }
            
	    if (temp_local > T_sw) {
                // total fluid cells of QGP emission
                // #pragma omp atomic update
                ncells++;
            }

            // indices for the output in (T, tau)
            int idx_T = (int)std::floor((temp_local - T_cutlow)/dT_cut + eps);
            int idx_tau = (int)std::floor((tau_local - tau_cut_low)/dtau_cut + eps);

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

            double utau = sqrt(1. + ux*ux + uy*uy + ueta*ueta);

            flow_u_mu_low[0] = utau;
            flow_u_mu_low[1] = -ux;
            flow_u_mu_low[2] = -uy;
            flow_u_mu_low[3] = -ueta;
            

            if( hydro_flag == 2)
            {
            double cosh_eta = cosh(eta_local);
            double sinh_eta = sinh(eta_local);

            // 4 velocity in Minkowski in lab frame
            flow_u_mu_Min[0] = cosh_eta*utau + sinh_eta*ueta;
            flow_u_mu_Min[1] = ux;
            flow_u_mu_Min[2] = uy;
            flow_u_mu_Min[3] = sinh_eta*utau + cosh_eta*ueta;
            
            
            }
 
            if( hydro_flag == 22)
            {
             //wxy    
            // 4 velocity in Minkowski in lab frame
            flow_u_mu_Min[0] = utau ;
            flow_u_mu_Min[1] = ux;
            flow_u_mu_Min[2] = uy;
            flow_u_mu_Min[3] = ueta;  
            }

            

            // Lorentz boost matrix from lab to local rest frame
            std::vector<std::vector<double>> lambda_munu(4, std::vector<double>(4, 0.0));
            lorentz_boost_matrix(lambda_munu, flow_u_mu_Min[0], 
                flow_u_mu_Min[1], flow_u_mu_Min[2], flow_u_mu_Min[3]);
            
            
            // Wij is reduced variables Wij/(e+P), Wieta = tau*Wieta
            double pi11 = fluidCellptr.pi11;
            double pi12 = fluidCellptr.pi12;
            double pi13 = fluidCellptr.pi13;
            double pi22 = fluidCellptr.pi22;
            double pi23 = fluidCellptr.pi23;
            // reconstruct all other components of the shear stress tensor
            double pi01 = (ux*pi11 + uy*pi12 + ueta*pi13)/utau;
            double pi02 = (ux*pi12 + uy*pi22 + ueta*pi23)/utau;
            double pi33 = (utau*(ux*pi01 + uy*pi02) - utau*utau*(pi11 + pi22)
                           + ueta*(ux*pi13 + uy*pi23))/(utau*utau - ueta*ueta);
            double pi00 = pi11 + pi22 + pi33;
            double pi03 = (ux*pi13 + uy*pi23 + ueta*pi33)/utau;

            double bulkPi_local = fluidCellptr.bulkPi;

            // qeta = tau*qeta, qi is qi/kappa_hat
            double qx = fluidCellptr.qx;
            double qy = fluidCellptr.qy;
            double qeta = fluidCellptr.qz;
            double qtau = (ux*qx + uy*qy + ueta*qeta)/utau;

            // prefactors in the dissipative correction
            double Cq = 0.97; //https://journals.aps.org/prc/pdf/10.1103/PhysRevC.89.034904 under eq.6
            double prefactor_diff = 1./(4.*pow(2*M_PI, 5)) * temp_inv/pow(hbarC, 4);
            double prefactor_visc = Cq*prefactor_diff;

            double spsfactor = 1.0;
            if(hydro_flag ==22 && test_code_flag != 1  ){
                //wxy
                double cell_t = tau0 + fluidCellptr.itau*dtau;
                double cell_z = -eta_max + fluidCellptr.ieta*deta;
                double cell_tau = 0.0;
                if( fabs(cell_t) >= fabs(cell_z))
                {
                    cell_tau = sqrt(cell_t*cell_t - cell_z*cell_z);
                }

                spsfactor = suppression_factor(cell_tau);
		spsfactor = spsfactor;
                // if(spsfactor > 0)
                // std::cout <<"t: "<< cell_t <<"  z: "<< cell_z <<" tau: "<< cell_tau<<" factor: "
                //           << spsfactor<< "Tem: "<< temp_local<<std::endl;
            }

         
            
            // photon momentum loops
            // #pragma omp parallel for collapse(4) private(p_lab_Min, p_lab_local)
            for (int k = 0; k < nrapidity; k++) {
                for (int m = 0; m < nphi; m++) {
                    for (int l = 0; l < np; l++) {
                        for (int j = 0; j < nm; j++) {

                            int i3 = (j+(l+(m+nphi * k) * np) * nm);
                            // p_q is p_T magnitude array
                            double M_T = sqrt(p_q[l]*p_q[l]+M_ll[j]*M_ll[j]); // transverse mass

                            y_q[k]=0.0;
                            //cos_phiq[m]=1.0;sin_phiq[m]=0.0;
                            double cosh_y = cosh(y_q[k]);
                            double sinh_y = sinh(y_q[k]);

                              
                            double cosh_y_minus_eta = cosh(y_q[k] - eta_local);
                            double sinh_y_minus_eta = sinh(y_q[k] - eta_local);
                             // from Minkowski to Milne in lab frame
                            p_lab_local[0] = M_T * cosh_y_minus_eta;
                            p_lab_local[1] = p_q[l]*cos_phiq[m];
                            p_lab_local[2] = p_q[l]*sin_phiq[m];
                            p_lab_local[3] = M_T * sinh_y_minus_eta; // tau*p^eta
                           

                            if(hydro_flag == 22)
                            {

                                p_lab_local[0] = M_T * cosh_y;
                                p_lab_local[1] = p_q[l]*cos_phiq[m];
                                p_lab_local[2] = p_q[l]*sin_phiq[m];
                                p_lab_local[3] = M_T * sinh_y; // tau*p^eta

                            }

                            

                            // Minkowski four momentum in lab frame
                            p_lab_Min[0] = M_T * cosh_y;
                            p_lab_Min[1] = p_q[l]*cos_phiq[m];
                            p_lab_Min[2] = p_q[l]*sin_phiq[m];
                            p_lab_Min[3] = M_T * sinh_y;
                            
                           

                            // Minkowski four momentum from Lab frame to LRF frame
                            double p_Min_lrf[4];
                            for (int j = 0; j < 4; j++) {
                                p_Min_lrf[j] = 0.;
                                for (int i = 0; i < 4; i++) {
                                    p_Min_lrf[j] += lambda_munu[j][i]*p_lab_Min[i];
                                }                            
                            }

                            
                            double pvec_lrf, pvec3, pvec5; // Minkowski spatial magnitude |vec p| and |vec p|^3
                            pvec_lrf = sqrt(p_Min_lrf[1]*p_Min_lrf[1] + p_Min_lrf[2]*p_Min_lrf[2] + p_Min_lrf[3]*p_Min_lrf[3]);
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
                            double diff_dot_p = qtau*p_lab_local[0] - qx*p_lab_local[1] 
                                                - qy*p_lab_local[2] - qeta*p_lab_local[3];
                            // validation setup
                            if(test_code_flag==1){
                                pi_photon = 0.0;
                                diff_dot_p = 0.0;
                            }

                            
                            double Eq_localrest_Tb = p_Min_lrf[0];

                            double M2 = M_ll[j] * M_ll[j];
                            double visc_fac = prefactor_visc * pi_photon * M2 / pvec5;
                            double diff_fac = prefactor_diff * diff_dot_p * M2 / pvec3;
                            double bulkPi_fac = bulkPi_local;

                            double dNd2pTdphidy_cell_eq = 0.;
                            double dNd2pTdphidy_cell_eqT = 0.;
                            double dNd2pTdphidy_cell_eqL = 0.;
                            double dNd2pTdphidy_cell_visc = 0.;
                            double dNd2pTdphidy_cell_diff = 0.;
                            double dNd2pTdphidy_cell_tot = 0.;

                            // begin to calculate thermal photon emission
                            if (hydro_flag==2 && temp_local > T_sw) {
                                // QGP emission
                                double QGP_fraction = 1.0;
                                muB_local = turn_on_muB_*muB_local;
                                rhoB_over_eplusp = turn_on_muB_*rhoB_over_eplusp;
                                dilepton_QGP_thermal->calThermalPhotonemission_3d(
                                    Eq_localrest_Tb, M_ll[j], visc_fac, bulkPi_fac, diff_fac,
                                    temp_local, muB_local, inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,
                                    dNd2pTdphidy_cell_eq, dNd2pTdphidy_cell_eqT, dNd2pTdphidy_cell_eqL, dNd2pTdphidy_cell_visc,
                                    dNd2pTdphidy_cell_diff, dNd2pTdphidy_cell_tot);
                            }
                            if (hydro_flag == 22 && temp_local > T_sw ){
                                double QGP_fraction = 1.0;
                                muB_local = turn_on_muB_*muB_local;
                                rhoB_over_eplusp = turn_on_muB_*rhoB_over_eplusp;
                                dilepton_QGP_thermal->calThermalPhotonemission_3d(
                                    Eq_localrest_Tb, M_ll[j], visc_fac, bulkPi_fac, diff_fac,
                                    temp_local, muB_local, inv_eplusp, rhoB_over_eplusp, volume, QGP_fraction,
                                    dNd2pTdphidy_cell_eq, dNd2pTdphidy_cell_eqT, dNd2pTdphidy_cell_eqL, dNd2pTdphidy_cell_visc,
                                    dNd2pTdphidy_cell_diff, dNd2pTdphidy_cell_tot);

                            }

                            // add contributions from QGP and Hadronic matter, etc
                            // these are dN/(MdM pTdpTdphi dy)
                            dNd2pTdphidy_eq_all[n + CORES * i3] += dNd2pTdphidy_cell_eq*spsfactor;
                            dNd2pTdphidy_eqT_all[n + CORES * i3] += dNd2pTdphidy_cell_eqT*spsfactor;
                            dNd2pTdphidy_eqL_all[n + CORES * i3] += dNd2pTdphidy_cell_eqL*spsfactor;
                            dNd2pTdphidy_visc_all[n + CORES * i3] += dNd2pTdphidy_cell_visc*spsfactor;
                            dNd2pTdphidy_diff_all[n + CORES * i3] += dNd2pTdphidy_cell_diff*spsfactor;
                            dNd2pTdphidy_tot_all[n + CORES * i3] += dNd2pTdphidy_cell_tot*spsfactor;

                            if (differential_flag == 1) {
                                dNd2pTdphidydTdtau_eq_all[idx_T][idx_tau][n+CORES*i3] += dNd2pTdphidy_cell_eq*spsfactor;
                                dNd2pTdphidydTdtau_visc_all[idx_T][idx_tau][n+CORES*i3] += dNd2pTdphidy_cell_visc*spsfactor;
                                dNd2pTdphidydTdtau_diff_all[idx_T][idx_tau][n+CORES*i3] += dNd2pTdphidy_cell_diff*spsfactor;
                                dNd2pTdphidydTdtau_tot_all[idx_T][idx_tau][n+CORES*i3] += dNd2pTdphidy_cell_tot*spsfactor;
                            }

                        } // M_ll
                    } // p_T
                } // phi_p
            } // y
        } // fluid cell
    }// cores

    #pragma omp parallel for collapse(4)
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

                    #pragma omp simd reduction(+:dN_pTdpTdphidy_eq_tmp,dN_pTdpTdphidy_eqT_tmp,dN_pTdpTdphidy_eqL_tmp,\
                            dN_pTdpTdphidy_visc_tmp,dN_pTdpTdphidy_diff_tmp,dN_pTdpTdphidy_tot_tmp)
                    for(long n = 0; n < CORES; n++)
                    {
                        dN_pTdpTdphidy_eq_tmp += dNd2pTdphidy_eq_all[n+CORES*i3];
                        dN_pTdpTdphidy_eqT_tmp += dNd2pTdphidy_eqT_all[n+CORES*i3];
                        dN_pTdpTdphidy_eqL_tmp += dNd2pTdphidy_eqL_all[n+CORES*i3];
                        dN_pTdpTdphidy_visc_tmp += dNd2pTdphidy_visc_all[n+CORES*i3];
                        dN_pTdpTdphidy_diff_tmp += dNd2pTdphidy_diff_all[n+CORES*i3];
                        dN_pTdpTdphidy_tot_tmp += dNd2pTdphidy_tot_all[n+CORES*i3];
                    } // sum over the cores

                    dNd2pTdphidy_eq[j][l][m][k] = dN_pTdpTdphidy_eq_tmp;
                    dNd2pTdphidy_eqT[j][l][m][k] = dN_pTdpTdphidy_eqT_tmp;
                    dNd2pTdphidy_eqL[j][l][m][k] = dN_pTdpTdphidy_eqL_tmp;
                    dNd2pTdphidy_visc[j][l][m][k] = dN_pTdpTdphidy_visc_tmp;
                    dNd2pTdphidy_diff[j][l][m][k] = dN_pTdpTdphidy_diff_tmp;
                    dNd2pTdphidy_tot[j][l][m][k] = dN_pTdpTdphidy_tot_tmp;

                    // distribution in (T, tau)
                    if (differential_flag == 1) {
                        for (int iT = 0; iT < nTcut; iT++) {
                             for (int it = 0; it < n_tau_cut; it++) {

                                double dN_pTdpTdphidydTdtau_eq_tmp = 0.0;
                                double dN_pTdpTdphidydTdtau_visc_tmp = 0.0;
                                double dN_pTdpTdphidydTdtau_diff_tmp = 0.0;
                                double dN_pTdpTdphidydTdtau_tot_tmp = 0.0;

                                #pragma omp simd reduction(+:dN_pTdpTdphidydTdtau_eq_tmp,dN_pTdpTdphidydTdtau_visc_tmp,\
                                    dN_pTdpTdphidydTdtau_diff_tmp,dN_pTdpTdphidydTdtau_tot_tmp)
                                for(long n = 0; n < CORES; n++) {
                                    dN_pTdpTdphidydTdtau_eq_tmp += dNd2pTdphidydTdtau_eq_all[iT][it][n+CORES*i3];
                                    dN_pTdpTdphidydTdtau_visc_tmp += dNd2pTdphidydTdtau_visc_all[iT][it][n+CORES*i3];
                                    dN_pTdpTdphidydTdtau_diff_tmp += dNd2pTdphidydTdtau_diff_all[iT][it][n+CORES*i3];
                                    dN_pTdpTdphidydTdtau_tot_tmp += dNd2pTdphidydTdtau_tot_all[iT][it][n+CORES*i3];
                                }

                              

                                dNd2pTdphidydTdtau_eq[iT][it][j][l][m][k] = dN_pTdpTdphidydTdtau_eq_tmp;
                                dNd2pTdphidydTdtau_visc[iT][it][j][l][m][k] = dN_pTdpTdphidydTdtau_visc_tmp;
                                dNd2pTdphidydTdtau_diff[iT][it][j][l][m][k] = dN_pTdpTdphidydTdtau_diff_tmp;
                                dNd2pTdphidydTdtau_tot[iT][it][j][l][m][k] = dN_pTdpTdphidydTdtau_tot_tmp;

                            }
                        }
                    }

                } // M_ll
            } // p_T
        } // phi_p
    } // y


    // Total number of elements that satisfy the condition
    // int total_count = ncells;
    printf("Cells above T_sw_high=%d...\n", ncells);

    // free memory
    free(dNd2pTdphidy_eq_all);
    free(dNd2pTdphidy_eqT_all);
    free(dNd2pTdphidy_eqL_all);
    free(dNd2pTdphidy_visc_all);
    free(dNd2pTdphidy_diff_all);
    free(dNd2pTdphidy_tot_all);
}


void PhotonEmission::calPhoton_SpvnpT_individualchannel() {
    // dilepton_QGP_thermal->calPhoton_SpvnpT_shell();
    if (differential_flag == 1) {
        dilepton_QGP_thermal->calPhoton_SpMatrix_dTdtau(dNd2pTdphidydTdtau_eq, 
            dNd2pTdphidydTdtau_visc, dNd2pTdphidydTdtau_diff, dNd2pTdphidydTdtau_tot);
        dilepton_QGP_thermal->calPhoton_Spectra_dTdtau();
        dilepton_QGP_thermal->calPhoton_Spvn_dTdtau();
    }
}


void PhotonEmission::outputPhotonSpvn_individualchannel() {
    // dilepton_QGP_thermal->outputPhoton_SpvnpT_shell(output_path);
    if (differential_flag == 1) {
        dilepton_QGP_thermal->outputPhoton_Spectra_dTdtau(output_path, T_cuthigh, T_cutlow, tau_cut_high, tau_cut_low);
        dilepton_QGP_thermal->outputPhoton_Spvn_dTdtau(output_path, T_cuthigh, T_cutlow, tau_cut_high, tau_cut_low);
    }
}


void PhotonEmission::calPhoton_total_Spvn() {
    // integrate over rapidity, azimuthal angle of momentum, and pT
    // calculate dN/(2pi dydM)

    double dy = dilepton_QGP_thermal->get_dy();
    for (int m = 0; m < nm; m++) {
        for (int i = 0; i < np; i++) {
            double p = dilepton_QGP_thermal->getPhotonp(i); //p_T
            double pweight = dilepton_QGP_thermal->getPhoton_pweight(i);
            for (int j = 0; j < nphi; j++) {
                double phi = dilepton_QGP_thermal->getPhotonphi(j); //phi
                double phi_weight = dilepton_QGP_thermal->getPhoton_phiweight(j);
                for (int k = 0; k < nrapidity; k++) {
                    double y_weight = dilepton_QGP_thermal->getPhoton_yweight(k); //y
                    // below dy/Dy to make everything rapidity density, Dy is the rapidity range
                    double weight = phi_weight*y_weight*dy/Dy;
                    // integrate over rapidity, azimuthal angle of momentum


                    dNd2pTd2M_eq[m][i]  += dNd2pTdphidy_eq[m][i][j][k]*weight;
                    dNd2pTd2M_eqT[m][i] += dNd2pTdphidy_eqT[m][i][j][k]*weight;
                    dNd2pTd2M_eqL[m][i] += dNd2pTdphidy_eqL[m][i][j][k]*weight;
                    dNd2pTd2M_visc[m][i] += dNd2pTdphidy_visc[m][i][j][k]*weight;
                    dNd2pTd2M_diff[m][i] += dNd2pTdphidy_diff[m][i][j][k]*weight;
                    dNd2pTd2M_tot[m][i] += dNd2pTdphidy_tot[m][i][j][k]*weight;
                    for (int order = 0; order < norder; order++) {
                        vnpT_cos_eq[order][m][i]  += (
                                dNd2pTdphidy_eq[m][i][j][k]*weight*cos(order*phi));
                        vnpT_cos_visc[order][m][i]  += (
                                dNd2pTdphidy_visc[m][i][j][k]*weight*cos(order*phi));
                        vnpT_cos_diff[order][m][i]  += (
                                dNd2pTdphidy_diff[m][i][j][k]*weight*cos(order*phi));
                        vnpT_cos_tot[order][m][i] += (
                                dNd2pTdphidy_tot[m][i][j][k]*weight*cos(order*phi));
                        vnpT_sin_eq[order][m][i]  += (
                                dNd2pTdphidy_eq[m][i][j][k]*weight*sin(order*phi));
                        vnpT_sin_visc[order][m][i]  += (
                                dNd2pTdphidy_visc[m][i][j][k]*weight*sin(order*phi));
                        vnpT_sin_diff[order][m][i]  += (
                                dNd2pTdphidy_diff[m][i][j][k]*weight*sin(order*phi));
                        vnpT_sin_tot[order][m][i] += (
                                dNd2pTdphidy_tot[m][i][j][k]*weight*sin(order*phi));
                    }
                }
            }

            for (int order=0; order < norder; order++) {
                // pT integrated flow, numerator
                vn_cos_eq[order][m]  += vnpT_cos_eq[order][m][i]*p*pweight;
                vn_sin_eq[order][m]  += vnpT_sin_eq[order][m][i]*p*pweight;
                vn_cos_visc[order][m]  += vnpT_cos_visc[order][m][i]*p*pweight;
                vn_sin_visc[order][m]  += vnpT_sin_visc[order][m][i]*p*pweight;
                vn_cos_diff[order][m]  += vnpT_cos_diff[order][m][i]*p*pweight;
                vn_sin_diff[order][m]  += vnpT_sin_diff[order][m][i]*p*pweight;
                vn_cos_tot[order][m] += vnpT_cos_tot[order][m][i]*p*pweight;
                vn_sin_tot[order][m] += vnpT_sin_tot[order][m][i]*p*pweight;

                // pT differential flow
                vnpT_cos_eq[order][m][i]  = vnpT_cos_eq[order][m][i]/dNd2pTd2M_eq[m][i];
                vnpT_cos_visc[order][m][i]  = vnpT_cos_visc[order][m][i]/dNd2pTd2M_visc[m][i];
                vnpT_cos_diff[order][m][i]  = vnpT_cos_diff[order][m][i]/dNd2pTd2M_diff[m][i];
                vnpT_cos_tot[order][m][i] = vnpT_cos_tot[order][m][i]/dNd2pTd2M_tot[m][i];
                vnpT_sin_eq[order][m][i]  = vnpT_sin_eq[order][m][i]/dNd2pTd2M_eq[m][i];
                vnpT_sin_visc[order][m][i]  = vnpT_sin_visc[order][m][i]/dNd2pTd2M_visc[m][i];
                vnpT_sin_diff[order][m][i]  = vnpT_sin_diff[order][m][i]/dNd2pTd2M_diff[m][i];
                vnpT_sin_tot[order][m][i] = vnpT_sin_tot[order][m][i]/dNd2pTd2M_tot[m][i];
            }

            // pT integrated spectra, dN/(MdMdy)
            dNd2Mdy_eq[m]  += dNd2pTd2M_eq[m][i]*p*pweight;
            dNd2Mdy_eqT[m] += dNd2pTd2M_eqT[m][i]*p*pweight;
            dNd2Mdy_eqL[m] += dNd2pTd2M_eqL[m][i]*p*pweight;
            dNd2Mdy_visc[m] += dNd2pTd2M_visc[m][i]*p*pweight;
            dNd2Mdy_diff[m] += dNd2pTd2M_diff[m][i]*p*pweight;
            dNd2Mdy_tot[m] += dNd2pTd2M_tot[m][i]*p*pweight;

            // pT differential spectra, dN/(2pi pTdpT MdM dy)
            dNd2pTd2M_eq[m][i]  = dNd2pTd2M_eq[m][i]/(2*M_PI); 
            dNd2pTd2M_eqT[m][i] = dNd2pTd2M_eqT[m][i]/(2*M_PI);
            dNd2pTd2M_eqL[m][i] = dNd2pTd2M_eqL[m][i]/(2*M_PI);
            dNd2pTd2M_visc[m][i] = dNd2pTd2M_visc[m][i]/(2*M_PI);
            dNd2pTd2M_diff[m][i] = dNd2pTd2M_diff[m][i]/(2*M_PI);
            dNd2pTd2M_tot[m][i] = dNd2pTd2M_tot[m][i]/(2*M_PI);
        }
    }
    
    // pT integrated flow
    for (int m = 0; m < nm; m++) {
        for (int order = 1; order < norder ; order++) {
            vn_cos_eq[order][m]  = vn_cos_eq[order][m]/dNd2Mdy_eq[m];
            vn_sin_eq[order][m]  = vn_sin_eq[order][m]/dNd2Mdy_eq[m];
            vn_cos_visc[order][m]  = vn_cos_visc[order][m]/dNd2Mdy_visc[m];
            vn_sin_visc[order][m]  = vn_sin_visc[order][m]/dNd2Mdy_visc[m];
            vn_cos_diff[order][m]  = vn_cos_diff[order][m]/dNd2Mdy_diff[m];
            vn_sin_diff[order][m]  = vn_sin_diff[order][m]/dNd2Mdy_diff[m];
            vn_cos_tot[order][m] = vn_cos_tot[order][m]/dNd2Mdy_tot[m];
            vn_sin_tot[order][m] = vn_sin_tot[order][m]/dNd2Mdy_tot[m];
        }
    }
}


void PhotonEmission::outputPhoton_total_SpMatrix_and_SpvnpT() {
    ostringstream filename_stream_eq_SpMatrix;
    ostringstream filename_stream_eq_Spvn;
    ostringstream filename_stream_eq_inte_Spvn;

    ostringstream filename_stream_eq_TL_SpMatrix;
    ostringstream filename_stream_eq_TL_Spvn;
    ostringstream filename_stream_eq_TL_inte_Spvn;

    ostringstream filename_stream_visc_SpMatrix;
    ostringstream filename_stream_visc_Spvn;
    ostringstream filename_stream_visc_inte_Spvn;

    ostringstream filename_stream_diff_SpMatrix;
    ostringstream filename_stream_diff_Spvn;
    ostringstream filename_stream_diff_inte_Spvn;

    ostringstream filename_stream_tot_SpMatrix;
    ostringstream filename_stream_tot_Spvn;
    ostringstream filename_stream_tot_inte_Spvn;

    string filename = "Total_dilepton";

    filename_stream_eq_SpMatrix << output_path << filename
                                << "_eq_SpMatrix.dat";
    filename_stream_eq_Spvn << output_path << filename 
                                << "_eq_Spvn.dat";
    filename_stream_eq_inte_Spvn << output_path << filename 
                                << "_eq_Spvn_inte.dat";

    filename_stream_eq_TL_SpMatrix << output_path << filename
                                << "_eq_TL_SpMatrix.dat";
    filename_stream_eq_TL_Spvn << output_path << filename 
                                << "_eq_TL_Sp.dat";
    filename_stream_eq_TL_inte_Spvn << output_path << filename 
                                << "_eq_TL_Sp_inte.dat";
    
    filename_stream_visc_SpMatrix << output_path << filename 
                                << "_visc_SpMatrix.dat";
    filename_stream_visc_Spvn << output_path << filename 
                                << "_visc_Spvn.dat";
    filename_stream_visc_inte_Spvn << output_path << filename 
                                << "_visc_Spvn_inte.dat";

    filename_stream_diff_SpMatrix << output_path << filename 
                                << "_diff_SpMatrix.dat";
    filename_stream_diff_Spvn << output_path << filename 
                                << "_diff_Spvn.dat";
    filename_stream_diff_inte_Spvn << output_path << filename 
                                << "_diff_Spvn_inte.dat";

    filename_stream_tot_SpMatrix << output_path << filename 
                                << "_tot_SpMatrix.dat";
    filename_stream_tot_Spvn << output_path << filename 
                                << "_tot_Spvn.dat";
    filename_stream_tot_inte_Spvn << output_path << filename 
                                << "_tot_Spvn_inte.dat";

    ofstream fphoton_eq_SpMatrix(filename_stream_eq_SpMatrix.str().c_str());
    ofstream fphoton_eq_Spvn(filename_stream_eq_Spvn.str().c_str());
    ofstream fphoton_eq_inte_Spvn(filename_stream_eq_inte_Spvn.str().c_str());

    ofstream fphoton_eq_TL_SpMatrix(filename_stream_eq_TL_SpMatrix.str().c_str());
    ofstream fphoton_eq_TL_Spvn(filename_stream_eq_TL_Spvn.str().c_str());
    ofstream fphoton_eq_TL_inte_Spvn(filename_stream_eq_TL_inte_Spvn.str().c_str());

    ofstream fphoton_visc_SpMatrix(filename_stream_visc_SpMatrix.str().c_str());
    ofstream fphoton_visc_Spvn(filename_stream_visc_Spvn.str().c_str());
    ofstream fphoton_visc_inte_Spvn(filename_stream_visc_inte_Spvn.str().c_str());

    ofstream fphoton_diff_SpMatrix(filename_stream_diff_SpMatrix.str().c_str());
    ofstream fphoton_diff_Spvn(filename_stream_diff_Spvn.str().c_str());
    ofstream fphoton_diff_inte_Spvn(filename_stream_diff_inte_Spvn.str().c_str());

    ofstream fphoton_tot_SpMatrix(filename_stream_tot_SpMatrix.str().c_str());
    ofstream fphoton_tot_Spvn(filename_stream_tot_Spvn.str().c_str());
    ofstream fphoton_tot_inte_Spvn(filename_stream_tot_inte_Spvn.str().c_str());

    double dy = dilepton_QGP_thermal->get_dy();
    for (int m = 0; m < nm; m++) {
        for (int i=0; i < nphi; i++) {
            double phi = dilepton_QGP_thermal->getPhotonphi(i);
            fphoton_eq_SpMatrix << phi << "  ";
            fphoton_eq_TL_SpMatrix << phi << "  ";
            fphoton_visc_SpMatrix << phi << "  ";
            fphoton_diff_SpMatrix << phi << "  ";
            fphoton_tot_SpMatrix << phi << "  ";
            for (int j = 0; j < np; j++) {
                double temp_eq  = 0.0;
                double temp_eqT = 0.0;
                double temp_eqL = 0.0;
                double temp_visc = 0.0;
                double temp_diff = 0.0;
                double temp_tot = 0.0;
                for (int k = 0; k < nrapidity; k++) {
                    double y_weight = dilepton_QGP_thermal->getPhoton_yweight(k);
                    // below dy/Dy to make everything rapidity density, Dy is the rapidity range
                    double weight = y_weight*dy/Dy;
                    temp_eq  += dNd2pTdphidy_eq[m][j][i][k]*weight;
                    temp_eqT += dNd2pTdphidy_eqT[m][j][i][k]*weight;
                    temp_eqL += dNd2pTdphidy_eqL[m][j][i][k]*weight;
                    temp_visc += dNd2pTdphidy_visc[m][j][i][k]*weight;
                    temp_diff += dNd2pTdphidy_diff[m][j][i][k]*weight;
                    temp_tot += dNd2pTdphidy_tot[m][j][i][k]*weight;
                }
                fphoton_eq_SpMatrix << scientific << setprecision(6) << setw(16)
                                    << temp_eq << "  ";
                fphoton_eq_TL_SpMatrix << scientific << setprecision(6) << setw(16)
                                    << temp_eqT << "  " << temp_eqL << "  ";
                fphoton_visc_SpMatrix << scientific << setprecision(6) << setw(16)
                                    << temp_visc << "  ";
                fphoton_diff_SpMatrix << scientific << setprecision(6) << setw(16)
                                    << temp_diff << "  ";
                fphoton_tot_SpMatrix << scientific << setprecision(6) << setw(16)
                                    << temp_tot << "  ";
            }
            fphoton_eq_SpMatrix << endl;
            fphoton_eq_TL_SpMatrix << endl;
            fphoton_visc_SpMatrix << endl;
            fphoton_diff_SpMatrix << endl;
            fphoton_tot_SpMatrix << endl;
        }
    }

    // pT differential, dN/(2pi pTdpT MdM dy) and vn(M, pT)
    for (int m = 0; m < nm; m++) {
        for (int i = 0; i < np; i++) {
            double M_ll = dilepton_QGP_thermal->getDileptonMass(m);
            double pT = dilepton_QGP_thermal->getPhotonp(i);
            fphoton_eq_Spvn << scientific << setprecision(6) << setw(16)
                            << M_ll << "  "<< pT << "  " << dNd2pTd2M_eq[m][i] << "  ";
            fphoton_eq_TL_Spvn << scientific << setprecision(6) << setw(16)
                            << M_ll << "  "<< pT << "  " << dNd2pTd2M_eqT[m][i] << "  " << dNd2pTd2M_eqL[m][i];
            fphoton_visc_Spvn << scientific << setprecision(6) << setw(16)
                            << M_ll << "  "<< pT << "  " << dNd2pTd2M_visc[m][i] << "  ";
            fphoton_diff_Spvn << scientific << setprecision(6) << setw(16)
                            << M_ll << "  "<< pT << "  " << dNd2pTd2M_diff[m][i] << "  ";
            fphoton_tot_Spvn << scientific << setprecision(6) << setw(16)
                            << M_ll << "  "<< pT << "  " << dNd2pTd2M_tot[m][i] << "  ";

            for (int order=1; order < norder; order++) {
                fphoton_eq_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vnpT_cos_eq[order][m][i] << "  "
                                << vnpT_sin_eq[order][m][i] << "  "
                                << sqrt(pow(vnpT_cos_eq[order][m][i], 2)
                                        + pow(vnpT_sin_eq[order][m][i], 2)) << "  ";
                fphoton_visc_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vnpT_cos_visc[order][m][i] << "  "
                                << vnpT_sin_visc[order][m][i] << "  "
                                << sqrt(pow(vnpT_cos_visc[order][m][i], 2)
                                        + pow(vnpT_sin_visc[order][m][i], 2)) << "  ";
                fphoton_diff_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vnpT_cos_diff[order][m][i] << "  "
                                << vnpT_sin_diff[order][m][i] << "  "
                                << sqrt(pow(vnpT_cos_diff[order][m][i], 2)
                                        + pow(vnpT_sin_diff[order][m][i], 2)) << "  ";
                fphoton_tot_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vnpT_cos_tot[order][m][i] << "  "
                                << vnpT_sin_tot[order][m][i] << "  "
                                << sqrt(pow(vnpT_cos_tot[order][m][i], 2)
                                        + pow(vnpT_sin_tot[order][m][i], 2)) << "  ";
            }
            fphoton_eq_Spvn << endl;
            fphoton_eq_TL_Spvn << endl;
            fphoton_visc_Spvn << endl;
            fphoton_diff_Spvn << endl;
            fphoton_tot_Spvn << endl;
        }
    }

    // pT integrated
    for (int m = 0; m < nm; m++) {
        double M_ll = dilepton_QGP_thermal->getDileptonMass(m);
        // get dN/dMdy from dN/MdMdy
        fphoton_eq_inte_Spvn << scientific << setprecision(6) << setw(16)
                        << M_ll << "  " << dNd2Mdy_eq[m]*M_ll << "  "; 
        fphoton_eq_TL_inte_Spvn << scientific << setprecision(6) << setw(16)
                        << M_ll << "  " << dNd2Mdy_eqT[m]*M_ll << "  " << dNd2Mdy_eqL[m]*M_ll;
        fphoton_visc_inte_Spvn << scientific << setprecision(6) << setw(16)
                        << M_ll << "  " << dNd2Mdy_visc[m]*M_ll << "  ";
        fphoton_diff_inte_Spvn << scientific << setprecision(6) << setw(16)
                        << M_ll << "  " << dNd2Mdy_diff[m]*M_ll << "  ";
        fphoton_tot_inte_Spvn << scientific << setprecision(6) << setw(16)
                        << M_ll << "  " << dNd2Mdy_tot[m]*M_ll << "  ";

        for (int order = 0; order < norder; order++) {
            fphoton_eq_inte_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vn_cos_eq[order][m] << "   "
                                << vn_sin_eq[order][m] << "   "
                                << sqrt(pow(vn_cos_eq[order][m], 2)
                                        + pow(vn_sin_eq[order][m], 2)) << "  ";
            fphoton_visc_inte_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vn_cos_visc[order][m] << "   "
                                << vn_sin_visc[order][m] << "   "
                                << sqrt(pow(vn_cos_visc[order][m], 2)
                                    + pow(vn_sin_visc[order][m], 2)) << "  ";
            fphoton_diff_inte_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vn_cos_diff[order][m] << "   "
                                << vn_sin_diff[order][m] << "   "
                                << sqrt(pow(vn_cos_diff[order][m], 2)
                                    + pow(vn_sin_diff[order][m], 2)) << "  ";
            fphoton_tot_inte_Spvn << scientific << setprecision(6) << setw(16)
                                << order << "   " << vn_cos_tot[order][m] << "   "
                                << vn_sin_tot[order][m] << "   "
                                << sqrt(pow(vn_cos_tot[order][m], 2)
                                    + pow(vn_sin_tot[order][m], 2)) << "  ";
        }
        fphoton_eq_inte_Spvn << endl;
        fphoton_eq_TL_inte_Spvn << endl;
        fphoton_visc_inte_Spvn << endl;
        fphoton_diff_inte_Spvn << endl;
        fphoton_tot_inte_Spvn << endl;
    }
}
