/////////////////////////////////////////////////////////////////////////
//                      hydrodynamics analysis
//                          photon emission
//
//              author: Chun Shen <shen@mps.ohio-state.edu>
//              copyright: Chun Shen
//
//  This program calculates the photon emission from the relativistic
//  heavy ion collision.
//
//  To do in the future:
/////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <memory>
//#include <omp.h>
#ifndef _OPENMP
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif

#include "./PhotonEmission.h"
#include "./Hydroinfo_h5.h"
#include "./Hydroinfo_MUSIC.h"
#include "./Stopwatch.h"
#include "./Arsenal.h"
#include "./ParameterReader.h"
#include "./gauss_quadrature.h"

using namespace std;

int main(int argc, char** argv) {
    Stopwatch sw;

    sw.tic();
    cout << "----------------------------------------" << endl;
#ifdef _OPENMP   
    double start_time, end_time;
    start_time = omp_get_wtime();
    printf("OpenMP acceleration is on...\n");
#endif
    std::shared_ptr<ParameterReader> paraRdr(new ParameterReader());
    paraRdr->readFromFile("parameters.dat");
    paraRdr->readFromArguments(argc, argv);

    // create integration grid along eta direction for boost-invariant medium
    int neta = paraRdr->getVal("neta");
    double eta_i = paraRdr->getVal("eta_i");
    double eta_f = paraRdr->getVal("eta_f");
    bool flag_hydro = paraRdr->getVal("flag_hydro");
    bool flag_prehydro = paraRdr->getVal("flag_prehydro");
    
    double* eta_ptr = new double[neta];
    double* etaweight_ptr = new double[neta];
    gauss_quadrature(neta, 1, 0.0, 0.0, eta_i, eta_f, eta_ptr, etaweight_ptr);

    PhotonEmission thermalPhotons(paraRdr);
    

    // initialize hydro medium
    int hydro_flag = paraRdr->getVal("hydro_flag");
    if(flag_hydro){

    
    if (hydro_flag == 2 ) {
        
        Hydroinfo_MUSIC* hydroinfo_ptr = new Hydroinfo_MUSIC();
	    int hydro_mode = 12;
	    int nskip_tau = 1;
        hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau);
        // calculate thermal photons from the hydro medium
        thermalPhotons.calPhotonemission_3d(hydroinfo_ptr);
        delete hydroinfo_ptr;

         // sum up all channels and compute thermal photon spectra and vn
        thermalPhotons.calPhoton_SpvnpT_individualchannel();
        thermalPhotons.calPhoton_total_Spvn();

         // output results
        thermalPhotons.outputPhotonSpvn_individualchannel();
        thermalPhotons.outputPhoton_total_SpMatrix_and_SpvnpT(hydro_mode);


    } else {
        cout << "main: unrecognized hydro_flag = " << hydro_flag << endl;
        exit(1);
    }

   

    }
    if( flag_prehydro ){

        if (hydro_flag == 2 ) {

         int hydro_mode = 22;
         int nskip_tau = 1;
	 //double sig_lambda_array[3] = {1.249,0.833,0.0};
	 double sig_lambda_array[1] = {0.0};
	 //double sig_lambda_array[1] = {1.249};
	 for(int isig = 0; isig < 3; isig++){
         for(int isuppress_order = 0; isuppress_order < 3; isuppress_order++)
	 {

	  paraRdr->setVal("sig_lambda",sig_lambda_array[isig]);
          paraRdr->setVal("suppress_order",isuppress_order);

         PhotonEmission thermalPhotons_prehydro(paraRdr); 
         Hydroinfo_MUSIC* hydroinfo_ptr_prehydro = new Hydroinfo_MUSIC(); 
         hydroinfo_ptr_prehydro->readHydroData(hydro_mode, nskip_tau);
         thermalPhotons_prehydro.calPhotonemission_3d(hydroinfo_ptr_prehydro,hydro_mode);
         delete hydroinfo_ptr_prehydro;

        // sum up all channels and compute thermal photon spectra and vn
        thermalPhotons_prehydro.calPhoton_SpvnpT_individualchannel();
        thermalPhotons_prehydro.calPhoton_total_Spvn();

         // output results
        thermalPhotons_prehydro.outputPhotonSpvn_individualchannel();
        thermalPhotons_prehydro.outputPhoton_total_SpMatrix_and_SpvnpT(hydro_mode);

        if(flag_hydro){
            thermalPhotons_prehydro.calPhoton_total_Spvn_sum(thermalPhotons);
            thermalPhotons_prehydro.outputPhoton_total_SpMatrix_and_SpvnpT();
        }



        }
        }
	}

    }





#ifdef _OPENMP
    end_time = omp_get_wtime();
    double total_time = end_time - start_time;
    printf("Total wall time: %f seconds.\n", total_time);
#endif

    sw.toc();
    cout << "Total CPU time: " << sw.takeTime() << " seconds. Bye!" << endl;

    // clean up
    delete [] eta_ptr;
    delete [] etaweight_ptr;

    return(0);
}

