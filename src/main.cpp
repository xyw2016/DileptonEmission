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
    double* eta_ptr = new double[neta];
    double* etaweight_ptr = new double[neta];
    gauss_quadrature(neta, 1, 0.0, 0.0, eta_i, eta_f, eta_ptr, etaweight_ptr);

    PhotonEmission thermalPhotons(paraRdr);

    // initialize hydro medium
    int hydro_flag = paraRdr->getVal("hydro_flag");

    if (hydro_flag == 2 || hydro_flag == 22  ) {
        Hydroinfo_MUSIC* hydroinfo_ptr = new Hydroinfo_MUSIC();
	int hydro_mode = 12;
	//if(hydro_flag == 22){
	//    hydro_mode = 13;
	//}
        
	int nskip_tau = 1;
        hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau);
        // calculate thermal photons from the hydro medium
        thermalPhotons.calPhotonemission_3d(hydroinfo_ptr);
        delete hydroinfo_ptr;
    } else {
        cout << "main: unrecognized hydro_flag = " << hydro_flag << endl;
        exit(1);
    }

    // sum up all channels and compute thermal photon spectra and vn
    thermalPhotons.calPhoton_SpvnpT_individualchannel();
    thermalPhotons.calPhoton_total_Spvn();

    // output results
    thermalPhotons.outputPhotonSpvn_individualchannel();
    thermalPhotons.outputPhoton_total_SpMatrix_and_SpvnpT();

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

