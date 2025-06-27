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

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
//#include <omp.h>
#ifndef _OPENMP
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#else
#include <omp.h>
#endif

#include "./Arsenal.h"
#include "./Hydroinfo_MUSIC.h"
#include "./Hydroinfo_h5.h"
#include "./ParameterReader.h"
#include "./PhotonEmission.h"
#include "./Stopwatch.h"
#include "./gauss_quadrature.h"
#include "data_struct.h"

using namespace std;

int main(int argc, char **argv) {
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

  double *eta_ptr = new double[neta];
  double *etaweight_ptr = new double[neta];
  gauss_quadrature(neta, 1, 0.0, 0.0, eta_i, eta_f, eta_ptr, etaweight_ptr);
    
  int dilepton_type = paraRdr->getVal("dilepton_type");
  if (dilepton_type == 0) {
        PhysConsts::me = PhysConsts::me_electron;
  } else if (dilepton_type == 1) {
        PhysConsts::me = PhysConsts::me_muon;
  } else {
        std::cerr << "Invalid input." << std::endl;
        return 1;
  }
  std::cout << "Lepton mass: " << PhysConsts::me << " GeV\n";
  PhotonEmission thermalPhotons(paraRdr);

  // initialize hydro medium
  int hydro_flag = paraRdr->getVal("hydro_flag");
  bool USE_2D_mode = paraRdr->getVal("USE_2D_mode");
  int hydro_info_tz_flag = paraRdr->getVal("hydro_info_tz");
  if (flag_hydro) {

    if (hydro_flag == 2) {

      Hydroinfo_MUSIC *hydroinfo_ptr = new Hydroinfo_MUSIC();
      int hydro_mode = 12;
      int nskip_tau = 1;
      hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau,hydro_info_tz_flag);
      //     // calculate thermal photons from the hydro medium
      if(USE_2D_mode){
          thermalPhotons.calPhotonemission_2d(hydroinfo_ptr,hydro_mode);
      }
      else{
          thermalPhotons.calPhotonemission_3d(hydroinfo_ptr, hydro_mode);
      }
      delete hydroinfo_ptr;

      //      // sum up all channels and compute thermal photon spectra and vn
          thermalPhotons.calPhoton_SpvnpT_individualchannel();
          thermalPhotons.calPhoton_total_Spvn();

      //      // output results
          thermalPhotons.outputPhotonSpvn_individualchannel("hydro");
          thermalPhotons.outputPhoton_total_SpMatrix_and_SpvnpT("hydro");
    }
    else {
        cout << "main: unrecognized hydro_flag = " << hydro_flag << endl;
        exit(1);
    }
  }
  if( flag_prehydro ){

      if (hydro_flag == 2 ) {

       int hydro_mode = 22;
       int nskip_tau = 1;
       double sig_lambda_array[1] = {0.0};
       //double sig_lambda_array[1] = {1.249};
       for(int isig = 0; isig < 1; isig++){
           for(int isuppress_order = 0; isuppress_order < 3; isuppress_order++)
       {

        paraRdr->setVal("sig_lambda",sig_lambda_array[isig]);
        paraRdr->setVal("suppress_order",isuppress_order);
        std::cout<<" ============================= " <<std::endl;
        std::cout<<" sig_lambda " << sig_lambda_array[isig]<<std::endl;
        std::cout<<" suppress_order " << isuppress_order<<std::endl;
        std::cout<<" ============================= " <<std::endl;
        PhotonEmission thermalPhotons_prehydro(paraRdr);
        Hydroinfo_MUSIC* hydroinfo_ptr_prehydro = new Hydroinfo_MUSIC();
        hydroinfo_ptr_prehydro->readHydroData(hydro_mode, nskip_tau);
        if(USE_2D_mode){
          thermalPhotons_prehydro.calPhotonemission_2d(hydroinfo_ptr_prehydro,hydro_mode);
        }
        else{
          thermalPhotons_prehydro.calPhotonemission_3d(hydroinfo_ptr_prehydro,hydro_mode);
        }
        delete hydroinfo_ptr_prehydro;

      // sum up all channels and compute thermal photon spectra and vn
        thermalPhotons_prehydro.calPhoton_SpvnpT_individualchannel();
        thermalPhotons_prehydro.calPhoton_total_Spvn();
        
        ostringstream file_name_label;  
        file_name_label << "pre_hydro" << "_"<<sig_lambda_array[isig]<<"_"<< isuppress_order;
       // output results
        thermalPhotons_prehydro.outputPhotonSpvn_individualchannel(file_name_label.str());
        thermalPhotons_prehydro.outputPhoton_total_SpMatrix_and_SpvnpT(file_name_label.str());

      if(flag_hydro){
          file_name_label.str("");
          file_name_label.clear(); 
          file_name_label << "allstage" << "_"<<sig_lambda_array[isig]<<"_"<< isuppress_order;;
          thermalPhotons_prehydro.calPhoton_total_Spvn_sum(thermalPhotons);
          thermalPhotons_prehydro.outputPhoton_total_SpMatrix_and_SpvnpT(file_name_label.str());
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
  delete[] eta_ptr;
  delete[] etaweight_ptr;

  return (0);
}
