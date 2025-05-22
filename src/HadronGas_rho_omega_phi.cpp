
#include "HadronGas_rho_omega_phi.h"
#include "ThermalPhoton.h"
#include "data_struct.h"
#include <algorithm>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>

#include "Arsenal.h"

using namespace std;
using PhysConsts::alphaEM;
using PhysConsts::hbarC;

using ARSENAL::createA3DMatrix;
using ARSENAL::deleteA3DMatrix;

HadronGas_rho::HadronGas_rho(std::shared_ptr<ParameterReader> paraRdr_in,
  std::string emissionProcess)
: ThermalPhoton{paraRdr_in, emissionProcess} {

  ratePath_ = "ph_rates/";
  readInEmissionTables(emissionProcess);
  

}

double HadronGas_rho::nF(double x) {
  double e = exp(-x);
  return e / (1. + e);
}

double HadronGas_rho::nB(double x) {
  double e = exp(-x);
  return e / (1. - e);
};


void HadronGas_rho::readInEmissionTables(std::string emissionProcess) {
  std::ostringstream eqrate_filename_stream;
  eqrate_filename_stream << ratePath_ << "rate_" << emissionProcess
                         << "_eqrate.dat";
  std::cout << "reading in file: [" << eqrate_filename_stream.str() << "]"
            << std::endl;

  std::ifstream fin;
  fin.open(eqrate_filename_stream.str().c_str());
  //std::cout << eqrate_filename_stream.str() << " wwwww" <<endl;
  if (fin.is_open() == false) {
      std::cout << "DileptonQGPNLO::readInEmissionTables error: "
                << "the data file cannot be opened." << std::endl;
      exit(-1);
  }

  const int nTemp = 9;
  const int nM = 89;
  const int nK = 40;
  rateRho = createA3DMatrix(nTemp, nK, nM, 0.);
  
  double T_dump, M_dump, k_dump;

  for (int iT = 0; iT < nTemp; iT++) {
      for (int ik = 0; ik < nK; ik++) {
          for (int im = 0; im < nM; im++) {
                  fin >> T_dump >> M_dump >> k_dump
                      >> rateRho[iT][ik][im];
          
                  if (iT == 0 && ik == 0 )
                      M_list.push_back(M_dump);
                  if (ik == 0 && im == 0)
                      T_list.push_back(T_dump);
                  if (iT == 0 && im == 0 )
                      K_list.push_back(k_dump);
          }
      }
  }

  for (int iT = 0; iT < nTemp; iT++) {
    std::cout<<T_list[iT]<<" ";
  }
  std::cout<< std::endl;

  for (int ik = 0; ik < nK; ik++) {
    std::cout<<K_list[ik]<<" ";
  }
  std::cout<< std::endl;


  for (int im = 0; im < nM; im++) {
    std::cout<<M_list[im]<<" ";
  }
  std::cout<< std::endl;



  std::cout << " ... done!" << std::endl;
  fin.close();
}


int HadronGas_rho::getIdx(double xval, std::vector<double> &xTable) {
  // binary search
  if (xval > xTable[xTable.size() - 1]) return xTable.size() - 2;
  if (xval < xTable[0]) return 0;
  int iU = xTable.size() - 1;
  int iL = 0;
  int iM = (iU + iL) / 2;
  while (iU - iL > 1) {
      if (xval >= xTable[iM])
          iL = iM;
      else
          iU = iM;
      iM = (iU + iL) / 2;
  }
  return iL;
}

void HadronGas_rho::interp(
  const double T, const double M, const double K, double &resRho) {
  const int i = getIdx(T, T_list);
  const int j = getIdx(K, K_list);
  const int k = getIdx(M, M_list);
  
  double a = (T - T_list[i]) / (T_list[i + 1] - T_list[i]);
  double b =
      (K - K_list[j]) / (K_list[j + 1] - K_list[j]);
  double c =
      (M - M_list[k]) / (M_list[k + 1] - M_list[k]);

  // avoid overflows
  a = std::max(0., std::min(1., a));
  b = std::max(0., std::min(1., b));
  c = std::max(0., std::min(1., c));


  resRho = 0;

  for (int iT = 0; iT < 2; iT++) {
      for (int ik = 0; ik < 2; ik++) {
          for (int im = 0; im < 2; im++) {
            double wT = (1.0 - a)*(1 - iT) + a*(iT);
            double wK = (1.0 - b)*(1 - ik) + b*(ik);
            double wM = (1.0 - c)*(1 - im) + c*(im);
            
            double weight = wT * wK * wM;
            resRho += weight * rateRho[i + iT][j + ik][k + im];  
          }
      }
  }
}

void HadronGas_rho::getRateFromTable(const double E,
  const double T_local, const double k_local, const double M_local, double &rateTot,
  double &rateT, double &rateL) {
  
  double rho_app = 0.0;
  interp(T_local,  M_local, k_local,  rho_app);
  // double iMM = 0.45;
  // double iKK = 0.5;

  // double iEE = sqrt(iMM*iMM + iKK*iKK);
  
  
  // interp(T_local,  iMM, iKK,  rho_app);
  

  double ImDV = rho_app;
  double M2 = M_local * M_local;
  //double M2 = iMM * iMM;
  double mass_rho0 = 0.77;
  double gv2 = 2.54*4*M_PI; 
  
  double prefactor = - alphaEM*alphaEM/pow(M_PI, 3.)/M2/pow(hbarC, 4.) ;


  

  rateTot = prefactor * nB(E / T_local)*ImDV*pow(mass_rho0,4)/gv2;
  //rateTot = prefactor * nB(iEE / T_local)*ImDV*pow(mass_rho0,4)/gv2;

  //std::cout<< prefactor<<" "<<nB(E / T_local)<<" "<< ImDV<<" " <<T_local<<" "<<rateTot<<std::endl;
  
  rateT = 1./3. * rateTot;
  rateL = 1./3. * rateTot;
}




