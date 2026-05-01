
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
using PhysConsts::me;


using ARSENAL::createA3DMatrix;
using ARSENAL::deleteA3DMatrix;
using ARSENAL::createA2DMatrix;
using ARSENAL::deleteA2DMatrix;

HadronGas_rho_omega_phi::HadronGas_rho_omega_phi(std::shared_ptr<ParameterReader> paraRdr_in,
  std::string emissionProcess)
: ThermalPhoton{paraRdr_in, emissionProcess} {

  ratePath_ = "ph_rates/";
  eosPath_ = "EOS/";
  EOS_table_flag = static_cast<int>(paraRdr_in->getVal("EOS_table_flag", 0));

  if (EOS_table_flag == 0){
    eosfilename = "eos_HadronGas_rho_eqrate_200.dat";
    nTemp_EOS = 7;
    nele_EOS = 3;
    index_pi = 1;
    index_k = 2;
  }
  else if (EOS_table_flag == 1)
  {
    eosfilename="eos_HadronGas_rho_eqrate_17p3.dat";
    nTemp_EOS = 21;
    nele_EOS = 13;
    index_pi = 6;
    index_k = 7;
  }


  readInEmissionTables(emissionProcess);
  

  

}
HadronGas_rho_omega_phi::~HadronGas_rho_omega_phi(){
  deleteA3DMatrix(rateRho, nTemp, nK);
  deleteA3DMatrix(rateRho_omega, nTemp, nK);
  deleteA3DMatrix(rateRho_phi, nTemp, nK);
  deleteA2DMatrix(eos_table,nTemp_EOS);

 
}

double HadronGas_rho_omega_phi::nF(double x) {
  double e = exp(-x);
  return e / (1. + e);
}

double HadronGas_rho_omega_phi::nB(double x) {
  double e = exp(-x);
  return e / (1. - e);
};


void HadronGas_rho_omega_phi::readInEmissionTables(std::string emissionProcess) {
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

  rateRho = createA3DMatrix(nTemp, nK, nM, 0.);
  rateRho_omega = createA3DMatrix(nTemp, nK, nM, 0.);
  rateRho_phi = createA3DMatrix(nTemp, nK, nM, 0.);
  
  double T_dump, M_dump, k_dump;

  for (int iT = 0; iT < nTemp; iT++) {
      for (int ik = 0; ik < nK; ik++) {
          for (int im = 0; im < nM; im++) {
                  fin >> T_dump >> M_dump >> k_dump
                      >> rateRho[iT][ik][im]>> rateRho_omega[iT][ik][im]>> rateRho_phi[iT][ik][im];
          
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


  eos_table = createA2DMatrix(nTemp_EOS, nele_EOS, 0.);

  std::ostringstream eos_filename_stream;
  eos_filename_stream << eosPath_ << eosfilename;
  std::cout << "reading in file: [" << eos_filename_stream.str() << "]"
            << std::endl;

  std::ifstream fin_eos;
  fin_eos.open(eos_filename_stream.str().c_str());
  if (fin_eos.is_open() == false) {
      std::cout << "HadronGas_rho::readEOS error: "
                << "the data file cannot be opened." << std::endl;
      exit(-1);
  }

  std::string line;
  double val;
  int iT = 0;
  while (std::getline(fin_eos, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    for (int iele = 0; iele < nele_EOS; ++iele) {
        iss >> val;
        eos_table[iT][iele] = val;
        
        if (iele == 0)
            T_list_EOS.push_back(eos_table[iT][iele]);
    }
    ++iT;
    if (iT >= nTemp_EOS) break;
}
  fin_eos.close();
  
}


int HadronGas_rho_omega_phi::getIdx(double xval, std::vector<double> &xTable) {
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

void HadronGas_rho_omega_phi::interp_eos(const double T, double& fugacity_pi, double& fugacity_k){

  

  if (T > T_list_EOS[nTemp_EOS-1]) {
      fugacity_pi = 1.0;
      fugacity_k  = 1.0;
      return;
  }

  const int i_eos = getIdx(T, T_list_EOS);

  double a_eos = (T - T_list_EOS[i_eos]) / (T_list_EOS[i_eos + 1] - T_list_EOS[i_eos]);
  
  a_eos = std::max(0.0, std::min(1.0, a_eos));

  double mu_pi_eff = 0.0;
  double mu_k_eff  = 0.0;

  for (int iT = 0; iT < 2; ++iT) {
      double wT = (1.0 - a_eos) * (1 - iT) + a_eos * iT;
      mu_pi_eff += wT * eos_table[i_eos + iT][index_pi];
      mu_k_eff  += wT * eos_table[i_eos + iT][index_k];
  }

  fugacity_pi = exp(mu_pi_eff / T);
  fugacity_k  = exp(mu_k_eff  / T);

  // std::cout << std::fixed << std::setprecision(6) 
  //         << T << " " << mu_pi_eff << " " << mu_k_eff << " "<< fugacity_pi<<" "<<fugacity_k<< std::endl;
}

  




void HadronGas_rho_omega_phi::interp(
  const double T, const double M, const double K, double &resRho, double &resRho_omega, double &resRho_phi) {
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
  resRho_omega = 0;
  resRho_phi = 0;

  for (int iT = 0; iT < 2; iT++) {
      for (int ik = 0; ik < 2; ik++) {
          for (int im = 0; im < 2; im++) {
            double wT = (1.0 - a)*(1 - iT) + a*(iT);
            double wK = (1.0 - b)*(1 - ik) + b*(ik);
            double wM = (1.0 - c)*(1 - im) + c*(im);
            
            double weight = wT * wK * wM;
            resRho += weight * rateRho[i + iT][j + ik][k + im];  
            resRho_omega += weight * rateRho_omega[i + iT][j + ik][k + im];  
            resRho_phi += weight * rateRho_phi[i + iT][j + ik][k + im];  
          }
      }
  }
}

void HadronGas_rho_omega_phi::getRateFromTable(const double E,
  const double T_local, const double k_local, const double M_local, double &rateTot,
  double &rateT, double &rateL) {
  
  double rho_app = 0.0;
  double omega_app = 0.0;
  double phi_app = 0.0;

  interp(T_local,  M_local, k_local,  rho_app,omega_app,phi_app);
  
  double zk = 1;
  double zpi = 1;
  interp_eos(T_local,  zpi, zk);

  //NPA806(2008)339
  double rho_figucity_factor = zpi*zpi;
  double omega_figucity_factor = rho_figucity_factor*zpi;
  double phi_figucity_factor = zk*zk*0.75*0.75;

 
  
  

  double ImDV = rho_app;
  double ImDV_omega = omega_app;
  double ImDV_phi = phi_app;
  double M2 = M_local * M_local;
  double mass_rho0 = 0.853;
  double gv2 = 5.9*5.9;

  double mass_omega= 0.782;
  double mass_phi = 1.02;
  double gv2_omega = 20.5*4*M_PI;
  double gv2_phi = 11.7*4*M_PI;
  //rare mass for omega and phi
  
  double prefactor = - alphaEM*alphaEM/pow(M_PI, 3.)/M2/pow(hbarC, 4.) ;
  double factorLM = (1+2*me*me/M2 )*sqrt(1-4*me*me/M2 );

 
  //rateTot = factorLM*prefactor * nB(E / T_local)*(ImDV*pow(mass_rho0,4)*rho_figucity_factor/gv2
  //                                      +ImDV_omega*pow(mass_omega,4)*omega_figucity_factor/gv2_omega 
  //					                            +ImDV_phi*pow(mass_phi,4)*phi_figucity_factor/gv2_phi);
  rateTot = factorLM*prefactor * nB(E / T_local)*(ImDV*pow(mass_rho0,4)*rho_figucity_factor/gv2);
 
   
  rateT = 2./3. * rateTot;
  rateL = 1./3. * rateTot;
}




