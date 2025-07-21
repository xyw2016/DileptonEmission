
#include "HadronGas_4pi.h"
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

HadronGas_4pi::HadronGas_4pi(std::shared_ptr<ParameterReader> paraRdr_in,
  std::string emissionProcess)
: ThermalPhoton{paraRdr_in, emissionProcess} {

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
  
  Evec   = 1.150;
  delvec = 0.110;
  wpi6   = 6.0 * 0.13957;
  Nc     = 3.0;
  flav   = (4.0/9.0 + 1.0/9.0 + 1.0/9.0);


  readInEmissionTables(emissionProcess);

}

HadronGas_4pi::~HadronGas_4pi(){
  deleteA2DMatrix(eos_table,nTemp_EOS);
 
}



double HadronGas_4pi::nB(double x) {
  double e = exp(-x);
  return e / (1. - e);
};


void HadronGas_4pi::readInEmissionTables(std::string emissionProcess) {
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




int HadronGas_4pi::getIdx(double xval, std::vector<double> &xTable) {
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


void HadronGas_4pi::interp_eos(const double T, double& fugacity_pi, double& fugacity_k){

  

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

}


void HadronGas_4pi::computeSFEMcont(double M, double &resRho){

  double Theta_con = 0.0;
  double contvac   = 0.0;

    if (M > wpi6) {
        Theta_con = sqrt(1.0 - (wpi6 * wpi6) / (M * M));
        contvac   = 1.0 / (1.0 + exp((Evec - M) / delvec));
    }

    resRho= Theta_con * contvac * Nc * flav * M * M / (12.0 * M_PI);
  }

void HadronGas_4pi::getRateFromTable(const double E,
  const double T_local, const double k_local, const double M_local, double &rateTot,
  double &rateT, double &rateL) {
  
  double rho_app = 0.0;
  
  computeSFEMcont( M_local,rho_app );

  double zk = 1;
  double zpi = 1;
  interp_eos(T_local,  zpi, zk);
  //NPA806(2008)339
  double rho_figucity_factor = zpi*zpi*zpi*zpi;
  double ImDV  = rho_app*rho_figucity_factor;
  double M2 = M_local * M_local;
  
  
  double prefactor = alphaEM*alphaEM/pow(M_PI, 3.)/M2/pow(hbarC, 4.) ;
  double factorLM = (1+2*me*me/M2 )*sqrt(1-4*me*me/M2 );

  

  rateTot = factorLM*prefactor * nB(E / T_local)*ImDV;  
  rateT = 2.*rateTot/3.;
  rateL = 1.*rateTot/3.;

  
}




