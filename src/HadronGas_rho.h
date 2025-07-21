#ifndef SRC_HadronGas_rho_ININ17
#define SRC_HadronGas_rho_ININ17

#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "ThermalPhoton.h"

class HadronGas_rho: public ThermalPhoton {
public:
  HadronGas_rho(std::shared_ptr<ParameterReader> paraRdr_in,
          std::string emissionProcess);
  ~HadronGas_rho();
  std::string ratePath_;
  std::string eosPath_;

  void readInEmissionTables(std::string emissionProcess);
  void interp(const double T, const double M, const double K, double &resRho,double &resRhoL,double &resRhoT);
  void interp_eos(const double T, double& fugacity_pi, double& fugacity_k);

  int getIdx(double xval, std::vector<double> &xTable);
  void getRateFromTable(const double E,
    const double T_local, const double k_local, const double M_local, double &rateTot,
    double &rateT, double &rateL);

  double nF(double x);
  
  double nB(double x); 
  
  std::vector<double> T_list;
  std::vector<double> M_list;
  std::vector<double> K_list;

  double ***rateRho;
  double ***rateRhoL;
  double ***rateRhoT;
  
  double **eos_table;

  static const int nTemp = 10;
  static const int nM = 75;
  static const int nK = 40;
  
  int nTemp_EOS;
  int nele_EOS;

  std::vector<double> T_list_EOS;

  std::string eosfilename;
  int  index_pi ;
  int  index_k ;
  int EOS_table_flag;  
 
};

#endif
