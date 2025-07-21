#ifndef SRC_HadronGas_4PI_V
#define SRC_HadronGas_4PI_V

#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "ThermalPhoton.h"

class HadronGas_4pi: public ThermalPhoton {
public:
  HadronGas_4pi(std::shared_ptr<ParameterReader> paraRdr_in,
          std::string emissionProcess);
  ~HadronGas_4pi();
  std::string eosPath_;

  void computeSFEMcont(double M, double &resRho);

   void getRateFromTable(const double E,
    const double T_local, const double k_local, const double M_local, double &rateTot,
    double &rateT, double &rateL);

  void readInEmissionTables(std::string emissionProcess);
  void interp_eos(const double T, double& fugacity_pi, double& fugacity_k);

  int getIdx(double xval, std::vector<double> &xTable);
  double nB(double x); 
  
  double **eos_table;

  
  int nTemp_EOS;
  int nele_EOS;

  std::vector<double> T_list_EOS;

  std::string eosfilename;
  int  index_pi ;
  int  index_k ;
  int EOS_table_flag;  

  double Evec   ;
  double delvec ;
  double wpi6   ;
  double Nc     ;
  double flav   ;
 
};

#endif
