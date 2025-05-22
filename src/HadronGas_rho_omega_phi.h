#ifndef SRC_HadronGas_rho
#define SRC_HadronGas_rho

#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "ThermalPhoton.h"

class HadronGas_rho: public ThermalPhoton {
public:
  HadronGas_rho(std::shared_ptr<ParameterReader> paraRdr_in,
          std::string emissionProcess);
  ~HadronGas_rho() {}
  std::string ratePath_;
  void readInEmissionTables(std::string emissionProcess);
  void interp(const double T, const double M, const double K, double &resRho);
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
  
 
};

#endif
