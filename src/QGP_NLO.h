#ifndef SRC_QGPNLO
#define SRC_QGPNLO

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

class QGP_NLO : public ThermalPhoton {
 public:
    QGP_NLO(std::shared_ptr<ParameterReader> paraRdr_in,
                        std::string emissionProcess);
    ~QGP_NLO() {}

    //void initialize(std::string fname);
    // void FiniteBaryonRates(double T, double muB, double rhoB_over_eplusp, double Eq, 
    //     double M_ll, double &eqrate_ptr, double &diffrate_ptr, int include_diff_deltaf);

    // double integrand_J(double x,int n,double a,double z);
    // double gaussLegendre(double* xs, double* ws, int m, int n, double a, double b, double z);
    // double intJn(int n, double a, double b, double z) ;

    // // for diffusion rate
    // double Bfun(double x);
    // double heaviside(double x) ;
    // double cross_sec(double omega, double q, double qsq, double m_ell2);
    // double a1(double omega,double q, double qsq,double T,double muB,double m_ell2,double nB_o_epp);
    // double l1f(double x);
    // double nB(double x);

    // // for finite muB rate
    // double rhoV(double omega,double k, double ksq,double T,double muB);
    // double fmuB_rate(double omega,double q, double qsq,double T,double muB,double m_ell2);
};

#endif
