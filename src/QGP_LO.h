#ifndef SRC_QGPLO
#define SRC_QGPLO

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

namespace Dilepton{
class QGP_LO : public ThermalPhoton {
 public:
    QGP_LO(std::shared_ptr<ParameterReader> paraRdr_in,
                        std::string emissionProcess);
    ~QGP_LO() {}

    void FiniteBaryonRates(double T, double muB, double inv_eplusp, double rhoB_over_eplusp, double Eq, 
        double M_ll, double &eqrate_ptr, double &eqrateT_ptr, double &eqrateL_ptr, double &viscrate_ptr, 
        double &diffrate_ptr, int include_visc_deltaf, int include_diff_deltaf);

    double integrand_J(double x,int n,double a,double z);
    double gaussLegendre(double* xs, double* ws, int m, int n, double a, double b, double z);
    double intJn(int n, double a, double b, double z) ;

    // for diffusion rate
    double Bfun(double x);
    double heaviside(double x) ;
    double cross_sec(double omega, double q, double qsq, double m_ell2);
    double a1(double omega,double q, double qsq,double T,double muB,double m_ell2,double sigma,double nB_o_epp);
    double s2(double omega,double q,double qsq,double T,double muB,double m_ell2,double sigma);
    double l1f(double x);
    double l2f(double x);
    double l3f(double x);
    double nB(double x);

    // for finite muB rate
    double rhoV(double omega,double k, double ksq,double T,double muB);
    double rhoL(double o, double k, double K2, double T, double mu);
    void fmuB_rate(double omega,double q, double qsq,double T,double muB,double m_ell2, double &rV, double &rL, double &rT);
};
}
#endif
