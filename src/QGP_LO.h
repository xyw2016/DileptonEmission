#ifndef SRC_QGP2TO2TOTAL
#define SRC_QGP2TO2TOTAL

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

class QGP_LO : public ThermalPhoton {
 public:
    QGP_LO(std::shared_ptr<ParameterReader> paraRdr_in,
                        std::string emissionProcess);
    ~QGP_LO() {}
    void analyticRates(double T, std::vector<double> &Eq, double *M,
                        std::vector<double> &eqrate_ptr, int nm,
                        int np, int nphi, int nrapidity);
    void analyticRatesDiffusion(double T, double muB, double rhoB_over_eplusp,
        std::vector<double> &Eq, double *M_ll, std::vector<double> &diffusion_ptr, 
        int nm, int np, int nphi, int nrapidity);
    // void NetBaryonCorrection(double T, double muB, std::vector<double> &Eq,
    //                          std::vector<double> &eqrate_ptr);
    double integrand_J(double x,int n,double a,double z);

    double gaussLegendre(double* xs, double* ws, int m, int n, double a, double b, double z);

    double intJn(int n, double a, double b, double z) ;

    // for diffusion rate
    double Bfun(double x);

    double heaviside(double x) ;

    double cross_sec(double omega, double q, double m_ell2);

    double a1(double omega,double q,double T,double muB,double m_ell2,double nB_o_epp);



    double l1f(double x);

    double nB(double x);


    // for finite muB rate
    double rhoV(double omega,double k,double T,double muB);

    double fmuB_rate(double omega,double q,double T,double muB,double m_ell2);
};

#endif
