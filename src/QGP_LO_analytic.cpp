
#include "QGP_LO_analytic.h"
#include "data_struct.h"
#include <cmath>

using PhysConsts::hbarC;

QGP_LO_analytic::QGP_LO_analytic(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}


// PRC. 93, 044902, 2016
void QGP_LO_analytic::FiniteBaryonRates(double T, double muB, double rhoB_over_eplusp, double Eq, 
    double M_ll, double &eqrate_ptr, double &diffrate_ptr, int include_diff_deltaf) {

    const double aem = 1./137.;
    const double Qu  = 2./3.;
    const double Qd  = -1./3.;
    const double Qs  = -1./3.;

    double prefac = (Qu*Qu+Qd*Qd+Qs*Qs)*aem*aem/(2.*pow(M_PI, 4))/pow(hbarC, 4);

    double M = M_ll;
    double E = Eq;
    double p = sqrt(E*E - M*M);
    double x = E/T;
    double y = p/T;
    double fq = 1./(exp(x) - 1.);
    double log_r = log(cosh(0.25*(x+y))/cosh(0.25*(x-y)));

    eqrate_ptr = prefac/y*fq*log_r;
    diffrate_ptr = 0.0;

}

