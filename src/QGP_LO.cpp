
#include "QGP_LO.h"
#include "data_struct.h"
#include <cmath>

using PhysConsts::hbarC;

QGP_LO::QGP_LO(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}


// PRC. 93, 044902, 2016
void QGP_LO::analyticRates(
            double T, std::vector<double> &Eq, double *M_ll,
            std::vector<double> &eqrate_ptr, int nm,
            int np, int nphi, int nrapidity) {

    const double aem = 1./137.;
    const double Qu  = 2./3.;
    const double Qd  = -1./3.;
    const double Qs  = -1./3.;

    double prefac = (Qu*Qu+Qd*Qd+Qs*Qs)*aem*aem/(2.*pow(M_PI, 4))/pow(hbarC, 4);

    int idx = 0;
    for (int k = 0; k < nrapidity; k++) {
        for (int m = 0; m < nphi; m++) {
            for (int l = 0; l < np; l++) {
                for (int j = 0; j < nm; j++) {

                    double M = M_ll[j];
                    double E = Eq[idx];
                    double p = sqrt(E*E - M*M);
                    double x = E/T;
                    double y = p/T;
                    double fq = 1./(exp(x) - 1.);
                    double log_r = log(cosh(0.25*(x+y))/cosh(0.25*(x-y)));

                    eqrate_ptr[idx] = prefac/y*fq*log_r;
                    idx++;
                }
            }
        }
    }
}

// void QGP_LO::NetBaryonCorrection(
//         double T, double muB, std::vector<double> &Eq,
//         std::vector<double> &eqrate_ptr) {

//     muB = std::abs(muB);
//     if (muB < 1e-8) return;

//     // set the upper limit to 0.4 GeV for this parameterization
//     muB = std::min(0.4, muB);

//     const double muB_over_T = muB/T;
//     const double kmuB = 1. + 1./(M_PI*M_PI)*muB_over_T*muB_over_T;

//     for (unsigned int i = 0; i < Eq.size(); i++) {
//         eqrate_ptr[i] *= kmuB;
//     }
// }
