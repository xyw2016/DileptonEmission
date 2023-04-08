#include "QGP_LO.h"
#include "data_struct.h"
#include <gsl/gsl_sf_fermi_dirac.h>
#include <cmath>

using PhysConsts::hbarC;
using PhysConsts::alphaEM;
using PhysConsts::eps;

QGP_LO::QGP_LO(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}

// for diffusion rate
double QGP_LO::Bfun(double x){
    if(x>0.25){
        return 0;
    }
    else{
        return (1.+2.*x)*sqrt(1.-4.*x);
    }
}

double QGP_LO::cross_sec(double omega, double q, double qsq, double m_ell2){
    double C_EM = 2./3.;
    double Nc = 3.;

    return Nc*C_EM*(4.* M_PI/3.)*(pow(alphaEM,2)/qsq)*Bfun(m_ell2/qsq);
}

// for a_1 part of diff rate
double QGP_LO::integrand_J(double x,int n,double a,double z) { // define the integrand function
    return pow(x,n)/((exp(x-z)+1.)*(exp(a-x+z)+1.)*(exp(z-x)+1.));
}

double QGP_LO::gaussLegendre(double* xs, double* ws, int m, int n, double a, double b, double z) {
    double integral = 0.0;
    double c = a / 2.0;
    double d = b / 2.0;
    for (int i = 0; i < m; i++) {
        double t = d * xs[i] + c;
        integral += ws[i] * integrand_J(t,n,a,z);
    }
    return d * integral;
}

double QGP_LO::intJn(int n, double a, double b, double z) { // computing J_n(a,b,z) in eq.19; a=q^0/T, b=|q|/T, z=muB/3T
    const int m = 31; // m=31 points gauss integral
    double xs[m] = {9.97087482e-01,  9.84685910e-01,  9.62503925e-01,  9.30756998e-01,
        8.89760030e-01,  8.39920320e-01,  7.81733148e-01,  7.15776785e-01,
        6.42706723e-01,  5.63249161e-01,  4.78193782e-01,  3.88385902e-01,
        2.94718070e-01,  1.98121199e-01,  9.95553122e-02, -6.76318781e-18,
       -9.95553122e-02, -1.98121199e-01, -2.94718070e-01, -3.88385902e-01,
       -4.78193782e-01, -5.63249161e-01, -6.42706723e-01, -7.15776785e-01,
       -7.81733148e-01, -8.39920320e-01, -8.89760030e-01, -9.30756998e-01,
       -9.62503925e-01, -9.84685910e-01, -9.97087482e-01};
    double ws[m] = {0.00747083, 0.01731862, 0.02700902, 0.03643227, 0.04549371,
       0.05410308, 0.06217479, 0.06962858, 0.07639039, 0.08239299,
       0.08757674, 0.09189011, 0.09529024, 0.09774334, 0.09922501,
       0.09972054, 0.09922501, 0.09774334, 0.09529024, 0.09189011,
       0.08757674, 0.08239299, 0.07639039, 0.06962858, 0.06217479,
       0.05410308, 0.04549371, 0.03643227, 0.02700902, 0.01731862,
       0.00747083}; // x and weights for the gaussian integral

    double result = gaussLegendre(xs, ws, m, n, a, b,z);
    return result;
}

double QGP_LO::a1(double omega,double q,double qsq,double T,double muB,double m_ell2,double nB_o_epp){
    // omega = q^0, q = magnitude of 3-vec q, T = temperature, muB = Baryon chemical potential, m_ell2 = the lepton mass squared
    // nB_o_epp = n_B/(e+p) must be import from hydro
    // final rate for a1 is d4R/d4q = a1*q.V/(T^2*kappa-hat) where V is the diffusion current
    // note that the factor T*qsq/(4.*pow((2.*M_PI),5)*pow(q,3)) is not included in a1 here.
    
    double muq = muB/3.;
    double sigma = cross_sec(omega,q,qsq,m_ell2);
    //double ps1 = T*qsq*sigma/(4.*pow((2.*M_PI),5)*pow(q,3));
    double ps2 = 2*omega*nB_o_epp*pow(T,3)*intJn(1,omega/T,q/T,muq/T) - pow(T,2)*(qsq*nB_o_epp + 2.*omega/3.)*intJn(0,omega/T,q/T,muq/T) 
                    + qsq*T*intJn(-1,omega/T,q/T,muq/T)/3.;
    double ps3 = 2*omega*nB_o_epp*pow(T,3)*intJn(1,omega/T,q/T,-muq/T) - pow(T,2)*(qsq*nB_o_epp - 2.*omega/3.)*intJn(0,omega/T,q/T,-muq/T) 
                    - qsq*T*intJn(-1,omega/T,q/T,-muq/T)/3.;

    return sigma*(ps2 + ps3);
}

double QGP_LO::heaviside(double x) {
    if (x < 0.0) {
        return 0.0;
    } else if (x == 0.0) {
        return 0.5;
    } else {
        return 1.0;
    }
}

// double QGP_LO::l1f(double x){
//     return log(1. + exp(-x));
// }

double QGP_LO::l1f(double x) { return +gsl_sf_fermi_dirac_0(-x); }

double QGP_LO::nB(double x){
    return 1./(exp(x) - 1.);
}


// for finite muB rate
double QGP_LO::rhoV(double omega,double k,double ksq, double T,double muB){
    // k is the magnitude of 3-vec k

    double kplus = (omega + k)*0.5;
    double kminus = (omega - k)*0.5;
    double muq = muB/3.;
    double Nc = 3.;
    double ps1 = Nc*ksq/(4.*M_PI*k);
    double ps2 = l1f((kplus - muq)/T) - l1f((abs(kminus) - muq)/T) + l1f((kplus + muq)/T) - l1f((abs(kminus) + muq)/T);
    double ps3 = k*heaviside(kminus);

    return ps1*(T*ps2 + ps3);
}

double QGP_LO::fmuB_rate(double omega,double q,double qsq,double T,double muB,double m_ell2){
    // omega = q^0, q = magnitude of 3-vec q, T = temperature, muB = Baryon chemical potential, m_ell2 = the lepton mass squared
    double C_EM = 2./3.;

    double rate = C_EM*pow(alphaEM,2)*nB(omega/T)*Bfun(m_ell2/qsq)*rhoV(omega,q,qsq,T,muB)/(3.*pow(M_PI,3)*qsq);

    return rate;
}


void QGP_LO::FiniteBaryonRates(double T, double muB, double rhoB_over_eplusp, double Eq, 
    double M_ll, double &eqrate_ptr, double &diffrate_ptr, int include_diff_deltaf){

    double me = 0.0051;
    double m_ell2 = me * me;
    double prefac = 1./pow(hbarC, 4);

    double M2 = M_ll*M_ll;
    double p = sqrt(Eq*Eq - M2);// magnitude of 3-vec p in LRF

    // equilibrium rate
    eqrate_ptr= prefac*fmuB_rate(Eq,p,M2,T,muB,m_ell2);

    // diffusion correction
    if(include_diff_deltaf==1)
        diffrate_ptr = a1(Eq, p, M2, T, muB, m_ell2, rhoB_over_eplusp);
    else
        diffrate_ptr = 0.0;

}
