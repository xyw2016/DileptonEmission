// Adapted by L. Du from H. Gao
#include "QGP_LO.h"
#include "data_struct.h"
#include <gsl/gsl_sf_fermi_dirac.h>
#include <cmath>
#include <iostream>

using PhysConsts::hbarC;
using PhysConsts::alphaEM;
using PhysConsts::me;

#define OOFP 0.0795774715459476678844418816863 //1/(4*pi)
#define sgn(x) (double) ((x>0)-(x<0))

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

    return 4*Nc*C_EM*(4.*M_PI/3.)*(pow(alphaEM,2)/qsq)*Bfun(m_ell2/qsq); // 4=initial spin configurations not to be average. i.e., the x-sec used used be 4-times larger than p&s (5.13)
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

double QGP_LO::a1(double omega,double q,double qsq,double T,double muB,double m_ell2,double sigma,double nB_o_epp){
    // omega = q^0, q = magnitude of 3-vec q, T = temperature, muB = Baryon chemical potential, m_ell2 = the lepton mass squared
    // nB_o_epp = n_B/(e+p) must be import from hydro
    // final rate for a1 is d4R/d4q = a1*q.V/(T^2*kappa-hat) where V is the diffusion current
    // note that the factor T*qsq/(4.*pow((2.*M_PI),5)*pow(q,3)) is not included in a1 here. Moved to PhotonEmission.cpp
    
    double muq = muB/3.;
    // double sigma = cross_sec(omega,q,qsq,m_ell2);
    //double ps1 = T*qsq*sigma/(4.*pow((2.*M_PI),5)*pow(q,3));
    double ps2 = 2*omega*nB_o_epp*pow(T,3)*intJn(1,omega/T,q/T,muq/T) - pow(T,2)*(qsq*nB_o_epp + 2.*omega/3.)*intJn(0,omega/T,q/T,muq/T) 
                    + qsq*T*intJn(-1,omega/T,q/T,muq/T)/3.;
    double ps3 = 2*omega*nB_o_epp*pow(T,3)*intJn(1,omega/T,q/T,-muq/T) - pow(T,2)*(qsq*nB_o_epp - 2.*omega/3.)*intJn(0,omega/T,q/T,-muq/T) 
                    - qsq*T*intJn(-1,omega/T,q/T,-muq/T)/3.;

    return sigma*(ps2 + ps3);
}

double QGP_LO::s2(double omega,double q,double qsq,double T,double muB,double m_ell2,double sigma){
    // omega = q^0, q = magnitude of 3-vec q, qsq = (q^0)^2 - q^2, T = temperature, muB = Baryon chemical potential, m_ell2 = the lepton mass squared
    // The final DPR is given by
    // deltaDPR_shear = q^a q^b pi^ab *s2 /(T^2(E+P))
    // note that the factor Cq*T*qsq/(4.*pow((2.*M_PI),5)*pow(q,5)) is not included in s2 here. Moved to PhotonEmission.cpp

    double muq = muB/3.;

    double j2 = intJn(2,omega/T,q/T,-muq/T) + intJn(2,omega/T,q/T,muq/T);
    double j1 = intJn(1,omega/T,q/T,-muq/T) + intJn(1,omega/T,q/T,muq/T);
    double j0 = intJn(0,omega/T,q/T,-muq/T) + intJn(0,omega/T,q/T,muq/T);

    // double sigma = cross_sec(omega,q,qsq,m_ell2);
    double ps2 = (3.*pow(omega,2) - pow(q,2))*pow(T,2)*j2 - 3.*omega*qsq*T*j1 + 3.*pow(qsq,2)*j0/4.;

    return sigma*ps2;
}

double QGP_LO::b1(double omega,double q, double qsq,double T,double muB,double sigma){
    // omega = q^0, 
    // q = magnitude of 3-vec q, 
    // qsq = Mll^2, 
    // T = temperature
    // muB = Baryon chemical potential,
    
    // return sigma/alpha_em*( 2q0*T* J0(omega/T,q/T,muq/T) - 2q0*T* J0(omega/T,q/T,-muq/T)
    //               - q^2*( intJn(-1,omega/T,q/T,muq/T) - intJn(-1,omega/T,q/T,-muq/T) ) )

    double muq = muB/3.;
    double j0  = 2*omega*T * ( intJn(0,omega/T,q/T,muq/T) - intJn(0,omega/T,q/T,-muq/T) );
    double jm1 = qsq* ( intJn(-1,omega/T,q/T,muq/T) - intJn(-1,omega/T,q/T,-muq/T) );
    
    

    double ps2 = j0 - jm1;

    
    ps2 = ps2/alphaEM;
    double factor_pi = 4.0 * pow(2*M_PI,6);
    return sigma*ps2/factor_pi;
}

double QGP_LO::heaviside(double x) {
    if (x < 0.0) {
        return 0.0;
    } else {
        return 1.0;
    }
}

// double QGP_LO::l1f(double x){
//     return log(1. + exp(-x));
// }

double QGP_LO::l1f(double x) { return +gsl_sf_fermi_dirac_0(-x); }
double QGP_LO::l2f(double x) { return -gsl_sf_fermi_dirac_1(-x); }
double QGP_LO::l3f(double x) { return -gsl_sf_fermi_dirac_2(-x); }

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
    //double ps2 = l1f((kplus - muq)/T) - l1f((abs(kminus) - muq)/T) + l1f((kplus + muq)/T) - l1f((abs(kminus) + muq)/T);
    //double ps3 = k*heaviside(kminus);
    //double ps2 = log((cosh(kplus/T) + cosh(muq/T))/(cosh(kminus/T) + cosh(muq/T)));
    double ps2 = k/T - log((1+exp(-(kminus + muq)/T))/(1+exp(-(kplus + muq)/T))) - log((1+exp(-(kminus - muq)/T))/(1+exp(-(kplus - muq)/T)));
    // return ps1*(T*ps2 + ps3);
    return ps1*(T*ps2);
}


double QGP_LO::rhoL(double o, double k, double K2, double T, double muB) { 
      // leading order result, see (2.4) of 1910.09567
      // (mu dependence not actually published yet...)

    double T2 = (T*T);
    double mu = muB/3.;

    o /= T; k /= T; K2 /= T2; mu /= T; // make them in units of T

    double rL, r00;
    double kp = .5*(o+k), km = .5*fabs(o-k), somk = sgn(o-k);

    double ps1 = l3f(kp+mu) + l3f(kp-mu) - l3f(km+mu) - l3f(km-mu);
    double ps2 = l2f(kp+mu) + l2f(kp-mu) + somk*( l2f(km+mu) + l2f(km-mu) );
    double ps3 = .5*k*k*(somk+1.);
    r00 = -OOFP * ((ps1*2./k + ps2)*6. + ps3);

    rL =  - K2*r00/(k*k);

    return rL*T2;// put back the unit
}

void QGP_LO::fmuB_rate(double omega,double q,double qsq,double T,double muB,double m_ell2, double &rV, double &rL, double &rT){
    // omega = q^0, q = magnitude of 3-vec q, T = temperature, muB = Baryon chemical potential, m_ell2 = the lepton mass squared

    double C_EM = 2./3.;
    double prefacor = C_EM*pow(alphaEM,2)*nB(omega/T)*Bfun(m_ell2/qsq)/(3.*pow(M_PI,3)*qsq);

    rV = prefacor*rhoV(omega,q,qsq,T,muB);
    //rL = prefacor*rhoL(omega,q,qsq,T,muB);    
    //rT = .5*(rV-rL);

    //todo: gsl: fermi_dirac.c:1325: ERROR: underflow
    //Default GSL error handler invoked.
    //wxy:20240508

    rL = 0.0;    
    rT = 0.0;

}


void QGP_LO::FiniteBaryonRates(double T, double muB, double inv_eplusp, double rhoB_over_eplusp, double Eq, 
    double M_ll, double &eqrate_ptr, double &eqrateT_ptr, double &eqrateL_ptr, double &viscrate_ptr, double &diffrate_ptr,double &em_ptr,
    int include_visc_deltaf, int include_diff_deltaf,int include_EM_deltaf){

    double m_ell2 = me * me;
    
    double M2 = M_ll*M_ll;
    double p = sqrt(Eq*Eq - M2);// magnitude of 3-vec p in LRF

    // equilibrium rate
    double rV = 0.0;
    double rL = 0.0;
    double rT = 0.0;
    fmuB_rate(Eq,p,M2,T,muB,m_ell2,rV, rL, rT);

    double prefac = 1./pow(hbarC, 4); // prefac to give the right unit
    eqrate_ptr = prefac*rV;
    eqrateL_ptr= prefac*rL;
    eqrateT_ptr= prefac*rT;

    double sigma = cross_sec(Eq,p,M2,m_ell2);
    // shear correction
    if(include_visc_deltaf==1)
        viscrate_ptr = s2(Eq, p, M2, T, muB, m_ell2, sigma);
    else
        viscrate_ptr = 0.0;

    // diffusion correction
    if(include_diff_deltaf==1)
        diffrate_ptr = a1(Eq, p, M2, T, muB, m_ell2, sigma, rhoB_over_eplusp);
    else
        diffrate_ptr = 0.0;

    if(include_EM_deltaf == 1)
       em_ptr = b1(Eq,p, M2, T, muB, sigma);
    else
       em_ptr = 0.0;



}
