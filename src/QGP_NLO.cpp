// GJ, 31.03.2023
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "data_struct.h"
#include "ThermalPhoton.h"
#include "QGP_NLO.h"

using namespace std;
using PhysConsts::hbarC;
using PhysConsts::alphaEM;


QGP_NLO::QGP_NLO(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}

/*--------------------------------------------------------------------*/
// some definations

#define sgn(x) (double) ((x>0)-(x<0))
#define sz size()
#define loop(i,a,b) for(int (i)=(a);(i)<(b);(i)++)
#define OOFP 0.0795774715459476678844418816863

double Nc=3.;
double Cf=((Nc)*(Nc)-1)/(2.*Nc);

/*--------------------------------------------------------------------*/
// some pQCD results

double nF(double x)  { double e=exp(-x); return e/(1.+e); };
double nB(double x)  { double e=exp(-x); return e/(1.-e); };
double l1f(double x) { return +gsl_sf_fermi_dirac_0(-x); }
double l2f(double x) { return -gsl_sf_fermi_dirac_1(-x); }
double l3f(double x) { return -gsl_sf_fermi_dirac_2(-x); }

void rho_LO(double o, double k, double mu, double &rT, double &rL) { 
  // leading order result, see (2.4) of 1910.09567
  // (mu dependence not actually published yet...)
  double rV, r00;
  double K2 = o*o - k*k; // note: omega, k are in units of T!
  double kp = .5*(o+k), km = .5*fabs(o-k), somk = sgn(o-k);

  rV = l1f(kp+mu)+l1f(kp-mu)-l1f(km+mu)-l1f(km-mu);
  rV = rV/k + .5*(somk+1.);
  rV *= K2*3.*OOFP;

  r00 = l3f(kp+mu)+l3f(kp-mu)-l3f(km+mu)-l3f(km-mu);
  r00 = r00*2./k + l2f(kp+mu) + l2f(kp-mu) + somk*( l2f(km+mu) + l2f(km-mu) );
  r00 = r00*6.  + .5*k*k*(somk+1.) ; 
  r00 *= -OOFP;

  rL = -K2*r00/(k*k);
  rT = .5*(rV-rL);
}

double rhoV_AMY(double alpha, double mu, double k) { // (1.9) and (1.10) of hep-ph/0111107
  double nf=3.; 
  double minf2, fcol, f2t2;
  minf2 = 4.* M_PI*alpha*(Cf/4.)*( 1. + 1.*mu*mu/(M_PI*M_PI) )  ;
  f2t2  = .041/k - .3615 + 1.01*exp(-1.35*k);
  fcol  = sqrt(1.+nf/6.)*( .548*log(12.28+1./k)*pow(k,-1.5) + .133*k/sqrt(1.+k/16.27) );

  return .5*(Nc/M_PI)*minf2*(1.-nF(k-mu)-nF(k+mu))*( -.5*log(minf2) + .5*log(2.*k) + f2t2 + fcol );
}

void rho_OPE(double o, double k, double alpha, double mu, double &rT, double &rL) {
  // operator product expansion, see appendix D of 1910.07552
  // (mu dependence not actually published yet...)
  double rV, r00;
  double K2 = o*o - k*k; // note: omega, k are in units of T!
  double coeff = 1. + 96.*mu*mu*OOFP*OOFP + 768.*pow(mu*OOFP,4) ;
            // = T^4 + (6/pi^2).T^2.mu^2 + (3/pi^4).mu^4 )
  rT =  OOFP*K2;
  rL =  OOFP*K2;
  rT += alpha*( 4.*OOFP*OOFP*K2 + coeff/(OOFP*OOFP*27.*K2) );
  rL += alpha*( 4.*OOFP*OOFP*K2 + coeff*(o*o+k*k)/(OOFP*OOFP*27.*K2*K2) );
}

/*--------------------------------------------------------------------*/
// grid structure (for look-up table)

int id(double _x, double *x, int nx) { // locate _x (assume x[nx] ordered)

  // Check if _x is within the range of x[] array
  if (_x < x[0] || _x > x[nx-1]) {
    // handle out of bounds error
    return -1;
  }

  int iU=nx-1, iM, iL=0;
  while (iU-iL > 1) {
    iM = (iU+iL) >> 1; // midpoint
    if (_x >= x[iM]) iL=iM;
    else iU=iM;
  }

  // Check if iL and iL+1 are within the range of the x[] array
  if (iL < 0 || iL >= nx-1 || iL+1 < 0 || iL+1 >= nx) {
    // handle out of bounds error
    return -1;
  }

  return iL;
}


ThermalPhoton::Table::Table(std::vector<double>  &_x, std::vector<double>  &_y, std::vector<double>  &_z, std::vector<double>  &_w, std::vector<double>  &_F) 
    : nx(_x.size()), ny(_y.size()), nz(_z.size()), nw(_w.size()), x(&_x[0]), y(&_y[0]), z(&_z[0]), w(&_w[0])
{
  // Implementation of Table constructor...

  F = new double***[nx];
  auto val = &_F[0];
  loop(i,0,nx) {
    F[i] = new double**[ny];
    loop(j,0,ny) {
    F[i][j] = new double*[nz];
      loop(k,0,nz) {
        F[i][j][k] = new double[nw];
        loop(l,0,nw) {
          F[i][j][k][l] = (val++)[0]; // populate grid
        }
      }
    }
  }

  x_min = *min_element(_x.begin(),_x.end());
  x_max = *max_element(_x.begin(),_x.end());
  y_min = *min_element(_y.begin(),_y.end());
  y_max = *max_element(_y.begin(),_y.end());
  z_min = *min_element(_z.begin(),_z.end());
  z_max = *max_element(_z.begin(),_z.end());
  w_min = *min_element(_w.begin(),_w.end());
  w_max = *max_element(_w.begin(),_w.end());
}


double ThermalPhoton::Table::interp(double _x, double _y, double _z, double _w) {
  // Implementation of Table member function interp...

  int i=id(_x,x,nx), j=id(_y,y,ny), k=id(_z,z,nz), l=id(_w,w,nw);
  double res, a, b, c, d;
  a = (_x-x[i])/(x[i+1]-x[i]);
  b = (_y-y[j])/(y[j+1]-y[j]);
  c = (_z-z[k])/(z[k+1]-z[k]);
  d = (_w-w[l])/(w[l+1]-w[l]);

  if (i < 0 || i >= nx-1 || j < 0 || j >= ny-1 
    || k < 0 || k >= nz-1 || l < 0 || l >= nw-1) { 

    printf("Warning: input values are out of the boundaries of the emission rate table!\n");

  }

  // a,b,c,d (each value between 0 and 1)
  res = (1.-a)*(1.-b)*(1.-c)*(1.-d)*F[i][j][k][l] 
      + (a)*(1.-b)*(1.-c)*(1.-d)*F[i+1][j][k][l] 
      + (1.-a)*(b)*(1.-c)*(1.-d)*F[i][j+1][k][l] 
      + (1.-a)*(1.-b)*(c)*(1.-d)*F[i][j][k+1][l] 
      + (1.-a)*(1.-b)*(1.-c)*(d)*F[i][j][k][l+1] 
      + (a)*(b)*(1.-c)*(1.-d)*F[i+1][j+1][k][l] 
      + (a)*(1.-b)*(c)*(1.-d)*F[i+1][j][k+1][l] 
      + (a)*(1.-b)*(1.-c)*(d)*F[i+1][j][k][l+1] 
      + (1.-a)*(b)*(c)*(1.-d)*F[i][j+1][k+1][l] 
      + (1.-a)*(b)*(1.-c)*(d)*F[i][j+1][k][l+1] 
      + (1.-a)*(1.-b)*(c)*(d)*F[i][j][k+1][l+1] 
      + (1.-a)*(b)*(c)*(d)*F[i][j+1][k+1][l+1] 
      + (a)*(1.-b)*(c)*(d)*F[i+1][j][k+1][l+1] 
      + (a)*(b)*(1.-c)*(d)*F[i+1][j+1][k][l+1] 
      + (a)*(b)*(c)*(1.-d)*F[i+1][j+1][k+1][l] 
      + (a)*(b)*(c)*(d)*F[i+1][j+1][k+1][l+1] 
      ; // hypercube
  return res;
}

void ThermalPhoton::initialize(string fname, std::vector<double> &a_list, std::vector<double> &B_list, std::vector<double> &M_list, 
        std::vector<double> &k_list, std::vector<double> &rhoT_list, std::vector<double> &rhoL_list) {

  cout << "\n reading in file: [" << fname << "]" << endl;

  ifstream fin;
  fin.open(fname);
  fin.ignore(256,'\n');

  double alpha, muB, M, k, rhoT, rhoL;
  int count=0;

  fin >> alpha >> muB >> M >> k >> rhoT >> rhoL; count++;
  
  a_list.push_back(alpha); // alpha_s values
  B_list.push_back(muB); // chemical potential (muB/T)
  M_list.push_back(M); // invariant mass (units of T)
  k_list.push_back(k); // 3-momentum (k/T), defined in the local rest frame
  rhoT_list.push_back(rhoT);
  rhoL_list.push_back(rhoL);

  while (true) {
    fin >> alpha >> muB >> M >> k >> rhoT >> rhoL; count++;
    if (fin.eof()) break;

    rhoT_list.push_back(rhoT);
    rhoL_list.push_back(rhoL);

    if (a_list.back()<alpha)  a_list.push_back(alpha);
    if (B_list.back()<muB)    B_list.push_back(muB);
    if (M_list.back()<M)      M_list.push_back(M);
    if (k_list.back()<k)      k_list.push_back(k);

  }

  cout << " ... done!" << endl;
  cout << endl;
  cout << "-> size of lattice ..."  << endl;
  cout << " number of alpha_i points: " << a_list.sz  << endl;
  cout << " number of muB points:     " << B_list.sz  << endl;
  cout << " number of M points:       " << M_list.sz  << endl;
  cout << " number of k points:       " << k_list.sz  << endl;
  cout << " nx.ny.nz.nw = " << rhoT_list.sz << endl;
  cout << endl;
}

/*--------------------------------------------------------------------*/
// interpolation & extrapolation

int approx_rho(void *input, // omega, k, alpha, muB
               struct ThermalPhoton::Table rhoT_list,
               struct ThermalPhoton::Table rhoL_list,
               double &rT, double &rL // out: T & L spectral functions
               ) {
  // int function returns: -1 if below LC (not needed for dileptons)
  //                       0 if interpolated (i.e. input values fall within grid)
  //                       1 if extrapolated (*)
  // 
  // (*) extrapolation will NOT work properly for large alpha_s, large muB, or
  // if k falls outside the grid! In these cases, just the boundary value is used.

  double o      = ((double*)input)[0];
  double k      = ((double*)input)[1];
  double alpha  = ((double*)input)[2];
  double muB    = ((double*)input)[3];

  double mu = muB/3.; // NB baryon vs quark chemical potential!

  if (o<k) { rT = 0.; rL = 0.; return -1; } // skip 'DIS region'

  double M2 = o*o - k*k;
  double M = sqrt(M2);
  double al_min = rhoT_list.x_min, al_max = rhoT_list.x_max;
  double M_min  = rhoT_list.z_min, M_max  = rhoT_list.z_max;
  double k_min  = rhoT_list.w_min, k_max  = rhoT_list.w_max;

  if (alpha<al_min-1e-5) { // use LO result to extrapolate if alpha < alpha_min
    //cout << "using LO...\n";
    double LO_rT, LO_rL;
    rho_LO(o,k,mu,LO_rT,LO_rL);
    double min_rT, min_rL;
    double _in[4] = {o,k,al_min,mu}; // use boundary of grid
    approx_rho(_in,rhoT_list,rhoL_list,min_rT,min_rL);

    double weight = alpha/al_min;
    rT = weight*min_rT+(1.-weight)*LO_rT;
    rL = weight*min_rL+(1.-weight)*LO_rL; return 1;
  }
  if (M>M_max+1e-5) { // use OPE expansion to extrapolate if M > M_max
    //cout << "using OPE...\n";
    double OPE_rT, OPE_rL;
    rho_OPE(o,k,alpha,mu,OPE_rT,OPE_rL);
    rT = OPE_rT;
    rL = OPE_rL; return 1;
  }
  else if (M<M_min-1e-5) { // use AMY approx. formula if M < M_min
    //cout << "using AMY...\n";
    double _in[4] = {sqrt(M_min*M_min+k*k),k,alpha,mu}; // use boundary of grid
    double min_rT, min_rL;
    approx_rho(_in,rhoT_list,rhoL_list,min_rT,min_rL);

    double weight = M/M_min;
    rT = weight*min_rT+(1.-weight)*rhoV_AMY(alpha,mu,k)/2.;
    rL = weight*min_rL; return 1;
  }

  rT = rhoT_list.interp(alpha,muB,M,k);
  rL = rhoL_list.interp(alpha,muB,M,k); return 0;
}

/*--------------------------------------------------------------------*/
// dimensionful emission rates

double B_(double x) { // kinematic factor
  if (4.*x>1.) 
    return 0.;
  else 
    return (1.+2.*x)*sqrt(1.-4.*x);
}

double ThermalPhoton::rate(struct ThermalPhoton::Table grid_T, struct ThermalPhoton::Table grid_L,
            double o, double k, double alpha_s, double muB, double T, double m_l) {
  // NB: quantities are defined in the local rest frame, i.e. u_\mu=(1,0,0,0)
  double M2 = o*o-k*k;

  double prefactor = B_(m_l*m_l/M2)*(2./3.)*pow(hbarC,-4.)*pow(alphaEM,2.); 
  // units: GeV^-4.fm^-4

  // m_l, T, muB, o and k must all be given in units of GeV!
  double in[4] = {o/T, k/T, alpha_s, muB/T}; // spf is a function of o/T etc.
  double rhoT_app, rhoL_app;
  approx_rho(in,grid_T,grid_L,rhoT_app,rhoL_app);

  return prefactor*nB(o/T)*pow(T,2.)*( 2.*rhoT_app + rhoL_app )/(3.*pow(M_PI,3.)*M2);
}


