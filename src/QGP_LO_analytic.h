#ifndef SRC_QGP_LO_A
#define SRC_QGP_LO_A

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

class QGP_LO_analytic : public ThermalPhoton {
 public:
    QGP_LO_analytic(std::shared_ptr<ParameterReader> paraRdr_in,
                        std::string emissionProcess);
    ~QGP_LO_analytic() {}
    void FiniteBaryonRates(double T, double muB, double rhoB_over_eplusp, double Eq, 
        double M_ll, double &eqrate_ptr, double &eqrateT_ptr, double &eqrateL_ptr, double &diffrate_ptr, int include_diff_deltaf);
};

#endif
