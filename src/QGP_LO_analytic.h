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
    // void NetBaryonCorrection(double T, double muB, std::vector<double> &Eq,
    //                          std::vector<double> &eqrate_ptr);
};

#endif
