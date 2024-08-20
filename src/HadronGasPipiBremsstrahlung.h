#ifndef SRC_HADRONGASPIPIBREMSSTRAHLUNG
#define SRC_HADRONGASPIPIBREMSSTRAHLUNG

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

namespace Photon{
class HadronGasPipiBremsstrahlung : public ThermalPhoton {
 public:
    HadronGasPipiBremsstrahlung(std::shared_ptr<ParameterReader> paraRdr_in,
                                std::string emissionProcess);
    ~HadronGasPipiBremsstrahlung() {}
    void analyticRates(double T, std::vector<double> &Eq,
                        std::vector<double> &eqrate_ptr);
};
}
#endif
