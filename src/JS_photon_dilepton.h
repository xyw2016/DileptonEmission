#ifndef JS_PHOTON_H_
#define JS_PHOTON_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <memory>

#include "./PhotonEmission.h"
#include "./Hydroinfo_h5_em.h"
#include "./Hydroinfo_MUSIC_em.h"
#include "./Stopwatch.h"
#include "./Arsenal.h"
#include "./ParameterReader.h"
#include "./gauss_quadrature.h"


namespace Photon_dilepton{

class JS_photon_dilepton{
    private:
        const std::string inputfile_;
    public:
	JS_photon_dilepton(std::string inputfile);
	~JS_photon_dilepton();

        std::shared_ptr<ParameterReader> paraRdr;
        Hydroinfo_MUSIC* hydroinfo_ptr;
        std::vector<float> photon_sp;




      std::vector<float> run(const std::vector<float>& bulkdata);

};

}

#endif
