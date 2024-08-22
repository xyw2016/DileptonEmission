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
#include "./Hydroinfo_h5_photon.h"
#include "./Hydroinfo_MUSIC_photon.h"
#include "./Stopwatch.h"
#include "./Arsenal.h"
#include "./ParameterReader.h"
#include "./gauss_quadrature.h"


namespace Photon{

class JS_photon{
    private:
        const std::string inputfile_;
    public:
	JS_photon(std::string inputfile);
	~JS_photon();

        std::shared_ptr<ParameterReader> paraRdr;
        Hydroinfo_MUSIC* hydroinfo_ptr;
        std::vector<float> photon_sp;




      std::vector<float> run(const std::vector<float>& bulkdata);

};

}

#endif
