#ifndef JS_DILEPTON_H_
#define JS_DILEPTON_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <memory>
//#include <omp.h>
#ifndef _OPENMP
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif

#include "./Hydroinfo_h5_dilepton.h"
#include "./Hydroinfo.h"
#include "./Stopwatch.h"
#include "./Arsenal.h"
#include "./ParameterReader.h"
#include "./gauss_quadrature.h"


namespace Dilepton{

class JS_dilepton{
   private:
      const std::string inputfile_;
      
   
   public:
      JS_dilepton(std::string inputfile);
      ~JS_dilepton();

      std::shared_ptr<ParameterReader> paraRdr;
      Hydroinfo* hydroinfo_ptr;
      std::vector<float> dilepton_sp;




      std::vector<float> run(const std::vector<float>& bulkdata);

};


}
#endif   // SRC_PHOTONEMISSION_H_

