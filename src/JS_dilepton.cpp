#include "JS_dilepton.h"
#include "./PhotonEmission.h"


namespace Dilepton{

    JS_dilepton::JS_dilepton(std::string inputfile,std::string table_path):
        inputfile_(inputfile), table_path_(table_path){

            paraRdr = std::make_shared<ParameterReader>();
            paraRdr->readFromFile(inputfile_);
    }
    
    JS_dilepton::~JS_dilepton() {
    }

    void JS_dilepton::run()
    {
        //PhotonEmission thermalPhotons(paraRdr);
        PhotonEmission thermalPhotons(paraRdr);

        hydroinfo_ptr = new Hydroinfo();

        int nskip_tau = 1;
        // to_do, input_hydro_info_from_js
        //hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau);
        thermalPhotons.calPhotonemission_3d(hydroinfo_ptr);
        
        delete hydroinfo_ptr;



    }


}