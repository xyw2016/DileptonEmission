#include "JS_dilepton.h"
#include "./PhotonEmission.h"


namespace Dilepton{

    JS_dilepton::JS_dilepton(std::string inputfile):
        inputfile_(inputfile){

            paraRdr = std::make_shared<ParameterReader>();
            paraRdr->readFromFile(inputfile_);
    }
    
    JS_dilepton::~JS_dilepton() {
    }

    std::vector<float> JS_dilepton::run(const std::vector<float>& bulkdata)
    {
        //PhotonEmission thermalPhotons(paraRdr);
        PhotonEmission thermalPhotons(paraRdr);

        hydroinfo_ptr = new Hydroinfo();

        int nskip_tau = 1;
        // to_do, input_hydro_info_from_js
        hydroinfo_ptr->readHydroDatafromJS(bulkdata);
        thermalPhotons.calPhotonemission_3d(hydroinfo_ptr);
        dilepton_sp = thermalPhotons.PassdileptonspectratoJS();
         
        
        delete hydroinfo_ptr;
        return dilepton_sp;
        



    }


}