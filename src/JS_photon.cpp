#include "JS_photon.h"
#include "./PhotonEmission.h"


namespace Photon{
    
    JS_photon::JS_photon(std::string inputfile):
        inputfile_(inputfile){

            paraRdr = std::make_shared<ParameterReader>();
            paraRdr->readFromFile(inputfile_);
    }
    
    JS_photon::~JS_photon() {
    }


    std::vector<float> JS_photon::run(const std::vector<float>& bulkdata)
    {
        int neta = paraRdr->getVal("neta");
        double eta_i = paraRdr->getVal("eta_i");
        double eta_f = paraRdr->getVal("eta_f");
        double* eta_ptr = new double[neta];
        double* etaweight_ptr = new double[neta];
        gauss_quadrature(neta, 1, 0.0, 0.0, eta_i, eta_f, eta_ptr, etaweight_ptr);
        PhotonEmission thermalPhotons(paraRdr);
        int hydro_flag = paraRdr->getVal("hydro_flag");
        if(hydro_flag!=2){
	   std::cout<<" Error: wrong hydro_flag in JS_photon.cpp"<< std::endl;
	}

	hydroinfo_ptr = new Hydroinfo_MUSIC();
        hydroinfo_ptr->readHydroDatafromJS(bulkdata);
        // calculate thermal photons from the hydro medium
        if (hydroinfo_ptr->isBoostInvariant()) {
            thermalPhotons.calPhotonemission(hydroinfo_ptr, eta_ptr,
                                             etaweight_ptr);
        } else {
            thermalPhotons.calPhotonemission_3d(hydroinfo_ptr);
        }
        
         
        // sum up all channels and compute thermal photon spectra and vn
        thermalPhotons.calPhoton_SpvnpT_individualchannel();
        thermalPhotons.calPhoton_total_SpMatrix();
        thermalPhotons.calPhoton_total_Spvn();
        int eventid = 0; 
        thermalPhotons.SavePhotoninJS(eventid);
        
        delete hydroinfo_ptr;
        //return dilepton_sp;
        return photon_sp;
        



    }

}
