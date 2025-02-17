#include "JS_photon_dilepton.h"
#include "./PhotonEmission.h"


namespace Photon_dilepton{
    
    JS_photon_dilepton::JS_photon_dilepton(std::string inputfile):
        inputfile_(inputfile){

            paraRdr = std::make_shared<ParameterReader>();
            paraRdr->readFromFile(inputfile_);
    }
    
    JS_photon_dilepton::~JS_photon_dilepton() {
    }


    void JS_photon_dilepton::run(const std::vector<float>& bulkdata,const int info_in_memory,const std::string ID)
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
	if (info_in_memory == 1){
        	hydroinfo_ptr->readHydroDatafromJS(bulkdata);
	}
	else{
        	int hydro_mode = 12;
        	int nskip_tau = 1;
		hydroinfo_ptr->set_evo_path("./evolution_all_xyeta_MUSIC.dat");
        	hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau);

	}

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
        std::ostringstream output_reset;
	output_reset << "ph_dlep_results/"<<ID<<"/";
        thermalPhotons.reset_output_path(output_reset.str());
        thermalPhotons.outputPhotonSpvn();
        
        delete hydroinfo_ptr;
        //return dilepton_sp;
        // clean up
        delete[] eta_ptr;
        delete[] etaweight_ptr;
        



    }

}
