/*!
 * class Configurator implementation file
 *
 *
 * Author: Jian Zhou
 * Date: 11-21-2013
 *
 */
#include <petsys_cscanner.h>

CScanner::CScanner() :
m_projector(0),
//m_ipsfmodel(0),
m_regularizer(0)
//m_sensitivity(0)
{
}

CScanner::~CScanner()
{
}

void CScanner::initialize(const char* opt_file)
{
    Configurator cfg;
    cfg.processOptionFile(opt_file, m_sysparms);
}

int CScanner::getNumberOfCrystalRings() const
{
    return m_sysparms.crystal_array_size.a * m_sysparms.number_of_detector_modules.a;
}

float CScanner::getCrystalRingPitchSize() const
{
    return m_sysparms.crystal_size.axial + m_sysparms.crystal_gap_size.axial;
}

float CScanner::getCrystalTransaxialPitchSize() const
{
    return m_sysparms.crystal_size.transaxial_front +
           m_sysparms.crystal_gap_size.transaxial_front;
}

int CScanner::getNumberOfRadialBins() const
{
	return m_sysparms.number_of_radial_bins;
}

Projector* CScanner::createProjector()
{
	// choose projector based on data format
	SystemLog::write("creating projector ...\n");
	
	// listmode based projector
	if (m_sysparms.input_raw_data_format_type == Projector::LM_PROJ) {
		
		SystemLog::write("creating listmode-based projector ...\n");
		
    	ListModeProjector* projector = new ListModeProjector;

    	projector->initialize(getNumberOfCrystalRings(),
    					 	  getCrystalRingPitchSize(),
    					 	  m_sysparms.number_of_detector_modules.t,
    					 	  m_sysparms.crystal_array_size.t,
    					 	  m_sysparms.detector_ring_diameter,
    					 	  getCrystalTransaxialPitchSize(),
    					 	  m_sysparms.crystal_size.depth,
    					 	  m_sysparms.image_size.i,
    					 	  m_sysparms.image_size.k,
    					 	  m_sysparms.voxel_size_i,
    					 	  m_sysparms.voxel_size_k,
    					 	  m_sysparms.tof_info.resolution,
    					 	  m_sysparms.tof_info.bin_size,
    					 	  m_sysparms.detector_module_angular_offset);
		
		// install iPSF model
		projector->initializeIPSFModel(m_sysparms.image_size,
									   m_sysparms.ipsf_transaxial.c_str(),
                	 				   m_sysparms.ipsf_axial.c_str());
		
		projector->initializeSensitivityImage(m_sysparms.sensitivity.c_str(), 
											  m_sysparms.image_size);

		// // install Kernel model
		// projector->initializeKernelModel(m_sysparms.image_size,
		// 							   m_sysparms.kernel_matrix.c_str());

		
		return projector;
		
	} else {
	
		// sinogram based projector
		if (m_sysparms.input_raw_data_format_type == Projector::SINO_PROJ) {
			//
			// create a sinogram-based projector
			// 	
			SystemLog::write("sinogram-based projector not supported yet.\n");
			return 0;
			
		} else {
		
			SystemLog::write("unknown projector type %d, not supported yet.\n",
							 m_sysparms.input_raw_data_format_type);
			return 0;
			
		}	
	
	}
}

/*
IPSFModel* CScanner::createIPSFModel()
{
	SystemLog::write("creating iPSF model ...\n");
	
    IPSFModel* ipsf = new IPSFModel(m_sysparms.image_size);
    if (!ipsf->initialize(m_sysparms.ipsf_transaxial.c_str(),
                	 m_sysparms.ipsf_axial.c_str())) {
    	SystemLog::write("error ... %s, %s\n", 
    		m_sysparms.ipsf_transaxial.c_str(), m_sysparms.ipsf_axial.c_str());
    }
	
    return ipsf;
}
*/

Regularizer* CScanner::createRegularizer()
{
	SystemLog::write("creating regularizer ...\n");

	if (m_sysparms.regularizer_model_type == Regularizer::PAIRWISE_MRF) {
		
		PairwiseMRFRegularizer* regularizer = new PairwiseMRFRegularizer();
		regularizer->initialize(m_sysparms.regularizer_potential_function_type,
								m_sysparms.regularizer_neighborhood_size,
								m_sysparms.regularizer_isotropic_or_not,
								m_sysparms.regularizer_buildin_parameters,
								m_sysparms.regularizer_number_of_buildin_parameters_detected);
		
		if (!m_sysparms.regularizer_spatial_variant_weight.empty()) {

			Image<float> spw(m_sysparms.image_size);
			if (!spw.read(m_sysparms.regularizer_spatial_variant_weight.c_str())) {
				SystemLog::write("reading spatial variant weight [in file `%s'] failed ... ignore\n",
					m_sysparms.regularizer_spatial_variant_weight.c_str());
			} else {
				regularizer->setSpatialVariantWeight(spw);
				SystemLog::write("spatial variant weight in file `%s' ...\n",
					m_sysparms.regularizer_spatial_variant_weight.c_str());
			}
		} else {
			SystemLog::write("no spatial variant weight specified ... ignore\n");
		}


		return regularizer;
	} else {
		
		if (m_sysparms.regularizer_model_type == Regularizer::PATCH) {
			
			PatchBasedRegularizer* regularizer = new PatchBasedRegularizer();
			regularizer->initialize(m_sysparms.regularizer_potential_function_type,
								m_sysparms.regularizer_neighborhood_size,
								m_sysparms.regularizer_isotropic_or_not,
								m_sysparms.regularizer_buildin_parameters,
								m_sysparms.regularizer_number_of_buildin_parameters_detected);
			
			if (!m_sysparms.regularizer_spatial_variant_weight.empty()) {

				Image<float> spw(m_sysparms.image_size);
				if (!spw.read(m_sysparms.regularizer_spatial_variant_weight.c_str())) {
					SystemLog::write("reading spatial variant weight [in file `%s'] failed ... ignore\n",
						m_sysparms.regularizer_spatial_variant_weight.c_str());
				} else {
					regularizer->setSpatialVariantWeight(spw);
					SystemLog::write("spatial variant weight in file `%s' ...\n",
						m_sysparms.regularizer_spatial_variant_weight.c_str());
				}
			} else {
				SystemLog::write("no spatial variant weight specified ... ignore\n");
			}

			return regularizer;	
				
		} else {
		
			SystemLog::write("unknown regularization model %d, not supported yet.\n",
				m_sysparms.regularizer_model_type);
			return 0;
			
		}
		
	}
	
}

const SYSTEMPARAMETERS& CScanner::getSystemParameters() const
{
	return (m_sysparms);
}

void CScanner::runRecon()
{
	// create projector, psfmodel, and regularizer
	m_projector = createProjector();
	if (!m_projector) {
		abort();
	}
	
	m_regularizer = createRegularizer();
	if (!m_regularizer) {
		abort();
	}

	// handle raw data
	m_projector->readRawData(m_sysparms.input_raw_data_file.c_str());
			
	// read in initial guess
	SystemLog::write("reading initial guess, file in '%s' ...\n", 
		m_sysparms.initial_guess.c_str());
	Image<float> reconstruction(m_sysparms.image_size);
	if (!reconstruction.read(m_sysparms.initial_guess.c_str())) {
		SystemLog::write("warning: read initial guess failed [%s], "
		"use default initial guess.\n", m_sysparms.initial_guess.c_str());
		reconstruction.set(1.0);
	}

	// warm up, ML reconstruction
	if (m_sysparms.warmup) { // for better convergence, PL requires a good initialization
	
		SystemLog::write("warmup ...\n");
		m_projector->makeSubsets(m_sysparms.number_of_subsets_warmup);
		EMBasedPLReconstructor::run(m_projector, 
									m_regularizer,
									m_sysparms.number_of_iterations_warmup,
									-1, 0.0, 
									m_sysparms.recon_output_folder.c_str(),
									m_sysparms.recon_output_filename_prefix.c_str(), 
									reconstruction);
		
		// reconstruction.write("warmup"); // test only!
        char str[512];
		sprintf(str, "%s/%s.warmup", m_sysparms.recon_output_folder.c_str(), m_sysparms.recon_output_filename_prefix.c_str());
		reconstruction.write(str);


	} else {
		SystemLog::write("no recon warm-up ...\n");
	}
	
	SystemLog::write("run PL recon ...\n");


		printf("Hello0.0\n");


	if (m_sysparms.warmup) {
		if (m_sysparms.number_of_subsets_warmup != m_sysparms.number_of_subsets)
			m_projector->makeSubsets(m_sysparms.number_of_subsets);
	} else {
		m_projector->makeSubsets(m_sysparms.number_of_subsets);
	}

	switch (m_sysparms.iterative_algorithm_type) {

		case EM:
			SystemLog::write("EM selected ...\n");
			EMBasedPLReconstructor::run(m_projector, 
										m_regularizer,
										m_sysparms.number_of_iterations,
										m_sysparms.stepsize_for_intermediate_result, 
										m_sysparms.regularizer_strength, 
										m_sysparms.recon_output_folder.c_str(),
										m_sysparms.recon_output_filename_prefix.c_str(), 
										reconstruction);
			break;
			
		case ADMM:
			SystemLog::write("ADMM selected ... ");
			ADMMBasedPLReconstructor::mu = m_sysparms.regularizer_strength * 2.0;
			SystemLog::write("beta = %.2e, mu = %.2e...\n", 
							 m_sysparms.regularizer_strength, 
							 ADMMBasedPLReconstructor::mu);
			ADMMBasedPLReconstructor::run(m_projector, 
										  m_regularizer,
										  m_sysparms.number_of_iterations,
										  m_sysparms.stepsize_for_intermediate_result, 
										  m_sysparms.regularizer_strength, 
										  m_sysparms.recon_output_folder.c_str(),
										  m_sysparms.recon_output_filename_prefix.c_str(), 
										  reconstruction);
			break;
		
		default:
			SystemLog::write("unknown iterative algorithm type [%d]\n,"
							 "valid type: 0 (EM) ...");
			abort();
			break;
	}
	
	printf("Hello1\n");

	// dump to disk
	char str[512];				 
	sprintf(str, "%s/%s.img", 
        	m_sysparms.recon_output_folder.c_str(),
        	m_sysparms.recon_output_filename_prefix.c_str());
    SystemLog::write("save to file '%s' ...\n", str);
    reconstruction.write(str);				 		
}




void CScanner::runRecon2()
{
	// create projector, psfmodel, and regularizer
	m_projector = createProjector();
	if (!m_projector) {
		abort();
	}
	
	m_regularizer = createRegularizer();
	if (!m_regularizer) {
		abort();
	}

	// handle raw data
	m_projector->readRawData(m_sysparms.input_raw_data_file.c_str());
			
	// read in initial guess
	SystemLog::write("reading initial guess, file in '%s' ...\n", 
		m_sysparms.initial_guess.c_str());
	Image<float> reconstruction(m_sysparms.image_size);
	if (!reconstruction.read(m_sysparms.initial_guess.c_str())) {
		SystemLog::write("warning: read initial guess failed [%s], "
		"use default initial guess.\n", m_sysparms.initial_guess.c_str());
		reconstruction.set(1.0);
	}

	// warm up, ML reconstruction
	if (m_sysparms.warmup) { // for better convergence, PL requires a good initialization
	
		SystemLog::write("warmup ...\n");
		m_projector->makeSubsets(m_sysparms.number_of_subsets_warmup);
		EMBasedPLReconstructor::run(m_projector, 
									m_regularizer,
									m_sysparms.number_of_iterations_warmup,
									-1, 0.0, 
									m_sysparms.recon_output_folder.c_str(),
									m_sysparms.recon_output_filename_prefix.c_str(), 
									reconstruction);
		
		reconstruction.write("warmup"); // test only!

	} else {
		SystemLog::write("no recon warm-up ...\n");
	}
	


	SystemLog::write("run PL recon ...\n");


		printf("Hello2.0\n");


	if (m_sysparms.warmup) {
		if (m_sysparms.number_of_subsets_warmup != m_sysparms.number_of_subsets)
			m_projector->makeSubsets(m_sysparms.number_of_subsets);
	} else {
		m_projector->makeSubsets(m_sysparms.number_of_subsets);
	}

	switch (m_sysparms.iterative_algorithm_type) {

		case EM:
			SystemLog::write("EM selected ...\n");
			EMBasedPLReconstructor::run(m_projector, 
										m_regularizer,
										m_sysparms.number_of_iterations,
										m_sysparms.stepsize_for_intermediate_result, 
										m_sysparms.regularizer_strength, 
										m_sysparms.recon_output_folder.c_str(),
										m_sysparms.recon_output_filename_prefix.c_str(), 
										reconstruction);
			break;
			
		case ADMM:
			SystemLog::write("ADMM selected ... ");
			ADMMBasedPLReconstructor::mu = m_sysparms.regularizer_strength * 2.0;
			SystemLog::write("beta = %.2e, mu = %.2e...\n", 
							 m_sysparms.regularizer_strength, 
							 ADMMBasedPLReconstructor::mu);
			ADMMBasedPLReconstructor::run(m_projector, 
										  m_regularizer,
										  m_sysparms.number_of_iterations,
										  m_sysparms.stepsize_for_intermediate_result, 
										  m_sysparms.regularizer_strength, 
										  m_sysparms.recon_output_folder.c_str(),
										  m_sysparms.recon_output_filename_prefix.c_str(), 
										  reconstruction);
			break;


		case KEM:
			SystemLog::write("KEM selected ...\n");
			KEMBasedPLReconstructor::run(m_projector, 
										m_regularizer,
										m_sysparms.number_of_iterations,
										m_sysparms.stepsize_for_intermediate_result, 
										m_sysparms.regularizer_strength, 
										m_sysparms.recon_output_folder.c_str(),
										m_sysparms.recon_output_filename_prefix.c_str(), 
										reconstruction);
			break;
			

		default:
			SystemLog::write("unknown iterative algorithm type [%d]\n,"
							 "valid type: 0 (EM) ...");
			abort();
			break;
	}
	
	printf("Hello2\n");
	/* save time and storage
	// dump to disk
	char str[512];				 
	sprintf(str, "%s/%s.img", 
        	m_sysparms.recon_output_folder.c_str(),
        	m_sysparms.recon_output_filename_prefix.c_str());
    SystemLog::write("save to file '%s' ...\n", str);
    reconstruction.write(str);	
	*/

}






// void CScanner::runRecon3()
// {
// 	// create projector, psfmodel, and regularizer
// 	m_projector = createProjector();
// 	if (!m_projector) {
// 		abort();
// 	}
	
// 	m_regularizer = createRegularizer();
// 	if (!m_regularizer) {
// 		abort();
// 	}

// 	// handle raw data
// 	m_projector->readRawData(m_sysparms.input_raw_data_file.c_str());
			
// 	// read in initial guess
// 	SystemLog::write("reading initial guess, file in '%s' ...\n", 
// 		m_sysparms.initial_guess.c_str());
// 	Image<float> reconstruction(m_sysparms.image_size);
// 	if (!reconstruction.read(m_sysparms.initial_guess.c_str())) {
// 		SystemLog::write("warning: read initial guess failed [%s], "
// 		"use default initial guess.\n", m_sysparms.initial_guess.c_str());
// 		reconstruction.set(1.0);
// 	}


//     Image<float> gb_CIPimg(m_sysparms.image_size);
    
//     if (m_sysparms.gb_alpha > 0)
//     {
//         SystemLog::write("reading gb_CIP, file in '%s' ...\n", m_sysparms.gb_CIP.c_str());
//         if (!gb_CIPimg.read(m_sysparms.gb_CIP.c_str()))
//         {
//             SystemLog::write("error: read gb_CIP failed in '%s'\n", m_sysparms.gb_CIP.c_str());
//             abort();
//         }
//     }
    
    
// 	// warm up, ML reconstruction
// 	if (m_sysparms.warmup) { // for better convergence, PL requires a good initialization
	
// 		SystemLog::write("warmup ...\n");
// 		m_projector->makeSubsets(m_sysparms.number_of_subsets_warmup);
// 		CIPEMBasedPLReconstructor::run(m_projector, 
// 									m_regularizer,
// 									m_sysparms.number_of_iterations_warmup,
// 									-1, 0.0, 
// 									m_sysparms.recon_output_folder.c_str(),
// 									m_sysparms.recon_output_filename_prefix.c_str(), 
// 									reconstruction,
//                                     0,  // warm up does not require gb_alpha. Set to 0. 
//                                     gb_CIPimg);
		
// 		reconstruction.write("warmup"); // test only!

// 	} else {
// 		SystemLog::write("no recon warm-up ...\n");
// 	}
	
// 	SystemLog::write("run PL recon ...\n");
// 	if (m_sysparms.warmup) {
// 		if (m_sysparms.number_of_subsets_warmup != m_sysparms.number_of_subsets)
// 			m_projector->makeSubsets(m_sysparms.number_of_subsets);
// 	} else {
// 		m_projector->makeSubsets(m_sysparms.number_of_subsets);
// 	}

// 	switch (m_sysparms.iterative_algorithm_type) {

// 		case EM:
// 			SystemLog::write("EM selected ...\n");
// 			EMBasedPLReconstructor::run(m_projector, 
// 										m_regularizer,
// 										m_sysparms.number_of_iterations,
// 										m_sysparms.stepsize_for_intermediate_result, 
// 										m_sysparms.regularizer_strength, 
// 										m_sysparms.recon_output_folder.c_str(),
// 										m_sysparms.recon_output_filename_prefix.c_str(), 
// 										reconstruction);
// 			break;
			
// 		case ADMM:
// 			SystemLog::write("ADMM selected ... ");
// 			ADMMBasedPLReconstructor::mu = m_sysparms.regularizer_strength * 2.0;
// 			SystemLog::write("beta = %.2e, mu = %.2e...\n", 
// 							 m_sysparms.regularizer_strength, 
// 							 ADMMBasedPLReconstructor::mu);
// 			ADMMBasedPLReconstructor::run(m_projector, 
// 										  m_regularizer,
// 										  m_sysparms.number_of_iterations,
// 										  m_sysparms.stepsize_for_intermediate_result, 
// 										  m_sysparms.regularizer_strength, 
// 										  m_sysparms.recon_output_folder.c_str(),
// 										  m_sysparms.recon_output_filename_prefix.c_str(), 
// 										  reconstruction);
// 			break;


// 		case KEM:
// 			SystemLog::write("KEM selected ...\n");
// 			KEMBasedPLReconstructor::run(m_projector, 
// 										m_regularizer,
// 										m_sysparms.number_of_iterations,
// 										m_sysparms.stepsize_for_intermediate_result, 
// 										m_sysparms.regularizer_strength, 
// 										m_sysparms.recon_output_folder.c_str(),
// 										m_sysparms.recon_output_filename_prefix.c_str(), 
// 										reconstruction);
// 			break;
			

// 		case CIPEM:
// 			SystemLog::write("KEM selected ...\n");
// 			CIPEMBasedPLReconstructor::run(m_projector, 
// 										m_regularizer,
// 										m_sysparms.number_of_iterations,
// 										m_sysparms.stepsize_for_intermediate_result, 
// 										m_sysparms.regularizer_strength, 
// 										m_sysparms.recon_output_folder.c_str(),
// 										m_sysparms.recon_output_filename_prefix.c_str(), 
// 										reconstruction,
//                                         m_sysparms.gb_alpha,
//                                         gb_CIPimg);
// 			break;
			


// 		default:
// 			SystemLog::write("unknown iterative algorithm type [%d]\n,"
// 							 "valid type: 0 (EM) ...");
// 			abort();
// 			break;
// 	}
	
// 	// dump to disk
// 	char str[512];				 
// 	sprintf(str, "%s/%s.img", 
//         	m_sysparms.recon_output_folder.c_str(),
//         	m_sysparms.recon_output_filename_prefix.c_str());
//     SystemLog::write("save to file '%s' ...\n", str);
//     reconstruction.write(str);				 		
// }









void CScanner::runfplinesection()
{
	// create projector, psfmodel, and regularizer
	m_projector = createProjector();
	if (!m_projector) {
		abort();
	}
	
	m_regularizer = createRegularizer();
	if (!m_regularizer) {
		abort();
	}

	// handle raw data
	m_projector->readRawData(m_sysparms.input_raw_data_file.c_str());
			
	// read in initial guess
	SystemLog::write("reading initial guess, file in '%s' ...\n", 
		m_sysparms.initial_guess.c_str());
	Image<float> reconstruction(m_sysparms.image_size);
	if (!reconstruction.read(m_sysparms.initial_guess.c_str())) {
		SystemLog::write("warning: read initial guess failed [%s], "
		"use default initial guess.\n", m_sysparms.initial_guess.c_str());
		reconstruction.set(1.0);
	}

	// warm up, ML reconstruction
	if (m_sysparms.warmup) { // for better convergence, PL requires a good initialization
	
		SystemLog::write("warmup ...\n");
		m_projector->makeSubsets(m_sysparms.number_of_subsets_warmup);
		EMBasedPLReconstructor::run(m_projector, 
									m_regularizer,
									m_sysparms.number_of_iterations_warmup,
									-1, 0.0, 
									m_sysparms.recon_output_folder.c_str(),
									m_sysparms.recon_output_filename_prefix.c_str(), 
									reconstruction);
		
		reconstruction.write("warmup"); // test only!

	} else {
		SystemLog::write("no recon warm-up ...\n");
	}
	
	SystemLog::write("run PL recon ...\n");



		printf("Hello3.0\n");



	if (m_sysparms.warmup) {
		if (m_sysparms.number_of_subsets_warmup != m_sysparms.number_of_subsets)
			m_projector->makeSubsets(m_sysparms.number_of_subsets);
	} else {
		m_projector->makeSubsets(m_sysparms.number_of_subsets);
	}

	switch (m_sysparms.iterative_algorithm_type) {

		case EM:
			SystemLog::write("EM selected ...\n");
			EMBasedPLReconstructor::run(m_projector, 
										m_regularizer,
										m_sysparms.number_of_iterations,
										m_sysparms.stepsize_for_intermediate_result, 
										m_sysparms.regularizer_strength, 
										m_sysparms.recon_output_folder.c_str(),
										m_sysparms.recon_output_filename_prefix.c_str(), 
										reconstruction);
			break;
			
		case ADMM:
			SystemLog::write("ADMM selected ... ");
			ADMMBasedPLReconstructor::mu = m_sysparms.regularizer_strength * 2.0;
			SystemLog::write("beta = %.2e, mu = %.2e...\n", 
							 m_sysparms.regularizer_strength, 
							 ADMMBasedPLReconstructor::mu);
			ADMMBasedPLReconstructor::run(m_projector, 
										  m_regularizer,
										  m_sysparms.number_of_iterations,
										  m_sysparms.stepsize_for_intermediate_result, 
										  m_sysparms.regularizer_strength, 
										  m_sysparms.recon_output_folder.c_str(),
										  m_sysparms.recon_output_filename_prefix.c_str(), 
										  reconstruction);
			break;
		
		default:
			SystemLog::write("unknown iterative algorithm type [%d]\n,"
							 "valid type: 0 (EM) ...");
			abort();
			break;
	}
	
	// dump to disk
	char str[512];				 
	sprintf(str, "%s/%s.img", 
        	m_sysparms.recon_output_folder.c_str(),
        	m_sysparms.recon_output_filename_prefix.c_str());
    SystemLog::write("save to file '%s' ...\n", str);
    reconstruction.write(str);				 		
}





void CScanner::runexplmfp()
{
	// create projector, psfmodel, and regularizer
	m_projector = createProjector();
	if (!m_projector) {
		abort();
	}
	
	m_regularizer = createRegularizer();
	if (!m_regularizer) {
		abort();
	}

	// handle raw data
	m_projector->readRawData(m_sysparms.input_raw_data_file.c_str());
			
	// read in initial guess
	SystemLog::write("reading initial guess, file in '%s' ...\n", 
		m_sysparms.initial_guess.c_str());
	Image<float> attnimg(m_sysparms.image_size);
	if (!attnimg.read(m_sysparms.initial_guess.c_str())) {
		SystemLog::write("warning: read initial guess failed [%s], "
		"use default initial guess.\n", m_sysparms.initial_guess.c_str());
		attnimg.set(1.0);
	}

	// // warm up, ML attnimg
	// if (m_sysparms.warmup) { // for better convergence, PL requires a good initialization
	
	SystemLog::write("runexplmfp ...\n");
	m_projector->makeSubsets(m_sysparms.number_of_subsets_warmup);
	ExpLineForwardProjector::run(m_projector, 
								m_regularizer,
								m_sysparms.number_of_iterations_warmup,
								-1, 0.0, 
								m_sysparms.recon_output_folder.c_str(),
								m_sysparms.recon_output_filename_prefix.c_str(), 
								attnimg);
		
	// 	explmfp.write("warmup"); // test only!

	// // } else {
	// // 	SystemLog::write("no recon warm-up ...\n");
	// // }
	
	// // dump to disk
	// char str[512];				 
	// sprintf(str, "%s/%s.sino", 
 //        	m_sysparms.recon_output_folder.c_str(),
 //        	m_sysparms.recon_output_filename_prefix.c_str());
 //    SystemLog::write("save to file '%s' ...\n", str);
 //    explmfp.write(str);				 		
}



void CScanner::runlmfp()
{
	// create projector, psfmodel, and regularizer
	m_projector = createProjector();
	if (!m_projector) {
		abort();
	}
	
	m_regularizer = createRegularizer();
	if (!m_regularizer) {
		abort();
	}

	// handle raw data
	m_projector->readRawData(m_sysparms.input_raw_data_file.c_str());
			
	// read in initial guess
	SystemLog::write("reading initial guess, file in '%s' ...\n", 
		m_sysparms.initial_guess.c_str());
	Image<float> attnimg(m_sysparms.image_size);
	if (!attnimg.read(m_sysparms.initial_guess.c_str())) {
		SystemLog::write("warning: read initial guess failed [%s], "
		"use default initial guess.\n", m_sysparms.initial_guess.c_str());
		attnimg.set(1.0);
	}

	// // warm up, ML attnimg
	// if (m_sysparms.warmup) { // for better convergence, PL requires a good initialization
	
	SystemLog::write("runlmfp ...\n");
	m_projector->makeSubsets(m_sysparms.number_of_subsets_warmup);
	LineForwardProjector::run(m_projector, 
								m_regularizer,
								m_sysparms.number_of_iterations_warmup,
								-1, 0.0, 
								m_sysparms.recon_output_folder.c_str(),
								m_sysparms.recon_output_filename_prefix.c_str(), 
								attnimg);
		
	// 	explmfp.write("warmup"); // test only!

	// // } else {
	// // 	SystemLog::write("no recon warm-up ...\n");
	// // }
	
	// // dump to disk
	// char str[512];				 
	// sprintf(str, "%s/%s.sino", 
 //        	m_sysparms.recon_output_folder.c_str(),
 //        	m_sysparms.recon_output_filename_prefix.c_str());
 //    SystemLog::write("save to file '%s' ...\n", str);
 //    explmfp.write(str);				 		
}



