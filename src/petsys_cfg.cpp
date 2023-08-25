/*!
 * class Configurator implementation file
 *
 *
 * Author: Jian Zhou
 * Date: 11-21-2013
 *
 */

#include <petsys_cfg.h>

Configurator::~Configurator()
{
}

Configurator::Configurator()
{
}

void Configurator::processOptionFile(const char* filename, SYSTEMPARAMETERS& sysparms)
{
	SystemLog::write("processing configuration file %s ...\n", filename);

    std::ifstream input(filename);

    if (!input) {

		SystemLog::write("open option file failed\n");
		abort();

    } else {

        int lnum = 1;
        std::string option, op_temp, key, value;
        std::pair<std::string, int> vinfo;

        while (getline(input, option)) {

			SystemLog::write("checking line %d ... \n", lnum);
			
            Configurator::OP_state ops = isValidOption(option, key, value);

            switch (ops) {
            case OP_COMMENT_ONLY:
                // show comment
                SystemLog::write("comment: %s\n", option.c_str());
				break;
				
            case OP_EMPTY:
				SystemLog::write("nothing\n");
                break;

            case OP_BAD_FORMAT:
            	SystemLog::write("bad format, check '%s' at line %d\n", 
            		option.c_str(), lnum);
                abort();
                break;

            case OP_UNKNOWN:
                SystemLog::write("unknown option '%s'\n",
                            key.c_str());
                break;

            case OP_OK:
                // put key_value pairs in map
                SystemLog::write("found option '%s'\n", key.c_str());
                if (m_key_map.find(key) != m_key_map.end()) {
                    SystemLog::write("error, '%s' is duplicated\n",
                                key.c_str());
                    abort();
                }
                vinfo = std::make_pair(value, lnum);
                m_key_map[key] = vinfo;
                break;

            default:
                SystemLog::write("missing case at line %d\n", lnum);
                abort();
                break;
            }

            lnum ++;
        }
    }

    // check if option is missing or not
    SystemLog::write("validating minimum requirement ...\n");
    std::vector<std::string> sys_opts_mandatory;
    separateString(g_sysopts_mandatory, sys_opts_mandatory, Configurator::va_delim);
    bool missing = false;

    for (std::vector<std::string>::iterator it = sys_opts_mandatory.begin();
            it != sys_opts_mandatory.end(); it ++) {

        if (m_key_map.find((*it)) == m_key_map.end()) {
            SystemLog::write("missing option '%s'\n",
                        (*it).c_str());
            missing = true;
        }

    }

    if (missing) {
        SystemLog::write("error, configuration is incomplete!\n");
        abort();
    }

    // set default parameters and
    // then overwrite by the values read from configuration file
	SystemLog::write("extracting parameters ...\n");    
    setDefaultParms(sysparms);
    std::vector<std::string> sys_opts_full;
    separateString(g_sysopts, sys_opts_full, Configurator::va_delim);

    for (std::vector<std::string>::iterator it = sys_opts_full.begin();
            it != sys_opts_full.end(); it ++) {

        // find
        if (m_key_map.find((*it)) != m_key_map.end()) {
			
			SystemLog::write("checking '%s' at line %d ... ", 
				(*it).c_str(), m_key_map[(*it)].second);
            if (!fillSysParms((*it), m_key_map[(*it)].first, sysparms)) {
                abort();
            }
			SystemLog::write("ok.\n");
        } else {
            // some options can have default values, so just ignore the missing check
            //SystemLog::write("option '%s' not found, ignore ...\n", (*it).c_str());
        }
    }

}

void Configurator::setDefaultParms(SYSTEMPARAMETERS& sysparms)
{
    sysparms.crystal_gap_size.transaxial_front = 0.0;
    sysparms.crystal_gap_size.transaxial_back = 0.0;
    sysparms.crystal_gap_size.axial = 0.0;
    sysparms.tof_info.resolution = -1;
    sysparms.tof_info.bin_size = -1;
    sysparms.doi_info.resolution = -1;
    sysparms.doi_info.number_of_bins = -1;
    sysparms.detector_module_angular_offset = 0.0;
    sysparms.recon_output_folder = "";
    sysparms.recon_output_filename_prefix = "";
    sysparms.sensitivity = "";
    sysparms.initial_guess = "";
    sysparms.ipsf_transaxial = "";
    sysparms.ipsf_axial = "";
    sysparms.input_raw_data_format_type = -1;
    sysparms.warmup = false;
    sysparms.number_of_subsets_warmup = -1;
    sysparms.number_of_iterations_warmup = -1;
    sysparms.input_raw_data_file = "";
    sysparms.input_raw_data_format_type = -1;
    sysparms.iterative_algorithm_type = -1;
    sysparms.number_of_subsets = -1;
    sysparms.number_of_iterations = -1;
    sysparms.stepsize_for_intermediate_result = -1;
    sysparms.regularizer_model_type = -1;
    sysparms.number_of_radial_bins = -1;
    sysparms.regularizer_spatial_variant_weight = "";

    // sysparms.gb_alpha = 0;
    // sysparms.gb_CIP = "";
}

bool Configurator::isValidKey(const std::string& key,
                              const std::string& syskeys,
                              const char syskey_delim)
{
    std::vector<std::string> str_key_buff;
    separateString(syskeys, str_key_buff, syskey_delim);

    for (size_t i = 0; i < str_key_buff.size(); i ++) {
        if (str_key_buff[i].compare(key) == 0) {
            return true;
        }
    }

    return false;
}

template<typename T>
bool Configurator::isValidValue(const std::string& value,
                                const char va_delim,
                                const int num_of_values_required,
                                std::vector<T>& t_buff)
{
    std::vector<std::string> str_val_buff;
    separateString(value, str_val_buff, va_delim);

    if (str_val_buff.size() < size_t(num_of_values_required)) {
        SystemLog::write("too few parameters\n");
        return false;
    } else {
        if (str_val_buff.size() > size_t(num_of_values_required)) {			
			SystemLog::write("too many parameters\n");
        }
    }

    t_buff.clear();
    for (int i = 0; i < num_of_values_required; i ++) {
        t_buff.push_back(Configurator::string_to_T<T>(str_val_buff[i]));
    }

    return true;
}

Configurator::OP_state Configurator::isValidOption(const std::string& option,
        std::string& key,
        std::string& value)
{
    std::string op_temp = option;
    eraseCharFromString(op_temp, " \t\n\b");

    if (op_temp.size() == 0) { // double check
        return OP_EMPTY;
    }

    if (op_temp[0] == Configurator::cm_delim) {
        return OP_COMMENT_ONLY;
    }

    size_t pos_sharp = option.find_first_of(Configurator::cm_delim);
    if (pos_sharp == std::string::npos) {
    //    SystemLog::write("no option comment is found\n");
    }

    op_temp = option.substr(0, pos_sharp);
    size_t pos_equal = op_temp.find_first_of(Configurator::eq_delim);

    if (pos_equal == std::string::npos) {
        SystemLog::write("can't tell if it is an option or not\n");
        return OP_BAD_FORMAT;
    }

    std::string op_temp2 = op_temp;
    eraseCharFromString(op_temp2, " \t\n\b");

    if (op_temp2[0] == eq_delim ||
        op_temp2[op_temp2.length() - 1] == Configurator::eq_delim) {
        SystemLog::write("option is incomplete\n");
        return OP_BAD_FORMAT;
    }

    std::vector<std::string> s_buff;
    separateString(op_temp, s_buff, eq_delim);

//	for (int i = 0; i < s_buff.size(); i ++) {
//		std::cout << "#" << i+1 << ":" << "["
//				  << s_buff[i].length() << "]>"
//				  << s_buff[i] << "<" << std::endl;
//	}

    if (s_buff.size() > 2) {
    	SystemLog::write("more than one '=' is found\n");
        return OP_BAD_FORMAT;
    }

    key = s_buff[0];
    eraseCharFromString(key, " \t\b\n");
    value = s_buff[1];
    eraseCharFromString(value, " \t\b\n");

    // check if it is a valid key
    if (!isValidKey(key)) {
        return OP_UNKNOWN;
    }

    return OP_OK;
}

void Configurator::eraseCharFromString(std::string& s, const char c)
{
    size_t pos_c = s.find(c);

    while (pos_c != std::string::npos) {
        s.erase(pos_c, 1);
        pos_c = s.find(c);
    }
}

void Configurator::eraseCharFromString(std::string& s, const std::string& clist)
{
    for (size_t i = 0; i < clist.length(); i ++) {
        eraseCharFromString(s, clist[i]);
    }
}

void Configurator::separateString(const std::string& s,
                                  std::vector<std::string>& s_buff,
                                  const char delim)
{
    std::string s_temp = s, s_sub;

    size_t pos_delim = s_temp.find(delim);

    // special case
    if (pos_delim == std::string::npos) { 
    	// meaning only one parameter
    	s_buff.push_back(s);
        return;
    }

    while (pos_delim != std::string::npos) {
        s_sub = s_temp.substr(0, pos_delim);
        //
        s_buff.push_back(s_sub);
        s_temp.erase(0, pos_delim + 1);
        pos_delim = s_temp.find(delim);
    }

    // process last substring
    if (s_temp.length() > 0) {
        s_buff.push_back(s_temp);
    }
}

bool Configurator::fillSysParms(const std::string& key,
                                const std::string& value,
                                SYSTEMPARAMETERS& sysparms)
{
    std::vector<std::string> s_buff;

    // image size
    if (key.compare("image_size") == 0) {
        std::vector<int> t_buff;

        if (isValidValue<int>(value, Configurator::va_delim, 3, t_buff)) {
            sysparms.image_size.i = t_buff[0];
            sysparms.image_size.j = t_buff[1];
            sysparms.image_size.k = t_buff[2];
        } else {
            SystemLog::write("invalid parameters: '%s', "
            				 "note: require 3 positive integers.\n",
                        value.c_str());
            return false;
        }

        if (sysparms.image_size.i <= 0 ||
            sysparms.image_size.j <= 0 ||
            sysparms.image_size.k <= 0) {
        	SystemLog::write("invalid parameters: '%s', "
            				 "note: should be positive integers.\n",
                        	 value.c_str());
            return false;
        }
    }

    // voxel size
    if (key.compare("voxel_size") == 0) {
        std::vector<double> t_buff;

        if (isValidValue<double>(value, Configurator::va_delim, 3, t_buff)) {
            sysparms.voxel_size_i = t_buff[0];
            sysparms.voxel_size_j = t_buff[1];
            sysparms.voxel_size_k = t_buff[2];
        } else {
            SystemLog::write("invalid parameters '%s', "
            				 "note: require 3 positive numbers (unit in mm).\n",
                        	 value.c_str());
            return false;
        }

        if (sysparms.voxel_size_i <= 0 ||
            sysparms.voxel_size_j <= 0 ||
            sysparms.voxel_size_k <= 0) {
            SystemLog::write("invalide parameters '%s', "
            				 "note: should be positive real numbers.\n",
                        	 value.c_str());
            return false;
        }
    }

    // ring diameter
    if (key.compare("detector_ring_diameter") == 0) {
        sysparms.detector_ring_diameter = Configurator::string_to_T<double>(value);

        if (sysparms.detector_ring_diameter < 0) {
            SystemLog::write("invalid parameter '%s', "
            				 "note: positive real number only\n",
                       key.c_str(), value.c_str());
            return false;
        }
    }

    // crystal size
    if (key.compare("crystal_size") == 0) {
        std::vector<double> t_buff;

        if (isValidValue<double>(value, Configurator::va_delim, 4, t_buff)) {
            sysparms.crystal_size.transaxial_front = t_buff[0];
            sysparms.crystal_size.transaxial_back = t_buff[1];
            sysparms.crystal_size.axial = t_buff[2];
            sysparms.crystal_size.depth = t_buff[3];
        } else {
            SystemLog::write("invalid parameters '%s', "
            				 "note: require 4 positive numbers (unit in mm).\n",
                        	 value.c_str());
            return false;
        }

        if (sysparms.crystal_size.transaxial_front < 0 ||
            sysparms.crystal_size.transaxial_back < 0 ||
            sysparms.crystal_size.axial < 0 ||
            sysparms.crystal_size.depth < 0) {

            SystemLog::write("invalid parameters '%s', "
            				 "note: positive real numbers only\n",
                        	 value.c_str());

            return false;
        }
    }

    // crystal array size
    if (key.compare("crystal_array_size") == 0) {
        std::vector<int> t_buff;

        if (isValidValue<int>(value, Configurator::va_delim, 2, t_buff)) {
            sysparms.crystal_array_size.t = t_buff[0];
            sysparms.crystal_array_size.a = t_buff[1];
        } else {
            SystemLog::write("invalid parameters '%s', "
            				 "note: require 2 positive integers.\n",
                             value.c_str());
            return false;
        }

        if (sysparms.crystal_array_size.t <= 0 ||
            sysparms.crystal_array_size.a <= 0) {

            SystemLog::write("invalid parameters '%s', "
            				 "note: positive integers only.\n",
                        	 value.c_str());
            return false;
        }
    }

    // gap size
    if (key.compare("crystal_gap_size") == 0) {
        std::vector<double> t_buff;

        if (isValidValue<double>(value, Configurator::va_delim, 3, t_buff)) {
            sysparms.crystal_gap_size.transaxial_front = t_buff[0];
            sysparms.crystal_gap_size.transaxial_back = t_buff[1];
            sysparms.crystal_gap_size.axial = t_buff[2];
        } else {
            SystemLog::write("invalid parameters '%s', "
            				 "note: require 3 nonnegative real numbers (unit in mm).\n",
                        	 value.c_str());
            return false;
        }

        if (sysparms.crystal_gap_size.transaxial_front < 0 ||
                sysparms.crystal_gap_size.transaxial_back < 0 ||
                sysparms.crystal_gap_size.axial < 0) {

            SystemLog::write("invalid parameters '%s', "
            				 "note: nonnegative real numbers.\n",
                        	 value.c_str());
            return false;

        }
    }

    // detector module arrangement
    if (key.compare("number_of_detector_modules") == 0) {
        std::vector<int> t_buff;

        if (isValidValue<int>(value, Configurator::va_delim, 2, t_buff)) {
            sysparms.number_of_detector_modules.t = t_buff[0];
            sysparms.number_of_detector_modules.a = t_buff[1];
        } else {
            SystemLog::write("invalid parameters '%s', "
            				 "note: require 2 positive integers.\n",
                        	 value.c_str());
            return false;
        }

        if (sysparms.number_of_detector_modules.t <= 0 ||
                sysparms.number_of_detector_modules.a <= 0) {
            SystemLog::write("invalid parameters '%s', "
            				 "note: positive integers only.\n",
                        	 value.c_str());
            return false;
        }
    }

    // doi information
    if (key.compare("DOI_information") == 0) {
        std::vector<double> t_buff;

        if (isValidValue<double>(value, Configurator::va_delim, 2, t_buff)) {
            sysparms.doi_info.resolution = t_buff[0];
            sysparms.doi_info.number_of_bins = int(t_buff[1]);
        } else {
            SystemLog::write("invalid parameters '%s', "
            				 "note: 2 parameters, 1st: resolution (unit in mm) and "
            				 "2nd: number of DOI bins (positive integer).\n",
                         	 value.c_str());
            return false;
        }

        if (sysparms.doi_info.resolution <= 0.0 ||
            sysparms.doi_info.number_of_bins <= 0) {
            SystemLog::write("invalid parameters '%s', "
            				 "note: either resolution or number "
            				 "of bins has to be positive.\n",
                             value.c_str());
            return false;
        }
    }

    // angle offset
    if (key.compare("detector_module_angular_offset") == 0) {
        sysparms.detector_module_angular_offset = Configurator::string_to_T<double>(value);
    }

    // initial image
    if (key.compare("initial_guess") == 0) {
        sysparms.initial_guess = value; // only take one parameter
    }

    // recon output
    if (key.compare("reconstruction_output_setting") == 0) {
        std::vector<std::string> t_buff;

        if (isValidValue<std::string>(value, Configurator::va_delim, 2, t_buff)) {
            sysparms.recon_output_folder = t_buff[0];
            sysparms.recon_output_filename_prefix = t_buff[1];
        }
    }

    // sensitivity
    if (key.compare("sensitivity") == 0) {
        sysparms.sensitivity = value; // only take one parameter
    }

    // tof_info
    if (key.compare("TOF_information") == 0) {
        std::vector<double> t_buff;

        if (isValidValue<double>(value, Configurator::va_delim, 2, t_buff)) {
            sysparms.tof_info.resolution = t_buff[0];
            sysparms.tof_info.bin_size = t_buff[1];
        } else {
            SystemLog::write("invalid parameters '%s', "
            				 "note: 2 parameters, 1st: resolution (unit in ps), "
            				 "and 2nd: TOF bin size (unit in ps).\n",
                             value.c_str());
            return false;
        }

        if (sysparms.tof_info.resolution <= 0.0 ||
            sysparms.tof_info.bin_size <= 0) {
            SystemLog::write("invalid parameters '%s', "
            				 "note: either resolution and bin size has to be positive.\n",
                        	 value.c_str());
            return false;
        }
    }

    // image domain blurring matrix
    if (key.compare("iPSF_model") == 0) {
        std::vector<std::string> t_buff;

        if (isValidValue<std::string>(value, Configurator::va_delim, 2, t_buff)) {
            sysparms.ipsf_transaxial = t_buff[0];
            sysparms.ipsf_axial = t_buff[1];
        }
        
        if (sysparms.ipsf_transaxial.empty() || 
            sysparms.ipsf_transaxial.empty()) {
            SystemLog::write("warning: incomplete parameters, ignored.");
        }
        
    }

    // regularizer
    if (key.compare("regularizer_strength") == 0) {
        std::vector<std::string> t_buff;
        sysparms.regularizer_strength = Configurator::string_to_T<double>(value);

        if (sysparms.regularizer_strength < 0.0) {
            SystemLog::write("invalid parameter '%s', "
                        	 "note: should be a positive real number.\n",
                        	 value.c_str());
            return false;
        }
    }   
    
    // regularizer setting
    if (key.compare("regularizer_model_type") == 0) {
    	sysparms.regularizer_model_type = Configurator::string_to_T<int>(value);   	
    }
    if (key.compare("regularizer_potential_function_type") == 0) {
    	sysparms.regularizer_potential_function_type = Configurator::string_to_T<int>(value);  
    	if (sysparms.regularizer_potential_function_type < 0 ||
    		sysparms.regularizer_potential_function_type > 999) {
    		SystemLog::write("invalid parameter '%s', "
    						 "note: the type should be an integer.\n");
    		return false;
    	}
    }
    if (key.compare("regularizer_neighborhood_properties") == 0) {
    	std::vector<int> t_buff;
    	if (isValidValue<int>(value, Configurator::va_delim, 2, t_buff)) {
    		sysparms.regularizer_neighborhood_size = t_buff[0];
    		sysparms.regularizer_isotropic_or_not = t_buff[1];
    	} else {
    		SystemLog::write("invalid parameter '%s', "
    						 "note: should be in the form of [ "
    						 "size, isotropic or not (0 or 1) ].\n", value.c_str());
    		return false;
    	}
    	
    	if (sysparms.regularizer_neighborhood_size < 0) {
    		SystemLog::write("invalid parameter '%s', "
    						 "neighborhood size should be larger than 0\n", value.c_str());
    		return false;
    	}
    	
    	if ((sysparms.regularizer_isotropic_or_not != 0) &&
    		(sysparms.regularizer_isotropic_or_not != 1)) {
    		SystemLog::write("invalid parameter '%s', "
    						 "warning: value for penalty "
    						 "form (isotropic or not) must be either 0 or 1.\n", value.c_str());
    		SystemLog::write("an invalid value is detected [%d], "
    						 "set it back to defaul value: 0.\n", 
    						 sysparms.regularizer_isotropic_or_not);
    		sysparms.regularizer_isotropic_or_not = 0;
    	}
    }
    if (key.compare("regularizer_buildin_parameter_list") == 0) {
    	std::vector<double> t_buff;
    	sysparms.regularizer_number_of_buildin_parameters_detected = 0;
    	SystemLog::write("scanning parameter list '%s' ...\n", value.c_str());
    	int n;
    	for (n = 1; n <=5; n ++) {
    		SystemLog::write("trying %d ... \n", n);
    		if (isValidValue<double>(value, Configurator::va_delim, n, t_buff)) {
    			for (int k = 0; k < n; k ++) {
    				sysparms.regularizer_buildin_parameters[k] = t_buff[k];
    			}
    		} else {
    			break;
    		}
    	}
    	sysparms.regularizer_number_of_buildin_parameters_detected = (n-1);
    	SystemLog::write("number of buildin parameters detected: %d\n", 
    		sysparms.regularizer_number_of_buildin_parameters_detected);
    	bool valid = true;
    	for (int k = 0; k < sysparms.regularizer_number_of_buildin_parameters_detected; k ++) {
    		if (sysparms.regularizer_buildin_parameters[k] < 0) {
    			valid = false;
    			sysparms.regularizer_buildin_parameters[k] = 0.0;				 
    		}
    	}
    	if (!valid) {
    		SystemLog::write("invalid parameter list '%s', "
    						 "note: buildin parameter should be nonnegative.\n", value.c_str());
    		SystemLog::write("set back to zero already.\n");				 
    	}						 
    }

	// regularizer spatial variant weight
	if (key.compare("regularizer_spatial_variant_weight") == 0) {
        sysparms.regularizer_spatial_variant_weight = value;
	}

    // raw data format
    if (key.compare("input_raw_data_format_type") == 0) {
        std::vector<std::string> t_buff;
        sysparms.input_raw_data_format_type = Configurator::string_to_T<int>(value);
    }
	
	// raw data file
	if (key.compare("input_raw_data_file") == 0) {
		sysparms.input_raw_data_file = value;		
	}

	// warm-up setting
	if (key.compare("warmup_setting") == 0) {
		std::vector<int> t_buff;
		if (isValidValue<int>(value, Configurator::va_delim, 2, t_buff)) {
            sysparms.number_of_subsets_warmup = t_buff[0];
            sysparms.number_of_iterations_warmup = t_buff[1];
        }
        sysparms.warmup = true;
        if ((sysparms.number_of_subsets_warmup < 0) ||
           (sysparms.number_of_iterations_warmup < 0)) {
           	SystemLog::write("invalid parameter '%s', "
           					 "note: should be 2 positive integers\n", value.c_str());
           	SystemLog::write("warning: warmup setting will be disabled.");
           	sysparms.warmup = false;
        }
	}

	// iterative algorithm
	if (key.compare("iterative_algorithm_type") == 0) {
		sysparms.iterative_algorithm_type = Configurator::string_to_T<int>(value);
	}
	
	// iteration_setting
	if (key.compare("iteration_setting") == 0) {
		std::vector<int> t_buff;
		if (isValidValue<int>(value, Configurator::va_delim, 3, t_buff)) {
			sysparms.number_of_subsets = t_buff[0];
			sysparms.stepsize_for_intermediate_result = t_buff[1];
			sysparms.number_of_iterations = t_buff[2];
		} else {
			SystemLog::write("invalid parameter '%s', "
							 "note: should be: [number of subsets, "
							 "stepsize for intermediate result, "
							 "number of iterations].\n", value.c_str());
			return false;
		}
		
		if ((sysparms.number_of_subsets <= 0) ||
			(sysparms.number_of_iterations <= 0)) {
			SystemLog::write("invalid parameter '%s', "
           					 "note: number of subsets and number of iterations"
           					 "should be positive integers.\n", value.c_str());
           	return false;				 
		}
	}
	
	// number of radial bins
	if (key.compare("number_of_radial_bins") == 0) {
		sysparms.number_of_radial_bins = Configurator::string_to_T<int>(value);
	}


//     // Guobao's alpha
//     if (key.compare("gb_alpha") == 0) {
//         std::vector<std::string> t_buff;
//         sysparms.gb_alpha = Configurator::string_to_T<double>(value);

//         if (sysparms.gb_alpha > 1 || sysparms.gb_alpha < 0 ) {
//             SystemLog::write("invalid parameter '%s', "
//                              "note: gb_alpha should be between 0 and 1.\n",
//                              value.c_str());
//             return false;
//         }
//     } 
    
//     // gb_CIP
//     if (key.compare("gb_CIP") == 0) {
//         sysparms.gb_CIP = value; // only take one parameter
//     }
    
    
    return true;
}
