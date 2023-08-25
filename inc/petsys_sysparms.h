/*!
 * Global structure of system parameter
 *
 % Modified from old structure
 %
 * Author: Jian Zhou
 * Date: 11-21-2013
 *
 */
#ifndef PETSYS_SYSPARMS_H
#define PETSYS_SYSPARMS_H

namespace UCDPETSYS
{
///
/// \brief Parameter keywords (strings).
///
/// They are options user can supply from outside.
/// Keywords have to be separated by comma.
///
/// Changed, use full name, keep options only related to listmode recon
///
static std::string g_sysopts("image_size,"
                             "voxel_size,"
                             "detector_ring_diameter,"
                             "crystal_size,"
                             "crystal_gap_size,"
                             "crystal_array_size,"
                             "number_of_detector_modules,"
                             "detector_module_axial_offsets,"
                             "DOI_information,"
                             "detector_module_angular_offset,"
                             "number_of_radial_bins,"
                             "sensitivity,"
                             "TOF_information,"
                             "iPSF_model,"
                             "initial_guess,"
                             "reconstruction_output_setting,"
                             "iterative_algorithm_type,"
                             "regularizer_strength,"
                             "regularizer_model_type,"
                             "regularizer_potential_function_type,"
                             "regularizer_neighborhood_properties,"
                             "regularizer_buildin_parameter_list,"
                             "regularizer_spatial_variant_weight,"
                             "warmup_setting,"
                             "iteration_setting,"
                             "input_raw_data_file,"
                             "input_raw_data_format_type,"
                             "kernel_matrix"
                             "gb_alpha,"
                             "gb_CIP");
                             // "input_raw_data_format_type");

//
// pick options that are necessary for a minimal requirement
//
static std::string g_sysopts_mandatory("image_size,"
                                       "voxel_size,"
                                       "detector_ring_diameter,"
                                       "crystal_size,"
                                       "crystal_array_size,"
                                       "number_of_detector_modules");

///
/// \brief Describe a crystal.
/// \see DetBlock
///
typedef struct tagCRYSTALSIZE {
    double transaxial_front; ///< transaxial size on front face (unit in mm).
    double transaxial_back;  ///< transaxial size on back face (unit in mm).
    double axial; ///< axial size (unit in mm).
    double depth; ///< depth size (as well as the length of crystal) (unit in mm).
} CRYSTALSIZE;

///
/// \brief Describe a gap.
///
typedef struct tagCRYSTALGAPSIZE {
    double transaxial_front; ///< transaxial size on front face (unit in mm).
    double transaxial_back; ///< transaxial size on back face (unit in mm).
    double axial; ///< axial size (unit in mm).
} CRYSTALGAPSIZE;

///
/// \brief Describe the size of an often-used array.
///
/// Examples are crystal array, crystal divisions which all use this size.
///
typedef struct tagSIZE {
    union {
        int t;
        int i;
    }; ///< transaxial size (or image dim i).
    union {
        int a;
        int j;
    }; ///< axial size (or image dim j).
    union {
        int d;
        int k;
    }; ///< size in depth direction (or image dim k).
} SIZE;

///
/// \brief Describe DOI information.
///
typedef struct tagDOIINFO {
    double resolution; ///< resolution (usually the FWHM, unit in mm).
    int number_of_bins; ///< number of DOI bins
} DOIINFO;

///
/// \brief Describe TOF information (for TOF recon)
///
typedef struct tagTOFINFO {
    double resolution; ///< resolution (FWHM, unit in peco second)
    double bin_size; ///< sampling spacing (unit in peco second)
} TOFINFO;

///
/// \brief A structure of system parameters
///
typedef struct tagSYSTEMPARAMETERS {
    SIZE image_size; ///< image volume size (i, j, k)
    double voxel_size_i; ///< voxel size i (unit in mm)
    double voxel_size_j; ///< voxel size j
    double voxel_size_k; ///< voxel size k
    double detector_ring_diameter; ///< detector ring diameter (unit in mm)
    SIZE crystal_array_size; ///< crystal array size (t, a)
    CRYSTALSIZE crystal_size; ///< crystal geometry inforamtion
    CRYSTALGAPSIZE crystal_gap_size; ///< gap information
    SIZE number_of_detector_modules; ///< detector block number (t, a)
    double detector_module_axial_offsets[10]; ///< detector block axial offsets (not used)
    DOIINFO doi_info; ///< DOI information
    double detector_module_angular_offset; ///< angle offset (when placing detector blocks on a ring)
    std::string recon_output_folder; ///< folder for saving reconstructions (used by recon)
    std::string recon_output_filename_prefix; ///< name prefix for recon result (used by recon)
    std::string initial_guess; ///< initial image
    std::string sensitivity; ///< sensitivity image (used by TOF LM-recon)
    TOFINFO tof_info; ///< TOF information
    std::string ipsf_transaxial; ///< transaxial image blurring kernel matrix
    std::string ipsf_axial; ///< axial image blurring kernel matrix
    double regularizer_strength;
    int regularizer_model_type;
    int regularizer_potential_function_type;
    int regularizer_neighborhood_size;
    bool regularizer_isotropic_or_not;    
    double regularizer_buildin_parameters[5]; // 5 at most
    int regularizer_number_of_buildin_parameters_detected;
    std::string regularizer_spatial_variant_weight;
    std::string input_raw_data_file;
    int input_raw_data_format_type;
    bool warmup;
    int number_of_subsets_warmup;
    int number_of_iterations_warmup;
    int iterative_algorithm_type;
    int number_of_subsets;    
    int number_of_iterations;
    int stepsize_for_intermediate_result;
    int number_of_radial_bins;

    std::string kernel_matrix; 

    float gb_alpha;
    std::string gb_CIP;

    
} SYSTEMPARAMETERS;

} // end of UCD_PET_SYS

#endif
