/*!
 * class DrawLine and ListModeProjector
 *
 * Bass class generalizing all projectors
 *
 * Author: Jian Zhou
 * Date: 11-22-2013
 *
 */
#ifndef PETSYS_PRJ_H
#define PETSYS_PRJ_H

#include <vector>
#include <petsys_image.h>
#include <petsys_log.h>
#include <petsys_ipsf.h>

namespace UCDPETSYS
{

class Projector
{
public:
    Projector() : m_ipsf_model(0), m_sensitivity(0) {};
    ~Projector() { 
    	if (m_ipsf_model != 0) {
    		delete m_ipsf_model; 
    	}
    	if (m_sensitivity != 0) {
    		delete m_sensitivity;
    	}
    };

public:
	enum {LM_PROJ = 0, SINO_PROJ};

public:
    virtual void readRawData(const char* raw_data_file) {};
    virtual void makeSubsets(int number_of_subsets) {};
    virtual int getNumberOfSubsets() const { 
        return 0;
    };
    virtual void doForwardProj(Image<float>&, const Image<float>&, int subset_id) const {};
    virtual void doBackProj(Image<float>&, const Image<float>&, int subset_id) const {};
    virtual Image<float>* allocateProjection(int subset_id) const {
        return 0;
    };
    virtual const std::vector<float>& getMultiplicativeFactor(int subset_id) const {
        return m_mul_fac_subsets[subset_id];
    };
    virtual const std::vector<float>& getAdditiveFactor(int subset_id) const {
        return m_add_fac_subsets[subset_id];
    };
    
    virtual const Image<float>* getSensitivity(int subset_id) const {};

public:
	void initializeIPSFModel(const SIZE& image_size, 
						 	 const char* trans_psf_file, 
						 	 const char* axial_psf_file);

    // void initializeKernelModel(const SIZE& image_size, 
    //                          const char* kernel_matrix_file);  
	
	typedef std::vector<std::pair<int, int> > XtalPair;
	XtalPair createCrystalPairs(const int num_of_detblocks_t,
        					    const int xtal_array_size_t,
        					    const int num_of_radial_bins_per_angle,
        					    const bool compr_enabled,
        					    const bool half_block_rotation); // old function	

protected:

	IPSFModel* m_ipsf_model;

    // IKernelModel* m_ikernel_model;

	Image<float>* m_sensitivity;
    std::vector<float> m_mul_fac_buff;
    std::vector<float> m_add_fac_buff;
    std::vector<std::vector<float> > m_mul_fac_subsets;
    std::vector<std::vector<float> > m_add_fac_subsets;

public:
    static int NUMBER_OF_THREADS_FP;
    static int NUMBER_OF_THREADS_BP;    
};

} // end of UCDPETSYS

#endif // petsys_prj.h
