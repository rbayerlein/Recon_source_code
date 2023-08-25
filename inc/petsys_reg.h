/*!
 * class Regularizer
 *
 *
 * Author: Jian Zhou
 * Date: 11-22-2013
 *
 */
#ifndef PETSYS_REG_H
#define PETSYS_REG_H

#include <petsys_image.h>
#include <petsys_log.h>
#include <cmath>
#include <vector>


namespace UCDPETSYS
{
class Regularizer
{

public:
    Regularizer();
    virtual ~Regularizer();

public:
	enum {
		PAIRWISE_MRF = 0,
		PATCH,
	};

public:
    virtual void initialize(int potential_function_type,
    		   				int neighborhood_size,
    		   				bool is_isotropic,
    		   				double* buildin_parms,
    		   				int number_of_buildin_parms_detected);
    		   		
    virtual void calculateGradientBasedOnOptimizationTransfer(const Image<float>& x,
            												 Image<float>& term1,
            												 Image<float>& term2);

    // virtual void gbcalg(const Image<float>& x,
    //         												 Image<float>& term1,
    //         												 Image<float>& term2,
    //                                                          const Image<float>& cip,
    //                                                          Image<int>* fd);

	virtual void calculateGradient(const Image<float>& x, Image<float>& grad);	
	
public:	
	virtual void setSpatialVariantWeight(const Image<float>& weight);

protected:
	enum POTENTIAL_FUNCTION_type{
		QUAD = 0,
		HYPERBOLA, // or TV
		FAIR,// or Lange
		HU, // Huber
		HL, // nonconvex, Hebert & Leahy
		HW, // nonconvex, Holland & Welsch
	};
	
	int m_potential_function_type;
	int m_neighborhood_size;
	bool m_is_isotropic;
	double* m_buildin_parms;
	Image<float>* m_spatial_variant_weight;


protected:
	double phi(double t);
	double dot_phi_over_t(double t);
	double dot_phi(double t);

public:
	static int NUMBER_OF_THREADS;
		
};

//
// class MeanBasedRegularizer
// to handle a simple prior of the form: \sum_j|x_j - x_mean_j|^2
//
class MeanBasedRegularizer : public Regularizer
{
public:
	MeanBasedRegularizer();
	~MeanBasedRegularizer();
	
public:
	void initialize(int potential_function_type,
    		   		int neighborhood_size,
    		   		bool is_isotropic,
    		   		double* buildin_parms,
    		   		int number_of_buildin_parms_detected);

	void setMeanImage(const Image<float>& mean_image);
  
	void calculateGradientBasedOnOptimizationTransfer(const Image<float>& x,
            										 Image<float>& term1,
            										 Image<float>& term2);   

    // //     @Guobao. similar function to calculate global penalty (gradient)
    // void gbcalg(const Image<float>& x,
    //         										   Image<float>& term1,
    //         										   Image<float>& term2,
    //                                                    const Image<float>& cip,
    //                                                    Image<int>* fd);

	void calculateGradient(const Image<float>& x, Image<float>& grad) {};	

private:
	Image<float>* m_mean_image;
};

//
// class PairwiseMRFRegularizer
//
class PairwiseMRFRegularizer : public Regularizer
{

public:
	PairwiseMRFRegularizer();
	~PairwiseMRFRegularizer();

public:
	void initialize(int potential_function_type,
    		   		int neighborhood_size,
    		   		bool is_isotropic,
    		   		double* buildin_parms,
    		   		int number_of_buildin_parms_detected);
		
	void calculateGradientBasedOnOptimizationTransfer(const Image<float>& x,
            										 Image<float>& term1,
            										 Image<float>& term2); // gradient of surrograte

    //     @Guobao. similar function to calculate global penalty (gradient)
    void gbcalg(const Image<float>& x,
            										   Image<float>& term1,
            										   Image<float>& term2,
                                                       const Image<float>& cip,
                                                       Image<int>* fd);    

	void calculateGradient(const Image<float>& x, Image<float>& grad); // gradient of original function

private:
//	void calculateNeighborhoodDifference(const Image<float>& x, 
//										 const Image<float>& neighbor_location, 
//										 Image<float>& diff);
	
	void createNeighborVoxelPositionOffset();
	void calculateTerm1AndTerm2Anisotropic(const Image<float>& x,
										   const Image<float>& neighbor_voxel_position_offset,
										   const float geometrical_weight,
										   Image<float>& term1,
										   Image<float>& term2);
	Image<float>* calculateWeightIsotropic(const Image<float>& x);
	
	void calculateTerm1AndTerm2Isotropic(const Image<float>& x,
										 const Image<float>& neighbor_voxel_position_offset,
										 const Image<float>& weight,
										 const float geometrical_weight,
										 Image<float>& term1,
										 Image<float>& term2);							  
	
	// for gradient (of orignal penalty function)
	void calculateGradientAnisotropic(const Image<float>& x,
									  const Image<float>& neighbor_voxel_position_offset,
									  const float geometrical_weight,
									  Image<float>& grad_k);
									  
	void calculateGradientIsotropic(const Image<float>& x,
									const Image<float>& neighbor_voxel_position_offset,
								    const Image<float>& weight,
									const float geometrical_weight,
									Image<float>& grad_k);

private:
	std::vector<Image<float>* > m_neighbor_voxel_position_offset;
	std::vector<float> m_geometrical_weight;
};

//
// class PatchBasedRegularizer
//
class PatchBasedRegularizer : public Regularizer
{
public:
	PatchBasedRegularizer();
	~PatchBasedRegularizer();

public:
	void initialize(int potential_function_type,
    		   		int neighborhood_size,
    		   		bool is_isotropic,
    		   		double* buildin_parms,
    		   		int number_of_buildin_parms_detected);

	void calculateGradientBasedOnOptimizationTransfer(const Image<float>& x,
            										 Image<float>& term1,
            										 Image<float>& term2);

//     @Guobao. similar function to calculate global penalty (gradient)
    void gbcalg(const Image<float>& x,
            										   Image<float>& term1,
            										   Image<float>& term2,
                                                       const Image<float>& cip,
                                                       Image<int>* fd);    

	void calculateGradient(const Image<float>& x, Image<float>& grad);

private:
	// added to compute isotropic form (jnzhou@04-21-2014)
	Image<float>* calculateWeightIsotropic(const Image<float>& x);

	float calculateWeightIsotropicSingleVoxel(const Image<float>& x, const int voxel_index);

	void createNeighborPatchPositionOffset();
	void calculateTerm1AndTerm2Anisotropic(const Image<float>& x,
										   Image<float>& patch_position_offset,
										   Image<float>& term1,
										   Image<float>& term2);											   
	void calculateGradientAnisotropic(const Image<float>& x,
									  Image<float>& patch_poisition_offset,
									  Image<float>& grad_k);

	void calculateTerm1AndTerm2Isotropic(const Image<float>& x,
										 const Image<float>& patch_position_offset,
										 const Image<float>& weight,
										 Image<float>& term1,
										 Image<float>& term2);

private:
	int m_patch_size;
	std::vector<Image<float>* > m_neighbor_patch_position_offset;
};

}// end of UCDPETSYS

#endif // petsys_reg.h

