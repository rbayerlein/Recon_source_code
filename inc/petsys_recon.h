#ifndef PETSYS_RECON_H
#define PETSYS_RECON_H

#include <petsys_image.h>
#include <petsys_prj.h>
#include <petsys_reg.h>
#include <petsys_timer.h>

//
// EM
//
class EMBasedPLReconstructor 
{
public:
	static void run(Projector* projector,
			   		Regularizer* regularizer,
               		const int number_of_iterations,
               		const int stepsize_for_intermediate_result,
               		const float beta, 
               		const char* recon_output_folder,
               		const char* recon_output_file_prefix,
               		Image<float>& output);
};

//
// ADMM
//
class ADMMBasedPLReconstructor
{
public:
	static float mu;
	static void run(Projector* projector,
			   		Regularizer* regularizer,
               		const int number_of_iterations,
               		const int stepsize_for_intermediate_result,
               		const float beta, 
               		const char* recon_output_folder,
               		const char* recon_output_file_prefix,
               		Image<float>& output);
};

#endif
