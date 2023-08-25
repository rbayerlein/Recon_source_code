/*!
 * class Configurator
 *
 * modified from old class Configurator
 *
 *
 * Author: Jian Zhou
 * Date: 11-21-2013
 *
 */

#ifndef PETSYS_CSCANNER_CUT_H
#define PETSYS_CSCANNER_CUT_H

#include <petsys_cfg.h>
#include <petsys_lmprj.h>
#include <petsys_ipsf.h>
#include <petsys_reg.h>
#include <petsys_timer.h>
#include <petsys_recon.h>

namespace UCDPETSYS
{

class CScanner
{

public:
    CScanner();
    ~CScanner();

public:
    void initialize(const char* config_file);
	void runRecon();
    void runRecon2();
    void runRecon2_cut();
    void runfplinesection(); //

public:
	enum {
		EM = 0, 
		ADMM
	};

public:
    int getNumberOfCrystalRings() const;
    float getCrystalRingPitchSize() const;
    float getCrystalTransaxialPitchSize() const;
    int getNumberOfRadialBins() const;
    const SYSTEMPARAMETERS& getSystemParameters() const;

protected:	
	Projector*		createProjector();
	Regularizer* 	createRegularizer();

//    int doEMBasedPLRecon(Image<float>& output,
//                         const int number_of_iterations,
//                         const int stepsize_for_intermediate_result,
//                         const float beta);

//    int doEMBasedPLReconWithMeanBasedRegularizer(Image<float>& x,
//							   					 const Image<float>& mean_image,
//                               					 const int number_of_iterations,
//                               					 const int stepsize_for_intermediate_result,
//                               					 const float beta);

//	int doADMMBasedPLRecon(Image<float>& output,
 //                          const int number_of_iterations,
   //                        const int stepsize_for_intermediate_result,
     //                      const float beta);

private:
    SYSTEMPARAMETERS 	m_sysparms;
    Projector* 			m_projector;
    Regularizer* 		m_regularizer;

};


}; // end of UCDPETSYS


#endif // petsys_cscanner.h
