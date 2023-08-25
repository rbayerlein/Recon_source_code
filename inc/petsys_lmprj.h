/*!
 * class DrawLine and ListModeProjector
 *
 * Modified based on old implementation.
 *
 * Author: Jian Zhou
 * Date: 11-22-2013
 *
 */
#ifndef PETSYS_LMPRJ_H
#define PETSYS_LMPRJ_H

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <petsys_image.h>
#include <petsys_prj.h>
#include <petsys_log.h>

namespace UCDPETSYS
{

/*
	\class DrawLine
	\brief do LOR-based on-the-fly forward and backprojection
 */
class DrawLine
{
public:
    explicit DrawLine(const std::size_t size_i,
                      const std::size_t size_j,
                      const std::size_t size_k,
                      const float vox_size_i,
                      const float vox_size_j,
                      const float vox_size_k,
                      const float tw_sigma = -1.0,
                      const float tw_spacing = -1.0);
    ~DrawLine();

public:
    // forward projection
    float fpBresenham(float p0x, float p0y, float p0z,
                      float p1x, float p1y, float p1z,
                      const Image<float>& img,
                      float weight, const int tbin_id);

    // backprojection
    void bpBresenham(float p0x, float p0y, float p0z,
                     float p1x, float p1y, float p1z,
                     Image<float>& img,
                     float weight, const int tbin_id);

    // backprojection (to dump the projection matrix, non-tof only)
    int bpBresenham(float p0x, float p0y, float p0z,
                    float p1x, float p1y, float p1z,
                    Image<float>& img,
                    float weight, int* index_buff);

    float fpSiddon(float p0x, float p0y, float p0z,
                   float p1x, float p1y, float p1z,
                   const Image<float>& img,
                   float weight, const int tbin_id);

    void bpSiddon(float p0x, float p0y, float p0z,
                  float p1x, float p1y, float p1z,
                  Image<float>& img,
                  float weight, const int tbin_id);

    void bpSiddonSquare(float p0x, float p0y, float p0z,
                  float p1x, float p1y, float p1z,
                  Image<float>& img,
                  float weight);

    int bpSiddon(float p0x, float p0y, float p0z,
                 float p1x, float p1y, float p1z,
                 Image<float>& img,
                 float weight, int* index_buff);

private:
    // check if ray hits cylinder
    bool hitCylFOV(float& tmin, float& tmax,
                   const float fov_radius,
                   const float p0x, const float dx,
                   const float p1x, const float dy);
    

    // v0: Jian's code for object's length longer than scanner
    // check if ray hits box, (x, y dir only, so it's a box without top and bottom faces)
    bool hitBoxFOV(float& tmin, float& tmax,
                   const float p0x, const float dx,
                   const float p0y, const float dy);


    //Xuezhu's modified code for EXPLORER (object's length maybe shorter than scanner)
    typedef float REAL;
    bool hitCheck(const REAL p0x, const REAL p0y, const REAL p0z,
                  const REAL p1x, const REAL p1y, const REAL p1z,
                  REAL& t_min, REAL& t_max);



    // calculate TOF LOR end points according to the given timing bin idx
    void calcTOFLOREndPoints(const int tbin_id,
                             float& p0x, float& p0y, float& p0z,
                             float& p1x, float& p1y, float& p1z,
                             float& nrm_dx, float& nrm_dy, float& nrm_dz);
//	public:
//		typedef union tagFLTYPE {
//			unsigned long i;
//			double f;
//		}FLTYPE;
public:
    static const float PI;
    static const float FLT_MIN;

private:
    std::size_t m_img_size_i;
    std::size_t m_img_size_j;
    std::size_t m_img_size_k;
    float m_vox_size_i;
    float m_vox_size_j;
    float m_vox_size_k;
    float m_bd_min_x;
    float m_bd_min_y;
    float m_bd_min_z;
    float m_fov_radius;

    // image planes
    std::vector<float> m_x_cuts;
    std::vector<float> m_y_cuts;
    std::vector<float> m_z_cuts;

    // gaussian window lookuptable
    float m_tw_sigma;
    float m_tw_spacing;
    std::vector<float> m_gw_lut;
    float m_gs_inv;

    // constants
    static const int GWSAMPLESIZE;
};

/*
	\class ListModeProjector
	\brief do list-mode-based forward and backprojection, a wrapper of DrawLine
 */
class ListModeProjector : public Projector
{

public:
    ListModeProjector();
    ~ListModeProjector();

public:
    /*
    	\brief list-mode event structure
    	modified and added two terms
     */
    typedef struct tagLMEVENT {
        short xtal_id_1; ///< crystal index for first crystal [0:MAX-1]
        short ring_id_1; ///< ring index for first crystal [0:MAX-1]
        short xtal_id_2; ///< crystal index for second crystal [0:MAX-1]
        short ring_id_2; ///< ring index for second crystal [0:MAX-1]
        short tbin_id; ///< timing bin index [-MAX : +MAX]
//       float mul_fac;
//       float add_fac;
    } LMEVENT;
    /*
    	\brief old structure to describe a crystal
    	this structure will be deprecated soon
     */
    typedef struct tagCRYSTAL {
        float pos[3]; ///< crystal spatial coord (only this information is required)
        float dummy[10];
    } CRYSTAL;

    /*
    	\brief 2d position
     */
    typedef struct tagCRYSTALXY {
        float x;
        float y;
    } CRYSTALXY;

public:
    /*
    	\brief setup a projector, new version, need more parameters
    	to figure out crystal spatial locations (in-plane)
    	for cylindrical scanner only !
     */
    bool initialize(const int num_of_rings,
               		const float ring_pitch,
               		const int num_of_block_t,
               		const int xtal_array_t,
               		const float ring_diameter,
               		const float xtal_pitch_t,
               		const float xtal_length,
               		const int img_size_ij,
               		const int img_size_k,
               		const float vox_size_ij,
               		const float vox_size_k,
               		const float tw_resolution = -1.0f,
               		const float tw_spacing = -1.0f,
               		const bool half_ang = false);

    void doForwardProj(Image<float>& proj,
                       const Image<float>& img,
                       const std::vector<LMEVENT>& lm) const;
    void doBackProj(Image<float>& img,
                    const Image<float>& proj,
                    const std::vector<LMEVENT>& lm) const;

	void initializeSensitivityImage(const char* file, const SIZE& image_size);

	// these two functions extract subset directly from the global buffer
	// save memory space
	void doFProj(Image<float>& proj, const Image<float>& img, int subset_id) const;
	void doBProj(Image<float>& img, const Image<float>& proj, int subset_id) const;

    // read events from file
//    static bool readEvents(const char* filename, std::vector<LMEVENT>& lm);
    virtual void readRawData(const char* raw_data_file);
    virtual void makeSubsets(int number_of_subsets);
    virtual int getNumberOfSubsets() const;
    virtual void doForwardProj(Image<float>& proj, const Image<float>& img, int subset_id) const;
    virtual void doBackProj(Image<float>& img, const Image<float>& proj, int subset_id) const;
    virtual Image<float>* allocateProjection(int subset_id) const;
    virtual const std::vector<float>& getMultiplicativeFactor(int subset_id) const;
    virtual const std::vector<float>& getAdditiveFactor(int subset_id) const;
	virtual const Image<float>* getSensitivity(int subset_id) const;

    // get file length in byte
    static long long getFileLength(const char* filename);

    DrawLine* getRayTracer();
    const std::vector<CRYSTALXY>& getXtalPosXY() const;
    const std::vector<float>& getRingZ() const;

private:
    std::vector<CRYSTALXY> m_xtal_pos_xy;
    std::vector<float> m_ring_z;
    DrawLine* m_proj;

protected:
    std::vector<LMEVENT> m_raw_lmdata_buff;
//    std::vector<std::vector<LMEVENT> > m_lmdata_subsets; // not necessary any more
    size_t m_number_of_subsets;   // int
    size_t m_max_subset_size;   // int 
// make geometric subsets ??
//    std::vector<int> m_angle_index_buff;    

public:
    enum {
        SIDDON = 0,
        BRESENHAM,
    };

    static int raytracer_type;
};

} // end of UCDPETSYS

#endif // petsys_lmprj.h
