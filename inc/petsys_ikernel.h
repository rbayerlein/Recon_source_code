/*!
 * class IKernelModel
 *
 * Modified based on old class ImageBlurSep.
 *
 * Author: Xuezhu Zhang
 * Date: 2014~2018
 *
 */
#ifndef PETSYS_IKERNEL_H
#define PETSYS_IKERNEL_H

#include <petsys_spmtx.h>
#include <petsys_sysparms.h>
#include <petsys_log.h>

namespace UCDPETSYS
{
class IKernelModel
{

public:
    IKernelModel(int img_size_i, int img_size_j, int img_size_k);
    IKernelModel(SIZE sz);
    ~IKernelModel();

public:
    enum BLUR_type {
        F_BLURRING = 0,
        B_BLURRING,
    };

public:
    bool initialize(const char* bt_filename, const char* ba_filename);
    void blur(Image<float>& img_blur, const Image<float>& img, BLUR_type type) const;


private:
    /*
    	\brief this function is for the easy use of openmp
     */
    void blurAxialK(Image<float>& img_blur,
                    const PSpMatrix* blur_a,
                    const Image<float>& img,
                    const int k) const;

private:
    int m_img_size_i;
    int m_img_size_j;
    int m_img_size_k;
    Image<float>* m_img_tmp;
    PSpMatrix* m_blur_t;
    PSpMatrix* m_blur_tt; // transpose of m_blur_t
    PSpMatrix* m_blur_a;
    PSpMatrix* m_blur_at; // transpose of m_blur_a

};

}
#endif // petsys_ipsf.h

