/*!
 * class IPSFModel
 *
 * Modified based on old class ImageBlurSep.
 *
 * Author: Jian Zhou
 * Date: 11-22-2013
 *
 */
#ifndef PETSYS_IPSF_H
#define PETSYS_IPSF_H

#include <petsys_spmtx.h>
#include <petsys_sysparms.h>
#include <petsys_log.h>

namespace UCDPETSYS
{
class IPSFModel
{

public:
    IPSFModel(int img_size_i, int img_size_j, int img_size_k);
    IPSFModel(SIZE sz);
    ~IPSFModel();

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

