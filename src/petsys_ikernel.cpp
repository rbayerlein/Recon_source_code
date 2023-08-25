#include <petsys_ikernel.h>

IKernelModel::IKernelModel(int img_size_i, int img_size_j, int img_size_k) :
    m_img_size_i(img_size_i),
    m_img_size_j(img_size_j),
    m_img_size_k(img_size_k),
    m_blur_t(0), m_blur_tt(0),
    m_blur_a(0), m_blur_at(0)
{
    m_img_tmp = new Image<float>(img_size_i, img_size_j, img_size_k);
}

IKernelModel::IKernelModel(SIZE sz) :
    m_img_size_i(sz.i),
    m_img_size_j(sz.j),
    m_img_size_k(sz.k),
    m_blur_t(0), m_blur_tt(0),
    m_blur_a(0), m_blur_at(0)
{
    m_img_tmp = new Image<float>(sz);
}


IKernelModel::~IKernelModel()
{
    if (!m_blur_t) {
        delete m_blur_t;
    }

    if (!m_blur_a) {
        delete m_blur_a;
    }

    if (!m_blur_tt) {
        delete m_blur_tt;
    }

    if (!m_blur_at) {
        delete m_blur_at;
    }

    if (!m_img_tmp) {
        delete m_img_tmp;
    }
}

bool IKernelModel::initialize(const char* kernel_matrix_filename)
{
	SystemLog::write("initializing iKernel model ...\n");

    // setup transaxial blurring matrix
    m_blur_t = new PSpMatrix();

    if (!m_blur_t->read(bt_filename)) {
        delete m_blur_t;
        m_blur_t = 0;
        //return false;
        m_blur_tt = 0;
    } else {
        m_blur_tt = m_blur_t->getTranspose();
    }

    // setup axial blurring matrix
    m_blur_a = new PSpMatrix();

    if (!m_blur_a->read(ba_filename)) {
        delete m_blur_a;
        m_blur_a = 0;
        return false;
    }

    m_blur_at = m_blur_a->getTranspose();

    return true;
}

void IKernelModel::blur(Image<float>& img_blur, const Image<float>& img, BLUR_type type) const
{
    int ni = img.getDimI();
    int nj = img.getDimJ();
    int nk = img.getDimK();

    if (ni != m_img_size_i ||
            nj != m_img_size_j ||
            nk != m_img_size_k) {
#if DEBUG
        std::printf("image blur failed: dimension mismatch");
#endif
        abort();
    }

    int nv_per_plane = ni * nj;
    int k, offset = 0;

    if (type == F_BLURRING) {
        if (m_blur_t) {
            // t-blur (always use OMP for acceleration)
            #pragma omp parallel for private(k, offset)
            for (k = 0; k < nk; k ++) {
                offset = k * nv_per_plane;
                m_blur_t->proj(img_blur.getPtr() + offset,
                               img.getPtr() + offset);
            }
        } else {
            // just pass without processing
            memcpy(img_blur.getPtr(), img.getPtr(), sizeof(float)*img.getSize());
        }

        if (!m_blur_a) {
            return;
        }

        // make a backup
        memcpy(m_img_tmp->getPtr(),
               img_blur.getPtr(),
               sizeof(float)*nv_per_plane * nk);
        img_blur.clear();
        // a-blur
        #pragma omp parallel for private(k)

        for (k = 0; k < nk; k ++) {
            blurAxialK(img_blur, m_blur_a, *m_img_tmp, k);
        }

    } else {
        if (m_blur_tt) {
            // t-blur (always use OMP for acceleration)
            #pragma omp parallel for private(k, offset)
            for (k = 0; k < nk; k ++) {
                offset = k * nv_per_plane;
                m_blur_tt->proj(img_blur.getPtr() + offset,
                                img.getPtr() + offset);
            }
        } else {
            // just pass without processing
            memcpy(img_blur.getPtr(), img.getPtr(), sizeof(float)*img.getSize());
        }

        if (!m_blur_at) {
            return;
        }

        // make a backup
        memcpy(m_img_tmp->getPtr(),
               img_blur.getPtr(),
               sizeof(float)*nv_per_plane * nk);
        img_blur.clear();
        // a-blur
        #pragma omp parallel for private(k)

        for (k = 0; k < nk; k ++) {
            blurAxialK(img_blur, m_blur_at, *m_img_tmp, k);
        }
    }
}

void IKernelModel::blurAxialK(Image<float>& img_blur,
                           const PSpMatrix* blur_a,
                           const Image<float>& img,
                           const int k) const
{
    int ni = img.getDimI();
    int nj = img.getDimJ();

    for (int l = blur_a->getRowPtr()[k]; l < blur_a->getRowPtr()[k + 1]; l++) {
        int c = blur_a->getCol()[l] - 1;
        float v = blur_a->getVal()[l];

        for (int j = 0; j < nj; j ++) {
            for (int i = 0; i < ni; i ++) {
                img_blur(i, j, k) += img(i, j, c) * v;
            }
        }
    }
}
