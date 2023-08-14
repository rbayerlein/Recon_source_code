#include <petsys_spmtx.h>

SSpMatrix::SSpMatrix(const int row_num, const int col_num, const long long nnz,
                     const float scaling_factor, const int* row_ptr, const PACKAGE* p) :
    m_row_num(row_num),
    m_col_num(col_num),
    m_nnz(nnz),
    m_scaling_factor(scaling_factor)
{
    if (row_num > 0) {
        m_row_ptr = new int[row_num + 1];
        memcpy(m_row_ptr, row_ptr, sizeof(int) * (row_num + 1));
    } else {
        m_row_ptr = 0;
    }

    if (nnz > 0) {
        m_p = new PACKAGE[nnz];
        memcpy(m_p, p, sizeof(PACKAGE)*nnz);
    } else {
        m_p = 0;
    }
}

SSpMatrix::SSpMatrix(const PSpMatrix& pmat, const SIZE& imgsize)
{
    m_row_num = pmat.getRowNum();
    m_col_num = pmat.getColNum();
    m_nnz = pmat.getNNZ();

    if (m_nnz == 0) {
        m_scaling_factor = 0.0;
        m_row_ptr = 0;
        m_p = 0;
    } else {
        // copy row ptr
        m_row_ptr = new int[m_row_num + 1];
        memcpy(m_row_ptr, pmat.getRowPtr(), sizeof(int) * (m_row_num + 1));

        //
        float max_val = 0.0;

        for (int i = 0; i < m_nnz; i ++) {
            if (max_val < pmat.getVal()[i]) {
                max_val = pmat.getVal()[i];
            }
        }

        //
#if USE_9BIT // use 9-bit to compress nonzero values including its index
        std::printf("use 9-bit compression ... ");
        m_scaling_factor = max_val / 511.0;
#else // else use 17-bit compression (up to 6 bytes for each PACKAGE)
        std::printf("use 17-bit compression ... ");
        m_scaling_factor = max_val / 131071.0;
#endif
        //
#if DEBUG
        std::printf("max nonzero: %e, scaling_factor: %e ... ", max_val, m_scaling_factor);
#endif
        // convert index to PACKAGE including nonzero values
        int img_dim_ij = imgsize.i * imgsize.j;

        m_p = new PACKAGE[m_nnz];

        for (int i = 0; i < m_nnz; i ++) {
            int idx = pmat.getCol()[i] - 1;

            int kk = (idx) / img_dim_ij;
            int rem = idx - kk * img_dim_ij;
            int ii = rem % imgsize.i;
            int jj = rem / imgsize.i;

            m_p[i].i = ii;
            m_p[i].j = jj;
            m_p[i].k = kk;
            m_p[i].w = int(pmat.getVal()[i] / m_scaling_factor);
        }

#if DEBUG
        std::printf("OK!\n");
#endif
    }
}

SSpMatrix::~SSpMatrix()
{
    if (m_row_ptr != 0) {
        delete [] m_row_ptr;
    }

    if (m_p != 0) {
        delete [] m_p;
    }
}

void SSpMatrix::read(const char* filename)
{
    std::ifstream input(filename);

    if (!input) { // file does not exist
        abort();
    } else {

        // follow the format of SMAT
        // 1. number of rows of a SMAT
        input.read((char*)&m_row_num, sizeof(int));

        if (input.fail()) {
#if DEBUG
            std::printf("bad sparse matrix format");
#endif
            abort();
        }

        // 2. number of columns
        input.read((char*)&m_col_num, sizeof(int));

        if (input.fail()) {
#if DEBUG
            std::printf("bad sparse matrix format");
#endif
            abort();
        }


        // 3. number of nonzeros
        input.read((char*)&m_nnz, sizeof(int));

        if (input.fail()) {
#if DEBUG
            std::printf("bad sparse matrix format");
#endif
            abort();
        }

        // 4. scaling factor
        input.read((char*)&m_scaling_factor, sizeof(float));

        if (input.fail()) {
#if DEBUG
            std::printf("bad sparse matrix format");
#endif
            abort();
        }

        // 5. row offset table
        if (m_row_ptr != 0)
            delete [] m_row_ptr;

        m_row_ptr = new int[m_row_num + 1];

        if (m_row_ptr == 0) {
#if DEBUG
            std::printf("allocate row_ptr failed");
#endif
            abort();
        }

        input.read((char*)m_row_ptr, sizeof(int) * (m_row_num + 1));

        if (input.fail()) {
#if DEBUG
            std::printf("bad sparse matrix format");
#endif
            abort();
        }

        // 6. voxel indices and (quantized) nonzero values
        if (m_p != 0)
            delete [] m_p;

        m_p = new PACKAGE[m_nnz];

        if (m_p == 0) {
#if DEBUG
            std::printf("allocate p failed");
#endif
            abort();
        }

        input.read((char*)m_p, sizeof(PACKAGE)*m_nnz);

        if (input.fail()) {
#if DEBUG
            std::printf("bad sparse matrix format");
#endif
            abort();
        }
    }
}

void SSpMatrix::write(const char* smat_filename)
{
    std::ofstream output(smat_filename);

    if (!output) { // no disk space, or wrong device
        abort();
    } else {

        output.write((char*)&m_row_num, sizeof(int));

        if (output.fail()) {
#if DEBUG
            std::printf("write row number failed.");
#endif
            abort();
        }

        output.write((char*)&m_col_num, sizeof(int));

        if (output.fail()) {
#if DEBUG
            std::printf("write col number failed.");
#endif
            abort();
        }

        output.write((char*)&m_nnz, sizeof(int));

        if (output.fail()) {
#if DEBUG
            std::printf("write nonzero number failed.");
#endif
            abort();
        }

        output.write((char*)&m_scaling_factor, sizeof(float));

        if (output.fail()) {
#if DEBUG
            std::printf("write scaling factor failed.");
#endif
            abort();
        }

        output.write((char*)m_row_ptr, sizeof(int) * (m_row_num + 1));

        if (output.fail()) {
#if DEBUG
            std::printf("write row ptr failed.");
#endif
            abort();
        }

        output.write((char*)m_p, sizeof(PACKAGE) * m_nnz);

        if (output.fail()) {
#if DEBUG
            std::printf("write package failed.");
#endif
            abort();
        }
    }

}

void SSpMatrix::proj(Image<float>& y, Image<float>& x, PROJ_type t)
{
    if (t == FORWARD) {
        PACKAGE p;

        for (int i = 0; i < m_row_num; i ++) {
            for (int j = m_row_ptr[i]; j < m_row_ptr[i + 1]; j ++) {
                p = m_p[j];
                y[i] += x[x.sub2ind(p.i, p.j, p.k)] * p.w;
            }

            y[i] *= m_scaling_factor;
        }
    }

    if (t == BACKWARD) {
        PACKAGE p;
        float y0;

        for (int i = 0; i < m_row_num; i ++) {
            y0 = y[i] * m_scaling_factor;

            for (int j = m_row_ptr[i]; j < m_row_ptr[i + 1]; j ++) {
                p = m_p[j];
                x[x.sub2ind(p.i, p.j, p.k)] += y0 * p.w;
            }
        }
    }
}

SSpMatrix* SSpMatrix::pickRows(const std::vector<int>& row_idx)
{
    if (row_idx.size() > 0) {
        if (m_p != 0 || m_row_ptr != 0) {
            int col_num = m_col_num;
            float scaling_factor = m_scaling_factor;
            int row_num = row_idx.size();
            int* row_nz = new int[row_num];
            int* row_ptr = new int[row_num + 1];
            int nnz = 0;
            row_ptr[0] = 0;

            for (size_t n = 0; n < row_idx.size(); n ++) {
                row_nz[n] = m_row_ptr[row_idx[n] + 1] - m_row_ptr[row_idx[n]];
                nnz += row_nz[n];
                row_ptr[n + 1] = row_ptr[n] + row_nz[n];
            }

            int offset = 0;
            PACKAGE* p = new PACKAGE[nnz];

            for (size_t n = 0; n < row_idx.size(); n ++) {
                memcpy(&p[offset], m_p + m_row_ptr[row_idx[n]], sizeof(PACKAGE)*row_nz[n]);
                offset += row_nz[n];
            }

            SSpMatrix* smat = new SSpMatrix(row_num, col_num, nnz, scaling_factor, row_ptr, p);

            if (smat != 0) {
                delete [] row_nz;
                delete [] row_ptr;
                delete [] p;
            }

            return smat;
        }

        return 0;
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// PSpMatrix (as PMAT)
//
////////////////////////////////////////////////////////////////////////////////////////////////////
PSpMatrix::PSpMatrix(const int row_num, const int col_num, const long long nnz,
                     const int* row_ptr, const int* col, const float* val) :
    m_row_num(row_num), m_col_num(col_num), m_nnz(nnz)
{
    if (row_num > 0 && row_ptr != 0) {
        m_row_ptr = new int[row_num + 1];
        memcpy(m_row_ptr, row_ptr, sizeof(int) * (row_num + 1));
    } else {
        m_row_ptr = 0;
    }

    if (nnz > 0 && col != 0) {
        m_col = new int[nnz];
        memcpy(m_col, col, sizeof(int)*nnz);
    } else {
        m_col = 0;
    }

    if (nnz > 0 && val != 0) {
        m_val = new float[nnz];
        memcpy(m_val, val, sizeof(float)*nnz);
    } else {
        m_val = 0;
    }
}

PSpMatrix::~PSpMatrix()
{
    //std::printf("release PSpMatrix ...\n");
    if (m_row_ptr != 0) {
        delete [] m_row_ptr;
    }

    if (m_col != 0) {
        delete [] m_col;
    }

    if (m_val != 0) {
        delete [] m_val;
    }
}

bool PSpMatrix::read(const char* filename)
{
    std::ifstream input(filename);

    if (!input) { // file does not exist
        return false;
    } else {

        // row #
        input.read((char*)&m_row_num, sizeof(int));

        if (input.fail()) {
#if DEBUG
            std::printf("read pmat row number failed");
#endif
            return false;
        }

        // col #
        input.read((char*)&m_col_num, sizeof(int));

        if (input.fail()) {
#if DEBUG
            std::printf("read pmat column number failed");
#endif
            return false;
        }

        // nz #
        input.read((char*)&m_nnz, sizeof(long long));

        if (input.fail()) {
#if DEBUG
            std::printf("read pmat nnz failed");
#endif
            return false;
        }

        // row ptr (nnz of each row)
        int* rptr = new int[m_row_num];

        if (!rptr) {
#if DEBUG
            std::printf("allocate row ptr failed");
#endif
            return false;
        }

        input.read((char*)rptr, sizeof(int)*m_row_num);

        if (input.fail()) {
#if DEBUG
            std::printf("read row ptr failed");
#endif
            return false;
        }

        // convert to nnz offset
        if (m_row_ptr) {
            delete [] m_row_ptr;
        }

        m_row_ptr = new int[m_row_num + 1];

        if (!m_row_ptr) {
#if DEBUG
            std::printf("allocate row ptr failed!!");
#endif
            return false;
        }

        m_row_ptr[0] = 0;

        for (int i = 1; i <= m_row_num; i ++) {
            m_row_ptr[i] = m_row_ptr[i - 1] + rptr[i - 1];
        }

        // col buff
        if (m_col) {
            delete [] m_col;
        }

        m_col = new int[m_nnz];

        if (!m_col) {
#if DEBUG
            std::printf("allocate col buff failed");
#endif
            return false;
        }

        input.read((char*)m_col, sizeof(int)*m_nnz);

        if (input.fail()) {
#if DEBUG
            std::printf("read col buff failed");
#endif
            return false;
        }

        // val buff
        if (m_val) {
            delete [] m_val;
        }

        m_val = new float[m_nnz];

        if (!m_val) {
#if DEBUG
            std::printf("allocate val buff failed");
#endif
            return false;
        }

        input.read((char*)m_val, sizeof(float)*m_nnz);

        if (input.fail()) {
#if DEBUG
            std::printf("read val buff failed");
#endif
            return false;
        }

        return true;
    }
}

bool PSpMatrix::write(const char* filename)
{
    std::ofstream output(filename);

    if (!output) { // no disk space, or wrong device
        return false;
    } else {

        // dimensions
        output.write((char*)&m_row_num, sizeof(int));

        if (output.fail()) {
#if DEBUG
            std::printf("write number of rows failed.");
#endif
            return false;
        }

        output.write((char*)&m_col_num, sizeof(int));

        if (output.fail()) {
#if DEBUG
            std::printf("write number of columns failed.");
#endif
            return false;
        }

        // nnz
        output.write((char*)&m_nnz, sizeof(long long));

        if (output.fail()) {
#if DEBUG
            std::printf("write number of nonzeros failed.");
#endif
            return false;
        }

        // row ptr
        int* rptr = new int[m_row_num];

        for (int n = 0; n < m_row_num; n ++) {
            rptr[n] = m_row_ptr[n + 1] - m_row_ptr[n];
        }

        output.write((char*)rptr, sizeof(int) * (m_row_num));
        delete [] rptr;

        if (output.fail()) {
#if DEBUG
            std::printf("write row ptr failed.");
#endif
            return false;
        }

        // col
        int* col = new int[m_nnz];

        for (int n = 0; n < m_nnz; n ++) {
            col[n] = m_col[n]; // because of old version
        }

        output.write((char*)col, sizeof(int) * m_nnz);

        if (output.fail()) {
#if DEBUG
            std::printf("write column index buffer failed.");
#endif
            return false;
        }

        // val
        output.write((char*)m_val, sizeof(float) * m_nnz);

        if (output.fail()) {
#if DEBUG
            std::printf("write nonzero value buffer failed.");
#endif
            return false;
        }

        return true;
    }
}

void PSpMatrix::proj(float* y, float* x, SSpMatrix::PROJ_type t)
{
    if (t == SSpMatrix::FORWARD) {
        for (int i = 0; i < m_row_num; i ++) {
            for (int j = m_row_ptr[i]; j < m_row_ptr[i + 1]; j ++) {
                y[i] += x[m_col[j] - 1] * m_val[j];
            }
        }

        return;
    }

    if (t == SSpMatrix::BACKWARD) {
        for (int i = 0; i < m_row_num; i ++) {
            for (int j = m_row_ptr[i]; j < m_row_ptr[i + 1]; j ++) {
                y[m_col[j] - 1] += x[i] * m_val[j];
            }
        }

        return;
    }
}

int PSpMatrix::getRowNum() const
{
    return m_row_num;
}

int PSpMatrix::getColNum() const
{
    return m_col_num;
}

long long PSpMatrix::getNNZ() const
{
    return m_nnz;
}

const int* PSpMatrix::getRowPtr() const
{
    return m_row_ptr;
}

const int* PSpMatrix::getCol() const
{
    return m_col;
}

const float* PSpMatrix::getVal() const
{
    return m_val;
}

PSpMatrix* PSpMatrix::getTranspose()
{
    // don't need to calc if the current sparse matrix is invalid
    if (m_row_num == 0 || m_col_num == 0 ||
            m_nnz == 0 || m_row_ptr == 0 || m_col == 0 || m_val == 0) {
        return 0;
    }

    // convert it to an IJK sparse matrix
    IJK* ijk_vec = new IJK[m_nnz];
    int c = 0;

    for (int i = 0; i < m_row_num; i ++) {
        IJK ijk;
        ijk.i = i;

        for (int j = m_row_ptr[i]; j < m_row_ptr[i + 1]; j ++) {
            ijk.j = m_col[j] - 1;
            ijk.k = m_val[j];
            ijk_vec[c] = ijk;
            c ++;
        }
    }

    // histogram on column index
    std::vector<std::vector<IJK> > ijk_col_major(m_col_num);

    for (long long i = 0; i < m_nnz; i ++) {
        IJK ijk = ijk_vec[i];
        ijk_col_major[ijk.j].push_back(ijk);
    }

    // recalc parameters required by a CSR sparse matrix
    int row_num_t = m_col_num;
    int col_num_t = m_row_num;
    long long nnz_t = m_nnz;
    int* row_ptr_t = new int[row_num_t + 1];
    row_ptr_t[0] = 0;

    for (int i = 1; i <= row_num_t; i ++) {
        row_ptr_t[i] = row_ptr_t[i - 1] + ijk_col_major[i - 1].size();
    }

    c = 0;
    int* col_t = new int[nnz_t];
    float* val_t = new float[nnz_t];

    for (size_t i = 0; i < ijk_col_major.size(); i ++) {
        for (size_t j = 0; j < ijk_col_major[i].size(); j ++) {
            col_t[c] = ijk_col_major[i][j].i + 1; // old format! 1based
            val_t[c] = ijk_col_major[i][j].k;
            c ++;
        }
    }

    // return it to outside
    PSpMatrix* pmat = new PSpMatrix(row_num_t, col_num_t, nnz_t,
                                    row_ptr_t, col_t, val_t);

    if (pmat != 0) {
        delete [] ijk_vec;
        delete [] row_ptr_t;
        delete [] col_t;
        delete [] val_t;
    }

    return pmat;
}

PSpMatrix* PSpMatrix::pickRows(const std::vector<int>& row_idx)
{
    if (row_idx.size() > 0) {
        if (m_col != 0 || m_val != 0 || m_row_ptr != 0) {
            int col_num = m_col_num;
            int row_num = row_idx.size();
            int* row_nz = new int[row_num];
            int* row_ptr = new int[row_num + 1];
            int nnz = 0;
            row_ptr[0] = 0;

            for (size_t n = 0; n < row_idx.size(); n ++) {
                row_nz[n] = m_row_ptr[row_idx[n] + 1] - m_row_ptr[row_idx[n]];
                nnz += row_nz[n];
                row_ptr[n + 1] = row_ptr[n] + row_nz[n];
            }

            int offset = 0;
            int* col = new int[nnz];
            float* val = new float[nnz];

            for (size_t n = 0; n < row_idx.size(); n ++) {
                memcpy(&col[offset], m_col + m_row_ptr[row_idx[n]], sizeof(int)*row_nz[n]);
                memcpy(&val[offset], m_val + m_row_ptr[row_idx[n]], sizeof(float)*row_nz[n]);
                offset += row_nz[n];
            }

            PSpMatrix* pmat = new PSpMatrix(row_num, col_num, nnz, row_ptr, col, val);

            if (pmat != 0) {
                delete [] row_nz;
                delete [] row_ptr;
                delete [] col;
                delete [] val;
            }

            return pmat;
        }

        return 0;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// PSpMatrixForBackProj
//
///////////////////////////////////////////////////////////////////////////////////////////////////
PSpMatrixForBackProj::PSpMatrixForBackProj(const PSpMatrix& pmat) :
    m_row_num_orig(0),
    m_nz_row_num(0),
    m_col_num(0),
    m_nnz(0),
    m_nz_row_mapping(0),
    m_nz_row_ptr(0),
    m_col(0),
    m_val(0)
{
    if (pmat.getNNZ() > 0) {
        m_row_num_orig = pmat.getRowNum();
        m_col_num = pmat.getColNum();
        m_nnz = pmat.getNNZ();
        m_col = new int[m_nnz];
        memcpy(m_col, pmat.getCol(), sizeof(int)*m_nnz);
        m_val = new float[m_nnz];
        memcpy(m_val, pmat.getVal(), sizeof(float)*m_nnz);

        std::vector<int> nz_row_ind;
        std::vector<int> row_nz;

        for (int n = 0; n < pmat.getRowNum(); n ++) {
            int nz = pmat.getRowPtr()[n + 1] - pmat.getRowPtr()[n];

            if (nz > 0) {
                nz_row_ind.push_back(n);
                row_nz.push_back(nz);
            }
        }

        m_nz_row_num = nz_row_ind.size();
        m_nz_row_mapping = new int[nz_row_ind.size()];
        memcpy(m_nz_row_mapping, &nz_row_ind[0], sizeof(int)*nz_row_ind.size());
        m_nz_row_ptr = new int[row_nz.size() + 1];
        m_nz_row_ptr[0] = 0;

        for (size_t n = 1; n < row_nz.size() + 1; n ++) {
            m_nz_row_ptr[n] = m_nz_row_ptr[n - 1] + row_nz[n - 1];
        }
    }
}

PSpMatrixForBackProj::~PSpMatrixForBackProj()
{
    //std::printf("release PSpMatrixForBackProj ...\n");
    if (m_nz_row_mapping != 0) {
        delete [] m_nz_row_mapping;
    }

    if (m_nz_row_ptr != 0) {
        delete [] m_nz_row_ptr;
    }

    if (m_col != 0) {
        delete [] m_col;
    }

    if (m_val != 0) {
        delete [] m_val;
    }
}

int PSpMatrixForBackProj::getRowNumOrig() const
{
    return m_row_num_orig;
}

int PSpMatrixForBackProj::getNZRowNum() const
{
    return m_nz_row_num;
}

int PSpMatrixForBackProj::getColNum() const
{
    return m_col_num;
}

long long PSpMatrixForBackProj::getNNZ() const
{
    return m_nnz;
}

const int* PSpMatrixForBackProj::getCol() const
{
    return m_col;
}

const float* PSpMatrixForBackProj::getVal() const
{
    return m_val;
}

const int* PSpMatrixForBackProj::getNZRowMapping() const
{
    return m_nz_row_mapping;
}

const int* PSpMatrixForBackProj::getNZRowPtr() const
{
    return m_nz_row_ptr;
}

