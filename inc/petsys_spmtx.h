/*!
 * \file petsys_spmtx.h
 *
 * \brief Define two sparse matrix classes (SSpMatrix (with data quantization) and PSpMatrix)
 *
 * \author Jian Zhou
 * \date 2011-10-05
 * \since 0.1
 *
 * Copyright (c) 2011, Jian Zhou. All rights reserved.
 */
#ifndef PETSYS_SPMTX_H
#define PETSYS_SPMTX_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <petsys_image.h>

namespace UCDPETSYS
{

class PSpMatrix;

class SSpMatrix
{

public:
    typedef struct __attribute__ ((__packed__)) tagPACKAGE {
        unsigned i : 11; ///< image i index
        unsigned j : 11; ///< image j index
#if USE_LONG_AXIAL_FOV
        unsigned k : 17; ///< image k index
#else
        unsigned k : 9;  ///< image k index
#endif

#if USE_9BIT
        unsigned w : 9; ///< quantized to 9 bits
#else
        unsigned w : 17;  ///< quantized to 17 bits
#endif
    } PACKAGE;

public:
    SSpMatrix(const int row_num = 0, const int col_num = 0, const long long nnz = 0,
              const float scaling_factor = 0, const int* row_ptr = 0, const PACKAGE* p = 0);

    explicit SSpMatrix(const PSpMatrix& pmat, const SIZE& imgsize);

    ~SSpMatrix();


public:
    void read(const char* smat_filename);
    void write(const char* smat_filename);
    int getRowNum() const;
    int getColNum() const;
    long long getNNZ() const;
    float getScalingFactor() const;
    const int* getRowPtr() const; // offset table, name from old version
    const PACKAGE* getP() const;
    SSpMatrix* pickRows(const std::vector<int>& row_idx);

public:
    enum PROJ_type {
        FORWARD = 0,
        BACKWARD
    };
    void proj(Image<float>& y, Image<float>& x, PROJ_type t = FORWARD);

private:
    int m_row_num;
    int m_col_num;
    long long m_nnz;
    float m_scaling_factor;
    int* m_row_ptr;
    PACKAGE* m_p;
};

class PSpMatrix
{
public:
    PSpMatrix(const int row_num = 0, const int col_num = 0, const long long nnz = 0,
              const int* row_ptr = 0, const int* col = 0, const float* val = 0);
    ~PSpMatrix();

public:
    bool read(const char* pmat_filename);
    bool write(const char* pmat_filename);
    void proj(float* y, float* x, SSpMatrix::PROJ_type t = SSpMatrix::FORWARD);
    int getRowNum() const;
    int getColNum() const;
    long long getNNZ() const;
    const int* getRowPtr() const;
    const int* getCol() const;
    const float* getVal() const;
    PSpMatrix* getTranspose();
    PSpMatrix* pickRows(const std::vector<int>& row_idx);

private:
    typedef struct tagIJK {
        int i;
        int j;
        float k;
    } IJK;

private:
    int m_row_num;
    int m_col_num;
    long long m_nnz;
    int* m_row_ptr;
    int* m_col;
    float* m_val;
};

class PSpMatrixForBackProj
{
public:
    PSpMatrixForBackProj(const PSpMatrix& pmat);
    ~PSpMatrixForBackProj();

public:
    int getRowNumOrig() const;
    int getNZRowNum() const;
    int getColNum() const;
    long long getNNZ() const;
    const int* getCol() const;
    const float* getVal() const;
    const int* getNZRowMapping() const;
    const int* getNZRowPtr() const;

private:
    int m_row_num_orig;
    int m_nz_row_num;
    int m_col_num;
    long long m_nnz;
    int* m_nz_row_mapping;
    int* m_nz_row_ptr;
    int* m_col;
    float* m_val;
};

};

#include "../src/details/petsys_spmtx.impl"

#endif // petsys_spmtx.h
