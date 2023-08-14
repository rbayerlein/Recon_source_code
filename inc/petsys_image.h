/*!
 * \file petsys_image.h
 *
 * \brief Define an Image class to manage images or matrices up to 3D dimension.
 *
 * \author Jian Zhou
 * \date 2011-10-05
 * \since 0.1
 *
 * Copyright (c) 2011, Jian Zhou. All rights reserved.
 */
#ifndef PETSYS_IMAGE_H
#define PETSYS_IMAGE_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <petsys_sysparms.h>

namespace UCDPETSYS
{

/* not implemented yet
	void copy(const T* source) { memcpy(m_data, source, sizeof(T)*getVoxelNum()); }
	void copy(const T* source, size_t size) { memcpy(m_data, source, sizeof(T)*size); }
*/
template <typename T>
class Image
{
public:
    Image(const size_t dim_i = 0, const size_t dim_j = 1, const size_t dim_k = 1);
    Image(const SIZE& sz);
    Image(const size_t sz[]);
    Image(const Image<T>& img);
    ~Image();

public:
    size_t getDimI() const;
    size_t getDimJ() const;
    size_t getDimK() const;
    size_t getSize() const;
    T* getPtr() const; ///< will be deprecated
    void set(const T& val);
    void set(const size_t i, const size_t j, const size_t k, const T& val);
    size_t sub2ind(const size_t i, const size_t j, const size_t k) const;
    void ind2sub(const size_t idx, size_t& i, size_t& j, size_t& k) const;
    bool read(const char* filename);
    bool write(const char* filename) const;
    void clear();

    Image<T>& operator =(const Image<T>& img);
    T& operator [](const size_t index);
    const T& operator [](const size_t index) const;
    T& operator ()(const size_t i = 0, const size_t j = 0, const size_t k = 0);
    const T& operator ()(const size_t i = 0, const size_t j = 0, const size_t k = 0) const;

private:
    void alloc(size_t sz);
    void release();

private:
    T* 	m_data;
    size_t	m_dim_i;
    size_t	m_dim_j;
    size_t	m_dim_k;
    size_t m_dim_ij; ///< = m_dim_i*m_dim*j
};

}

#include "../src/details/petsys_image.impl"

#endif // petsys_image.h

