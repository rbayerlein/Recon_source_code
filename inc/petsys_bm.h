#ifndef PETSYS_BM_H
#define PETSYS_BM_H

/*
	PatchOperator
*/
class BlockMatching {
public:
	BlockMatching(int block_half_size_ij, 
				  int block_half_size_k,
				  int searching_window_half_size_ij, 
				  int searching_window_half_size_k,
				  int sampling_step_ij, 
				  int sampling_step_k) : 
				  	m_block_half_size_ij(block_half_size_ij),
				  	m_block_half_size_k(block_half_size_k),
				  	m_searching_window_half_size_ij(searching_window_half_size_ij),
				  	m_searching_window_half_size_k(searching_window_half_size_k),
				  	m_sampling_step_ij(sampling_step_ij),
				  	m_sampling_step_k(sampling_step_k),
				  	m_block_sampling_index(0),
				  	m_block_sampling_weight(0) {

	};

	~BlockMatching() {
		if (m_block_sampling_weight != NULL) {
			delete m_block_sampling_weight;
		}

		if (m_block_sampling_index != NULL) {
			delete m_block_sampling_index;
		}
	};

public:

	void computeKernels(const Image<float>& image, const float h) {
		int ni = image.getDimJ();
		int nj = image.getDimJ();
		int nk = image.getDimK();
		int np = ni * nj * nk;

		int number_of_blocks_inside_searching_window_max = 
			(2*m_searching_window_half_size_ij+1)*(2*m_searching_window_half_size_ij+1)*(2*m_searching_window_half_size_k+1);

		Image<int> neighbor_block_index(number_of_blocks_inside_searching_window_max, np);
		Image<float> neighbor_block_weight(number_of_blocks_inside_searching_window_max, np);
		neighbor_block_index.set(-1);

		printf("it might take time ...\n");
		int n = 0;
#if USE_OMP
		#pragma omp parallel for private(n)
#endif
		for (n = 0; n < np; n ++) {
			computeKernelsForSingleBlock(image, n, h, 
										 neighbor_block_index.getPtr() + n * number_of_blocks_inside_searching_window_max,
										 neighbor_block_weight.getPtr() + n * number_of_blocks_inside_searching_window_max);
		}
		printf("OK!\n");

		// dump to disk
		neighbor_block_weight.write("ker_weight");
		neighbor_block_index.write("ker_ind");
	}

	Image<int>* computeKNearestNeighbors(const Image<float>& image, const int K) {
		int ni = image.getDimJ();
		int nj = image.getDimJ();
		int nk = image.getDimK();

		// sampling weight
		if (m_block_sampling_weight != NULL) {
			delete m_block_sampling_weight;
		}
		m_block_sampling_weight = new Image<float>(ni, nj, nk);

		// block sampling positions
		std::vector<int> index;
		for (int k = 0; k < nk; k += m_sampling_step_k) {
			for (int j = 0; j < nj; j += m_sampling_step_ij) {
				for (int i = 0; i < ni; i += m_sampling_step_ij) {
					index.push_back(image.sub2ind(i, j, k));
				}
			}
		}
		if (m_block_sampling_index != NULL) {
			delete m_block_sampling_index;
		}
		m_block_sampling_index = new Image<int>(index.size());
		memcpy(m_block_sampling_index->getPtr(), &index[0], sizeof(int) * index.size());

		int np = m_block_sampling_index->getSize();
#if DEBUG		
		printf("number of blocks in total: %d\n", np);
#endif

		Image<int>* knn = new Image<int>(K, np);

		int n = 0;
#if USE_OMP
		#pragma omp parallel for private(n)
#endif		
		for (n = 0; n < np; n ++) {
			computeKNearestNeighborsForSingleBlock(image, (*m_block_sampling_index)[n], K, knn->getPtr() + n * K);
		}

		return knn;
	};

	Image<float>* getBlockSamplingWeight() {
		return m_block_sampling_weight;
	}

	Image<int>* getBlockSamplingIndex() {
		return m_block_sampling_index;
	}

private:

	// for conventional neighborhood-based filters
	void computeKernelsForSingleBlock(const Image<float>& image,
									  const int current_block_index,
									  const float h, 
									  int* neighbor_block_index,
									  float* neighbor_block_weight) {
		// note: index is global voxel index over the image

		int ni = image.getDimI();
		int nj = image.getDimJ();
		int nk = image.getDimK();

		size_t i, j, k;
		image.ind2sub(current_block_index, i, j, k);

		// compute range of searching window
		int swi_min = (int(i) - m_searching_window_half_size_ij < 0) ? (0) : (i - m_searching_window_half_size_ij);
		int swi_max = (int(i) + m_searching_window_half_size_ij >= ni) ? (ni-1) : (i + m_searching_window_half_size_ij);
		int swj_min = (int(j) - m_searching_window_half_size_ij < 0) ? (0) : (j - m_searching_window_half_size_ij);
		int swj_max = (int(j) + m_searching_window_half_size_ij >= nj) ? (nj-1) : (j + m_searching_window_half_size_ij);
		int swk_min = (int(k) - m_searching_window_half_size_k < 0) ? (0) : (k - m_searching_window_half_size_k);
		int swk_max = (int(k) + m_searching_window_half_size_k >= nk) ? (nk-1) : (k + m_searching_window_half_size_k);

		// check patches one by one inside searching window
		int c = 0;
		for (int sk = swk_min; sk <= swk_max; sk ++) {
			for (int sj = swj_min; sj <= swj_max; sj ++) {
				for (int si = swi_min; si <= swi_max; si ++) {

					float d = 0.0;
					for (int kk = -m_block_half_size_k; kk <= m_block_half_size_k; kk ++) {
						int k1 = sk + kk;
						int k0 = k + kk;

						for (int jj = -m_block_half_size_ij; jj <= m_block_half_size_ij; jj ++) {
							int j1 = sj + jj;
							int j0 = j + jj;

							for (int ii = -m_block_half_size_ij; ii <= m_block_half_size_ij; ii ++) {
								int i1 = si + ii;
								int i0 = i + ii;

								float v1 = ((k1>=nk) || (k1 < 0) || 
											(j1>=nj) || (j1 < 0) ||
											(i1>=ni) || (i1 < 0)) ? 0.0 : image(i1, j1, k1);

								float v0 = ((k0>=nk) || (k0 < 0) || 
											(j0>=nj) || (j0 < 0) ||
											(i0>=ni) || (i0 < 0)) ? 0.0 : image(i0, j0, k0);

								d += (v0 - v1) * (v0 - v1); // l2-norm distance, simple!
							}
						}
					}

					neighbor_block_index[c] = image.sub2ind(si, sj, sk);
					neighbor_block_weight[c] = expf(-d / h);
					c ++;
					
				}
			}
		}

	}

	// for KNN-based Filter
	static bool comparison(std::pair<float, int> p1, std::pair<float, int> p2) { 
		return (p1.first < p2.first); 
	};

	void computeKNearestNeighborsForSingleBlock(const Image<float>& image,
												const int current_block_index,
												const int K,
												int* neighbor_block_index) {

		// note: index is global voxel index over the image

		int ni = image.getDimI();
		int nj = image.getDimJ();
		int nk = image.getDimK();

		size_t i, j, k;
		image.ind2sub(current_block_index, i, j, k);

		for (int kk = -m_block_half_size_k; kk <= m_block_half_size_k; kk ++) {
			int k0 = k + kk;
			for (int jj = -m_block_half_size_ij; jj <= m_block_half_size_ij; jj ++) {
				int j0 = j + jj;
				for (int ii = -m_block_half_size_ij; ii <= m_block_half_size_ij; ii ++) {
					int i0 = i + ii;
					if ((i0 >= 0) && (i0 < ni) &&
						(j0 >= 0) && (j0 < nj) &&
						(k0 >= 0) && (k0 < nk)) {  
						(*m_block_sampling_weight)(i0, j0, k0) += 1.0;
					}
				}
			}
		}

		// compute range of searching window
		int swi_min = (int(i) - m_searching_window_half_size_ij < 0) ? (0) : (i - m_searching_window_half_size_ij);
		int swi_max = (int(i) + m_searching_window_half_size_ij >= ni) ? (ni-1) : (i + m_searching_window_half_size_ij);
		int swj_min = (int(j) - m_searching_window_half_size_ij < 0) ? (0) : (j - m_searching_window_half_size_ij);
		int swj_max = (int(j) + m_searching_window_half_size_ij >= nj) ? (nj-1) : (j + m_searching_window_half_size_ij);
		int swk_min = (int(k) - m_searching_window_half_size_k < 0) ? (0) : (k - m_searching_window_half_size_k);
		int swk_max = (int(k) + m_searching_window_half_size_k >= nk) ? (nk-1) : (k + m_searching_window_half_size_k);

		// total number of voxels (as well as patches) inside searching window
//		int np_in_sw = (swi_max - swi_min + 1) * (swj_max - swj_min + 1) * (swk_max - swk_min + 1);
//		std::vector<std::pair<float, int> > diff_list(np_in_sw - 1);
		std::vector<std::pair<float, int> > diff_list;

		// check patches one by one inside searching window
		int c = 0;
		for (int sk = swk_min; sk <= swk_max; sk ++) {
			for (int sj = swj_min; sj <= swj_max; sj ++) {
				for (int si = swi_min; si <= swi_max; si ++) {

#if 0 // take all

#if 1
					// don't count the same patch again
					if ((si == i) && (sj == j) && (sk == k)) {
						continue;
					}
#else // don't count any block which has overlap with the target block
					if ((fabs(si - i) < (2*m_block_half_size_ij+1)) &&
						(fabs(sj - j) < (2*m_block_half_size_ij+1)) &&
						(fabs(sk - k) < (2*m_block_half_size_k+1)) ) {
						continue;
					}
#endif
#endif
					float d = 0.0;
					for (int kk = -m_block_half_size_k; kk <= m_block_half_size_k; kk ++) {
						int k1 = sk + kk;
						int k0 = k + kk;
						for (int jj = -m_block_half_size_ij; jj <= m_block_half_size_ij; jj ++) {
							int j1 = sj + jj;
							int j0 = j + jj;
							for (int ii = -m_block_half_size_ij; ii <= m_block_half_size_ij; ii ++) {
								int i1 = si + ii;
								int i0 = i + ii;

								float v1 = ((k1>=nk) || (k1 < 0) || 
											(j1>=nj) || (j1 < 0) ||
											(i1>=ni) || (i1 < 0)) ? 0.0 : image(i1, j1, k1);

								float v0 = ((k0>=nk) || (k0 < 0) || 
											(j0>=nj) || (j0 < 0) ||
											(i0>=ni) || (i0 < 0)) ? 0.0 : image(i0, j0, k0);

								d += (v0 - v1) * (v0 - v1); // l2-norm distance, simple!
							}
						}
					}

//					diff_list[c].first = d;
//					diff_list[c].second = image.sub2ind(si, sj, sk);
//					c ++;
					
					std::pair<float, int> bk;
					bk.first = d;
					bk.second = image.sub2ind(si, sj, sk);
					diff_list.push_back(bk);
				}
			}
		}

		// sort all patch differences, ascend
		std::sort(diff_list.begin(), diff_list.end(), comparison);

		// pick first K NN index
		for (int kk= 0; kk < K; kk ++) {
			neighbor_block_index[kk] = diff_list[kk].second;
		}

	};

public:

	void doBackwardFiltering(const Image<float>& image, 
						     const Image<int>& knn, 
						     const Image<int>& block_sampling_index,
						     const Image<float>& block_sampling_weight,						     
							 Image<float>& image_filtered) {

		int ni = image.getDimI();
		int nj = image.getDimJ();
		int nk = image.getDimK();
		int K = knn.getDimI();
		int np = knn.getDimJ();
		int n = 0;

		Image<float> image_temp(ni,nj,nk);
		memcpy(image_temp.getPtr(), image.getPtr(), sizeof(float)*image.getSize());

		// applying sampling weight first
		for (int n = 0; n < image.getSize(); n ++) {
			image_temp[n] /= block_sampling_weight[n];
		}

#if USE_OMP
		#pragma omp parallel for private(n)
#endif
		for (n = 0; n < np; n ++) {

			size_t i, j, k;
			image.ind2sub(block_sampling_index[n], i, j, k);

			for (int l = 0; l < K; l ++) {

				size_t ii, jj, kk;
				image.ind2sub(knn(l, n), ii, jj, kk);

				for (int pk=-m_block_half_size_k; pk<=m_block_half_size_k; pk ++) {
					int kkk = kk + pk;
					int k0 = k + pk;
					for (int pj=-m_block_half_size_ij; pj<=m_block_half_size_ij; pj ++) {
						int jjj = jj + pj;
						int j0 = j + pj;
						for (int pi=-m_block_half_size_ij; pi<=m_block_half_size_ij; pi ++) {
							int iii = ii + pi;
							int i0 = i + pi;

							float v = ((k0>=nk) || (k0 < 0) || 
									   (j0>=nj) || (j0 < 0) ||
									   (i0>=ni) || (i0 < 0))  ? 0.0 : image_temp(i0,j0,k0);

							if 	(!((kkk>=nk) || (kkk < 0) || 
								 (jjj>=nj) || (jjj < 0) ||
								 (iii>=ni) || (iii < 0))) {
#if USE_OMP
								#pragma omp atomic
#endif								
								image_filtered(iii,jjj,kkk) += v / K;
							} 

						}
					}		
				}

			}

		}

	};

	void doForwardFiltering(const Image<float>& image, 
						    const Image<int>& knn, 
						    const Image<int>& block_sampling_index,
						    const Image<float>& block_sampling_weight,
							Image<float>& image_filtered) {

		int ni = image.getDimI();
		int nj = image.getDimJ();
		int nk = image.getDimK();
		int K = knn.getDimI();
		int np = knn.getDimJ();
		int n = 0;

#if USE_OMP
		#pragma omp parallel for private(n)
#endif
		for (n = 0; n < np; n ++) {

			size_t i, j, k;
			image.ind2sub(block_sampling_index[n], i, j, k);

			for (int l = 0; l < K; l ++) {

				size_t ii, jj, kk;
				image.ind2sub(knn(l, n), ii, jj, kk);

				for (int pk=-m_block_half_size_k; pk<=m_block_half_size_k; pk ++) {
					int kkk = kk + pk;
					int k0 = k + pk;
					for (int pj=-m_block_half_size_ij; pj<=m_block_half_size_ij; pj ++) {
						int jjj = jj + pj;
						int j0 = j + pj;
						for (int pi=-m_block_half_size_ij; pi<=m_block_half_size_ij; pi ++) {
							int iii = ii + pi;
							int i0 = i + pi;

							float v = ((kkk>=nk) || (kkk < 0) || 
									   (jjj>=nj) || (jjj < 0) ||
									   (iii>=ni) || (iii < 0))  ? 0.0 : image(iii,jjj,kkk);

							if 	(!((k0>=nk) || (k0 < 0) || 
								 (j0>=nj) || (j0 < 0) ||
								 (i0>=ni) || (i0 < 0))) {
#if USE_OMP
								#pragma omp atomic
#endif								
								image_filtered(i0,j0,k0) += v / float(K);
							} 

						}
					}		
				}
			}

		}

		// applying sampling weight
		for (int n = 0; n < image_filtered.getSize(); n ++) {
			image_filtered[n] /= block_sampling_weight[n];
		}

	}

private:
	int m_block_half_size_ij;
	int m_block_half_size_k;
	int m_searching_window_half_size_ij;
	int m_searching_window_half_size_k; 
	int m_sampling_step_ij;
	int m_sampling_step_k;
	Image<float>* m_block_sampling_weight;
	Image<int>* m_block_sampling_index;
};

#endif