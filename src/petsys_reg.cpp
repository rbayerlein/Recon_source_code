#include <petsys_reg.h>

// int Regularizer::NUMBER_OF_THREADS = 32;
int Regularizer::NUMBER_OF_THREADS = 64;


Regularizer::Regularizer() :
m_potential_function_type(0),
m_neighborhood_size(0),
m_is_isotropic(false),
m_buildin_parms(0),
m_spatial_variant_weight(0)
{
}

Regularizer::~Regularizer()
{
	if (m_buildin_parms != 0) {
		delete m_buildin_parms;
	}
	
	if (m_spatial_variant_weight != 0) {
		delete m_spatial_variant_weight;
	}
}

void Regularizer::initialize(int potential_function_type,
    		   				 int neighborhood_size,
    		   				 bool is_isotropic,
    		   				 double* buildin_parms,
    		   				 int number_of_buildin_parms_detected)
{
	m_potential_function_type = potential_function_type;
	m_neighborhood_size = neighborhood_size;
	m_is_isotropic = is_isotropic;
	if (number_of_buildin_parms_detected > 0) {
		m_buildin_parms = new double[number_of_buildin_parms_detected];
		memcpy(m_buildin_parms, buildin_parms, 
			   sizeof(double) * number_of_buildin_parms_detected);
	}

	switch (m_potential_function_type) {
		case QUAD:
			SystemLog::write("using quadratic penalty ...\n");
			if (number_of_buildin_parms_detected > 0) {
				SystemLog::write("warning, qudratic penalty does not require buildin parameters. "
								 "ignored here.\n");
			}
			break;
		case HYPERBOLA:
			SystemLog::write("using hyperbola penalty ...\n");
			if (number_of_buildin_parms_detected == 0) {
				SystemLog::write("error, hyperbola requires 1 buildin parameter.\n");
				abort();
			} else {
				if (number_of_buildin_parms_detected > 1) {
					SystemLog::write("warning, too many buildin parameters for hyperbola, "
									 "will take the first parameter only.\n");
					if (m_buildin_parms[0] < 0.0) {
						SystemLog::write("invalid buildin parameter which should be positive.\n");
						abort();						
					}
				} 
				SystemLog::write("delta = %.2e\n", m_buildin_parms[0]);
			}
			break;
		case FAIR:
			SystemLog::write("using fair penalty ...\n");
			if (number_of_buildin_parms_detected == 0) {
				SystemLog::write("error, fair requires 1 buildin parameter.\n");
				abort();
			} else {
				if (number_of_buildin_parms_detected > 1) {
					SystemLog::write("warning, too many buildin parameters for hyperbola, "
									 "will take the first parameter only.\n");
					if (m_buildin_parms[0] < 0.0) {
						SystemLog::write("invalid buildin parameter which should be positive.\n");
						abort();						
					}
				} 
				SystemLog::write("delta = %.2e\n", m_buildin_parms[0]);
			}
			break;
		case HL:
			SystemLog::write("using hl penalty ...\n");
			if (number_of_buildin_parms_detected == 0) {
				SystemLog::write("error, hl requires 1 buildin parameter.\n");
				abort();
			} else {
				if (number_of_buildin_parms_detected > 1) {
					SystemLog::write("warning, too many buildin parameters for hyperbola, "
									 "will take the first parameter only.\n");
					if (m_buildin_parms[0] < 0.0) {
						SystemLog::write("invalid buildin parameter which should be positive.\n");
						abort();						
					}
				} 
				SystemLog::write("delta = %.2e\n", m_buildin_parms[0]);
			}
			break;
		case HW:
			SystemLog::write("using hw penalty ...\n");
			if (number_of_buildin_parms_detected == 0) {
				SystemLog::write("error, hw requires 1 buildin parameter.\n");
				abort();
			} else {
				if (number_of_buildin_parms_detected > 1) {
					SystemLog::write("warning, too many buildin parameters for hyperbola, "
									 "will take the first parameter only.\n");
					if (m_buildin_parms[0] < 0.0) {
						SystemLog::write("invalid buildin parameter which should be positive.\n");
						abort();						
					}
				} 
				SystemLog::write("delta = %.2e\n", m_buildin_parms[0]);
			}
			break;
		case HU:
			SystemLog::write("using Huber's penalty ...\n");
			if (number_of_buildin_parms_detected == 0) {
				SystemLog::write("error, Huber needs one buildon parameter.\n");
				abort();
			} else {
				if (number_of_buildin_parms_detected > 1) {
					SystemLog::write("warning, too many buildin parameters in list,"
									 "will take the first parameter only.\n");
					if (m_buildin_parms[0] < 0.0) {
						SystemLog::write("invalid buildin parameter which must be positive.\n");
						abort();
					}
				}
				SystemLog::write("delta = %.2e\n", m_buildin_parms[0]);
			}
			break;
		default:
			SystemLog::write("unknown potential function type [%d], not supported yet.\n",
				m_potential_function_type);
			abort();
			break;
	}

}

inline double Regularizer::phi(double t)
{
	double delta, h;
	switch (m_potential_function_type) {
		case QUAD:
			return (t * t);
		case HYPERBOLA:
			return (sqrt(t * t + m_buildin_parms[0]));
		case FAIR:
			delta = m_buildin_parms[0];
			h = t / delta;
			return (delta * (h - log(1.0 + fabs(h))));
		case HL:
			delta = m_buildin_parms[0];
			h = t / delta;
			return (delta * delta * 0.5 * log(1.0 + h * h));
		case HW:
			h = t / m_buildin_parms[0];
			return (1.0 - exp(- h * h));
		case HU:
			delta = m_buildin_parms[0];
			h = (fabs(t) <= delta) ? t*t*0.5 : delta*fabs(t)-delta*delta/2;
			return h;
	}
}

inline double Regularizer::dot_phi_over_t(double t)
{
	double h, delta;
	switch (m_potential_function_type) {
		case QUAD:
			return 2.0;
		case HYPERBOLA:
			return (1.0 / sqrt(t * t + m_buildin_parms[0]));
		case FAIR:
			return (1.0 / (fabs(t) + m_buildin_parms[0]));
		case HL:
			h = t / m_buildin_parms[0];
			return (1.0 / (1.0 + h * h));
		case HW:
			h = t / m_buildin_parms[0];
			return exp(-(h * h));
		case HU:
			delta = m_buildin_parms[0];
			h = (fabs(t)<=delta) ? 1.0 : delta / fabs(t);
			return h;
	}
}

inline double Regularizer::dot_phi(double t)
{
	double h, delta;
	switch (m_potential_function_type) {
		case QUAD:
			return 2.0 * t;
		case HYPERBOLA:
			return (t / sqrt(t * t + m_buildin_parms[0]));
		case FAIR:
			return (t / (fabs(t) + m_buildin_parms[0]));
		case HL:
			h = t / m_buildin_parms[0];
			return (t / (1.0 + h * h));
		case HW:
			h = t / m_buildin_parms[0];
			return (t * exp(-(h * h)));
		case HU:
			delta = m_buildin_parms[0];
			h = (fabs(t)<=delta) ? t : (t>0 ? delta : (t==0 ? 0 : -delta));
			return h;
	}
}

void Regularizer::setSpatialVariantWeight(const Image<float>& weight)
{
	if (m_spatial_variant_weight != 0) {
		delete m_spatial_variant_weight;
	}

	m_spatial_variant_weight = new Image<float>(weight.getDimI(), 
												weight.getDimJ(), 
												weight.getDimK());
	memcpy(m_spatial_variant_weight->getPtr(),
		   weight.getPtr(), sizeof(float) * weight.getSize());
}

void Regularizer::calculateGradientBasedOnOptimizationTransfer(const Image<float>& current_image,
        Image<float>& term1,
        Image<float>& term2)
{
}

// void Regularizer::gbcalg(const Image<float>& x,
//         Image<float>& t1,
//         Image<float>& t2,
//         const Image<float>& cip,
//         Image<int>* fd)
// {
// }

void Regularizer::calculateGradient(const Image<float>& x, Image<float>& grad)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
PairwiseMRFRegularizer::PairwiseMRFRegularizer()
{
	// Regularizer::Regularizer();
	Regularizer();
}

PairwiseMRFRegularizer::~PairwiseMRFRegularizer()
{
	for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
		if (m_neighbor_voxel_position_offset[n] != 0) {
			delete m_neighbor_voxel_position_offset[n];
		}
	}

//	if (m_spatial_variant_weight != 0) {
//		delete m_spatial_variant_weight;
//	}
}

void PairwiseMRFRegularizer::initialize(int potential_function_type,
    		   				 			int neighborhood_size,
    		   				 			bool is_isotropic,
    		   				 			double* buildin_parms,
    		   				 			int number_of_buildin_parms_detected)
{
	SystemLog::write("initializing pairwise MRF regularizer ...\n");

	Regularizer::initialize(potential_function_type,
				 			neighborhood_size,
				 			is_isotropic,
				 			buildin_parms,
				 			number_of_buildin_parms_detected);
	
	switch (neighborhood_size) {
		case 0:
			SystemLog::write("1st order neighborhood selected.\n");
			break;
		case 1:
			SystemLog::write("2nd order neighborhood selected.\n");
			break;
		default:
			SystemLog::write("larger neighborhood size selected. it takes time!\n");
//							 "valid value is 0 (1st order) or 1 (2nd order).\n");
//			abort();
			break;
	}
	
	SystemLog::write("creating offsets ...\n");
	createNeighborVoxelPositionOffset();
	
	SystemLog::write("[ns:%d, iso:%d], ok.\n", m_neighborhood_size, m_is_isotropic);	 			
}

void PairwiseMRFRegularizer::createNeighborVoxelPositionOffset()
{	
	if (m_neighborhood_size == 0) {
	
		// modified by jnzhou@04-21-2014
		// double the number of filters (because of assymmetric filter)
		// old version does affect the isotropic form very much
		m_neighbor_voxel_position_offset.resize(6); 
		for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			m_neighbor_voxel_position_offset[n] = new Image<float>(3);
			m_neighbor_voxel_position_offset[n]->clear();
		}
		
		// first-order neighbor (direct neighbor)
		(*m_neighbor_voxel_position_offset[0])(0) = 1; // i+1, j, k
		(*m_neighbor_voxel_position_offset[1])(1) = 1; // i, j+1, k
		(*m_neighbor_voxel_position_offset[2])(2) = 1; // i, j, k+1	
		(*m_neighbor_voxel_position_offset[3])(0) = -1; // i-1, j, k
		(*m_neighbor_voxel_position_offset[4])(1) = -1; // i, j-1, k
		(*m_neighbor_voxel_position_offset[5])(2) = -1; // i, j, k-1	
		
	} else {
		
		if (m_neighborhood_size == 1) {
			
			// modified by jnzhou@04-21-2014
			m_neighbor_voxel_position_offset.resize(26); 
			for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
				m_neighbor_voxel_position_offset[n] = new Image<float>(3);
				m_neighbor_voxel_position_offset[n]->clear();
			}
		
			// first-order neighbor
			(*m_neighbor_voxel_position_offset[0])(0) = 1; // i+1, j, k
			(*m_neighbor_voxel_position_offset[1])(1) = 1; // i, j+1, k
			(*m_neighbor_voxel_position_offset[2])(2) = 1; // i, j, k+1

			// including diagonal neighbors
			(*m_neighbor_voxel_position_offset[3])(0) = 1;
			(*m_neighbor_voxel_position_offset[3])(1) = 1;
			(*m_neighbor_voxel_position_offset[4])(0) = 1;
			(*m_neighbor_voxel_position_offset[4])(1) = -1;
			//
			(*m_neighbor_voxel_position_offset[5])(1) = 1;
			(*m_neighbor_voxel_position_offset[5])(2) = 1;
			(*m_neighbor_voxel_position_offset[6])(1) = 1;
			(*m_neighbor_voxel_position_offset[6])(2) = -1;
			//
			(*m_neighbor_voxel_position_offset[7])(2) = 1;
			(*m_neighbor_voxel_position_offset[7])(0) = 1;
			(*m_neighbor_voxel_position_offset[8])(2) = 1;
			(*m_neighbor_voxel_position_offset[8])(0) = -1;
		
			// more 
			(*m_neighbor_voxel_position_offset[9])(0) = 1;
			(*m_neighbor_voxel_position_offset[9])(1) = 1;
			(*m_neighbor_voxel_position_offset[9])(2) = 1;
			//
			(*m_neighbor_voxel_position_offset[10])(0) = 1;
			(*m_neighbor_voxel_position_offset[10])(1) = 1;
			(*m_neighbor_voxel_position_offset[10])(2) = -1;
			//
			(*m_neighbor_voxel_position_offset[11])(0) = 1;
			(*m_neighbor_voxel_position_offset[11])(1) = -1;
			(*m_neighbor_voxel_position_offset[11])(2) = 1;
			//
			(*m_neighbor_voxel_position_offset[12])(0) = 1;
			(*m_neighbor_voxel_position_offset[12])(1) = -1;
			(*m_neighbor_voxel_position_offset[12])(2) = -1;

			// symmetric	
			// first-order neighbor
			(*m_neighbor_voxel_position_offset[0 + 13])(0) = -1; // i+1, j, k
			(*m_neighbor_voxel_position_offset[1 + 13])(1) = -1; // i, j+1, k
			(*m_neighbor_voxel_position_offset[2 + 13])(2) = -1; // i, j, k+1

			// including diagonal neighbors
			(*m_neighbor_voxel_position_offset[3 + 13])(0) = -1;
			(*m_neighbor_voxel_position_offset[3 + 13])(1) = -1;
			(*m_neighbor_voxel_position_offset[4 + 13])(0) = -1;
			(*m_neighbor_voxel_position_offset[4 + 13])(1) = 1;
			//
			(*m_neighbor_voxel_position_offset[5 + 13])(1) = -1;
			(*m_neighbor_voxel_position_offset[5 + 13])(2) = -1;
			(*m_neighbor_voxel_position_offset[6 + 13])(1) = -1;
			(*m_neighbor_voxel_position_offset[6 + 13])(2) = 1;
			//
			(*m_neighbor_voxel_position_offset[7 + 13])(2) = -1;
			(*m_neighbor_voxel_position_offset[7 + 13])(0) = -1;
			(*m_neighbor_voxel_position_offset[8 + 13])(2) = -1;
			(*m_neighbor_voxel_position_offset[8 + 13])(0) = 1;
		
			// more 
			(*m_neighbor_voxel_position_offset[9 + 13])(0) = -1;
			(*m_neighbor_voxel_position_offset[9 + 13])(1) = -1;
			(*m_neighbor_voxel_position_offset[9 + 13])(2) = -1;
			//
			(*m_neighbor_voxel_position_offset[10 + 13])(0) = -1;
			(*m_neighbor_voxel_position_offset[10 + 13])(1) = -1;
			(*m_neighbor_voxel_position_offset[10 + 13])(2) = 1;
			//
			(*m_neighbor_voxel_position_offset[11 + 13])(0) = -1;
			(*m_neighbor_voxel_position_offset[11 + 13])(1) = 1;
			(*m_neighbor_voxel_position_offset[11 + 13])(2) = -1;
			//
			(*m_neighbor_voxel_position_offset[12 + 13])(0) = -1;
			(*m_neighbor_voxel_position_offset[12 + 13])(1) = 1;
			(*m_neighbor_voxel_position_offset[12 + 13])(2) = 1;
		
		} else {

			if (m_neighborhood_size >= 2) {

				// added by jnzhou@06-12-2014
				int hw = m_neighborhood_size;
				int ps = 2*m_neighborhood_size + 1;
				m_neighbor_voxel_position_offset.resize(ps * ps * ps - 1); 
				for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
					m_neighbor_voxel_position_offset[n] = new Image<float>(3);
					m_neighbor_voxel_position_offset[n]->clear();
				}

				int c = 0;
				for (int k = -hw; k <= hw; k ++) {
					for (int j = -hw; j <= hw; j ++) {
						for (int i = -hw; i <= hw; i ++) {
							if ((i == 0) && (j == 0) && (k == 0)) {
								continue;
							}
							(*m_neighbor_voxel_position_offset[c])(0) = i;
							(*m_neighbor_voxel_position_offset[c])(1) = j;							
							(*m_neighbor_voxel_position_offset[c])(2) = k;
							c ++;
						}	
					}					
				}

			} else {
				SystemLog::write("neighbor size not supported yet. "
							 "valid sizes are 0, 1 (3 x 3 x 3), 2 (5 x 5 x 5), ...\n");
				abort();
			}	

		}
	}
	
	m_geometrical_weight.resize(m_neighbor_voxel_position_offset.size());
	for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
		m_geometrical_weight[n] = 1.0 / sqrt((*m_neighbor_voxel_position_offset[n])(0) * 
							   		(*m_neighbor_voxel_position_offset[n])(0) +
							   		(*m_neighbor_voxel_position_offset[n])(1) * 
							   		(*m_neighbor_voxel_position_offset[n])(1) + 
							   		(*m_neighbor_voxel_position_offset[n])(2) * 
							   		(*m_neighbor_voxel_position_offset[n])(2)); 
	}
}

void PairwiseMRFRegularizer::calculateGradientAnisotropic(const Image<float>& x,
									  					  const Image<float>& neighbor_voxel_position_offset,
									  					  const float geometrical_weight,
									  					  Image<float>& grad)									  
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	//
	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// voxel j
				float voxel_j = x(i, j, k);
			
				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 :
						   	   (*m_spatial_variant_weight)(i, j, k);

				// voxel k, neighbor of voxel j
				int i0 = i + neighbor_voxel_position_offset(0);
				int j0 = j + neighbor_voxel_position_offset(1);
				int k0 = k + neighbor_voxel_position_offset(2);
				float voxel_k = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : x(i0, j0, k0);
				float voxel_k1 = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : 1.0;

/*				
				float kappa2 = ((m_spatial_variant_weight == 0) || (i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 1.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0);
*/
/*				float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
								( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0));
*/

				float kappa2 = ( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(i0, j0, k0) );					

				// for weight;
				float diff = voxel_j - voxel_k;
				float weight = dot_phi(diff) * geometrical_weight * kappa1 * kappa2;
				
				// insert back
				grad(i, j, k) += weight;
				if (!((i0<0) || (i0>=ni) || 
					  (j0<0) || (j0>=nj) || 
					  (k0<0) || (k0>=nk))) {
					grad(i0, j0, k0) -= weight;
				}
			}
		}
	}	
}			

void PairwiseMRFRegularizer::calculateGradientIsotropic(const Image<float>& x,
														const Image<float>& neighbor_voxel_position_offset,
								    					const Image<float>& image_weight,
														const float geometrical_weight,
														Image<float>& grad) 
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	//
	for (int k = 0; k < nk; k ++) {

		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// voxel j
				float voxel_j = x(i, j, k);

				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 :
						   	   (*m_spatial_variant_weight)(i, j, k);
			
				// voxel k, neighbor of voxel j
				int i0 = i + neighbor_voxel_position_offset(0);
				int j0 = j + neighbor_voxel_position_offset(1);
				int k0 = k + neighbor_voxel_position_offset(2);
				float voxel_k = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : x(i0, j0, k0);
				float voxel_k1 = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : 1.0;

/*
				float kappa2 = ((m_spatial_variant_weight == 0) || (i0<0) || (i0>=ni) || 
								(j0<0) || (j0>=nj) || 
								(k0<0) || (k0>=nk)) ? 1.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0);
*/
/*				float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
								( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0));				
*/
				float kappa2 = ( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(i0, j0, k0) );					

				
				// for weight;
				float diff = voxel_j - voxel_k;
				float weight = image_weight(i, j, k) * geometrical_weight * kappa1 * kappa2;
				
				// insert back
				grad(i, j, k) += weight * diff;

				if (!((i0<0) || (i0>=ni) || 
					  (j0<0) || (j0>=nj) || 
					  (k0<0) || (k0>=nk))) {
					grad(i0, j0, k0) -= weight * diff;
				}
			}
		}
	}	
}														

void PairwiseMRFRegularizer::calculateGradient(const Image<float>& x, Image<float>& grad)
{
#if USE_OMP	 
	 
	std::vector<Image<float>* > grad_list;
	grad_list.resize(m_neighbor_voxel_position_offset.size());

	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	if (!m_is_isotropic) { // anisotropic
	 
		size_t n = 0;
	 	#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			grad_list[n] = new Image<float>(ni, nj, nk);
			calculateGradientAnisotropic(x, *m_neighbor_voxel_position_offset[n], 
										 m_geometrical_weight[n], 
										 *grad_list[n]);
		}
			
	} else { // isotropic
	
		Image<float>* weight = calculateWeightIsotropic(x);

	 	size_t n = 0;
	 	#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			grad_list[n] = new Image<float>(ni, nj, nk);
			calculateGradientIsotropic(x, *m_neighbor_voxel_position_offset[n], 
									   *weight, m_geometrical_weight[n], 
									   *grad_list[n]);
		}
		
		delete weight;
		
	}
	
	for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
		int j = 0;
		#pragma omp parallel for private(j) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (j = 0; j < x.getSize(); j ++) {
			grad[j] += (*grad_list[n])[j];
		}
			
		delete grad_list[n];
	}
	
#else

	// anisotropic case only
	if (!m_is_isotropic) {
		for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			calculateGradientAnisotropic(x, *m_neighbor_voxel_position_offset[n], 
										 m_geometrical_weight[n], grad);
		}	
	} else { // isotropic
		
		Image<float>* weight = calculateWeightIsotropic(x);
				
		for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			calculateGradientIsotropic(x, *m_neighbor_voxel_position_offset[n], 
									   *weight, m_geometrical_weight[n], grad);
		}
		
		delete weight;
	} 
#endif	
}

void PairwiseMRFRegularizer::calculateTerm1AndTerm2Anisotropic(const Image<float>& x,
										   					   const Image<float>& neighbor_voxel_position_offset,
										   					   const float geometrical_weight,
										   					   Image<float>& term1,
										   					   Image<float>& term2)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	//
	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// voxel j
				float voxel_j = x(i, j, k);
			
				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 :
						   	   (*m_spatial_variant_weight)(i, j, k);

				// voxel k, neighbor of voxel j
				int i0 = i + neighbor_voxel_position_offset(0);
				int j0 = j + neighbor_voxel_position_offset(1);
				int k0 = k + neighbor_voxel_position_offset(2);
				float voxel_k = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : x(i0, j0, k0);
				float voxel_k1 = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : 1.0;

/*				
				float kappa2 = ((m_spatial_variant_weight == 0) || (i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 1.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0);
*/
/*								
				float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
								( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0));				
*/

				float kappa2 = ( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(i0, j0, k0) );					


				// for weight;
				float diff = voxel_j - voxel_k;
				float sum = 1.0 + voxel_k1;
				float weight = dot_phi_over_t(voxel_j - voxel_k) * 
								geometrical_weight * kappa1 * kappa2;
				
				// insert back
				term2(i, j, k) += diff * weight;
				term1(i, j, k) += sum * weight;		
				if (!((i0<0) || (i0>=ni) || 
					  (j0<0) || (j0>=nj) || 
					  (k0<0) || (k0>=nk))) {
					term2(i0, j0, k0) -= diff * weight;
					term1(i0, j0, k0) += sum * weight;
				}
			}
		}
	}
}

void PairwiseMRFRegularizer::calculateGradientBasedOnOptimizationTransfer(const Image<float>& x,
        																 Image<float>& t1,
        																 Image<float>& t2)
{ 
#if USE_OMP	 
	 
	std::vector<Image<float>* > t1_list;
	std::vector<Image<float>* > t2_list;
	t1_list.resize(m_neighbor_voxel_position_offset.size());
	t2_list.resize(m_neighbor_voxel_position_offset.size());

	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	if (!m_is_isotropic) {
	 
	 	size_t n = 0;
	 	#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			t1_list[n] = new Image<float>(ni, nj, nk);
			t2_list[n] = new Image<float>(ni, nj, nk);
			calculateTerm1AndTerm2Anisotropic(x, *m_neighbor_voxel_position_offset[n], 
											  m_geometrical_weight[n], 
											  *t1_list[n], *t2_list[n]);
		}
			
	} else { // isotropic
		Image<float>* weight = calculateWeightIsotropic(x);

	 	size_t n = 0;
	 	#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			t1_list[n] = new Image<float>(ni, nj, nk);
			t2_list[n] = new Image<float>(ni, nj, nk);
			calculateTerm1AndTerm2Isotropic(x, *m_neighbor_voxel_position_offset[n], 
											*weight, m_geometrical_weight[n], 
											*t1_list[n], *t2_list[n]);
		}
		
		delete weight;
	}
	
	for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
		int j = 0;
		#pragma omp parallel for private(j) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (j = 0; j < x.getSize(); j ++) {
			t1[j] += (*t1_list[n])[j];
			t2[j] += (*t2_list[n])[j];
		}
			
		delete t1_list[n];
		delete t2_list[n];
	}
	
#else
	// anisotropic case
	if (!m_is_isotropic) {
		for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			calculateTerm1AndTerm2Anisotropic(x, *m_neighbor_voxel_position_offset[n], 
											  m_geometrical_weight[n], t1, t2);
		}	
	} else { // isotropic
		
		Image<float>* weight = calculateWeightIsotropic(x);
				
		for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
			calculateTerm1AndTerm2Isotropic(x, *m_neighbor_voxel_position_offset[n], 
											*weight, m_geometrical_weight[n], t1, t2);
		}
		
		delete weight;
	}
#endif
	
}

Image<float>* PairwiseMRFRegularizer::calculateWeightIsotropic(const Image<float>& x)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	Image<float>* weight = new Image<float>(ni, nj, nk);
	
	//
	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// voxel j
				float voxel_j = x(i, j, k);
				
				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 : //?
						   	   (*m_spatial_variant_weight)(i, j, k);

				// voxel k, neighbor of voxel j
				float s = 0.0;
				for (size_t n = 0; n < m_neighbor_voxel_position_offset.size(); n ++) {
								
					int i0 = i + (*m_neighbor_voxel_position_offset[n])(0);
					int j0 = j + (*m_neighbor_voxel_position_offset[n])(1);
					int k0 = k + (*m_neighbor_voxel_position_offset[n])(2);
					
					float voxel_k = ((i0<0) || (i0>=ni) || 
									 (j0<0) || (j0>=nj) || 
									 (k0<0) || (k0>=nk)) ? 
									 0.0 : x(i0, j0, k0);
/*				
					float kappa2 = ((m_spatial_variant_weight == 0) || (i0<0) || (i0>=ni) || 
									 (j0<0) || (j0>=nj) ||  
									 (k0<0) || (k0>=nk)) ? 1.0 : // incorrect!
									(*m_spatial_variant_weight)(i0, j0, k0);
*/
/*									
					float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
									( ((i0<0) || (i0>=ni) || 
									 (j0<0) || (j0>=nj) ||  
									 (k0<0) || (k0>=nk)) ? 0.0 : 
									(*m_spatial_variant_weight)(i0, j0, k0));				
*/
				float kappa2 = ( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(i0, j0, k0) );					


					float diff = voxel_j - voxel_k;
					
					s += diff * diff * m_geometrical_weight[n] * kappa1 * kappa2;

				}

				(*weight)(i, j, k) = dot_phi_over_t(sqrt(s)); 

			}
		}
	}
	
	return weight;	
}								  

void PairwiseMRFRegularizer::calculateTerm1AndTerm2Isotropic(const Image<float>& x,
										 					 const Image<float>& neighbor_voxel_position_offset,
										 					 const Image<float>& image_weight,
										 					 const float geometrical_weight,
										 					 Image<float>& term1,
										 					 Image<float>& term2)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	//
	for (int k = 0; k < nk; k ++) {

		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// voxel j
				float voxel_j = x(i, j, k);

				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 :
						   	   (*m_spatial_variant_weight)(i, j, k);
			
				// voxel k, neighbor of voxel j
				int i0 = i + neighbor_voxel_position_offset(0);
				int j0 = j + neighbor_voxel_position_offset(1);
				int k0 = k + neighbor_voxel_position_offset(2);
				float voxel_k = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : x(i0, j0, k0);
				float voxel_k1 = ((i0<0) || (i0>=ni) || 
								 (j0<0) || (j0>=nj) || 
								 (k0<0) || (k0>=nk)) ? 
								 0.0 : 1.0;

/*				float kappa2 = ((m_spatial_variant_weight == 0) || (i0<0) || (i0>=ni) || 
								(j0<0) || (j0>=nj) || 
								(k0<0) || (k0>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0);
*/
				// fixed!
/*								
				float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
								( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(i0, j0, k0));				
*/
				float kappa2 = ( ((i0<0) || (i0>=ni) || 
									(j0<0) || (j0>=nj) ||  
									(k0<0) || (k0>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(i0, j0, k0) );	
				
				// for weight;
				float diff = voxel_j - voxel_k;
				float sum = 1.0 + voxel_k1;
				float weight = image_weight(i, j, k) * geometrical_weight * kappa1 * kappa2;
				
				// insert back
				term2(i, j, k) += diff * weight;
				term1(i, j, k) += sum * weight;		

				if (!((i0<0) || (i0>=ni) || 
					  (j0<0) || (j0>=nj) || 
					  (k0<0) || (k0>=nk))) {
					term2(i0, j0, k0) -= diff * weight;
					term1(i0, j0, k0) += sum * weight;
				}
			}
		}
	}	
}										 

/*
void PairwiseMRFRegularizer::createDifferenceFilter(int neighborhood_size)
{
	if (neighborhood_size == 0) {
		m_difference_filter.resize(3);
		m_filter_weight.resize(3);
		for (size_t n = 0; n < m_difference_filter.size(); n ++) {
			m_difference_filter[n] = new Image<float>(3,3,3);
			m_difference_filter[n]->clear();
			m_filter_weight[n] = 1.0;
		}
		
		// 3 first-order diff
		(*m_difference_filter[0])(1,1,1) = 1;
		(*m_difference_filter[0])(2,1,1) = -1;
		//
		(*m_difference_filter[1])(1,1,1) = 1;
		(*m_difference_filter[1])(1,2,1) = -1;
		//
		(*m_difference_filter[2])(1,1,1) = 1;
		(*m_difference_filter[2])(1,1,2) = -1;
	}
	
	if (neighborhood_size == 1) {
		m_difference_filter.resize(13);
		m_filter_weight.resize(13);
		for (size_t n = 0; n < m_difference_filter.size(); n ++) {
			m_difference_filter[n] = new Image<float>(3,3,3);
			m_difference_filter[n]->clear();
		}
		// 3 first-order diff
		(*m_difference_filter[0])(1,1,1) = 1;
		(*m_difference_filter[0])(2,1,1) = -1;
		m_filter_weight[0] = 1.0;
		//
		(*m_difference_filter[1])(1,1,1) = 1;
		(*m_difference_filter[1])(1,2,1) = -1;
		m_filter_weight[1] = 1.0;
		//
		(*m_difference_filter[2])(1,1,1) = 1;
		(*m_difference_filter[2])(1,1,2) = -1;
		m_filter_weight[2] = 1.0;
		
		// 
		(*m_difference_filter[3])(1,1,1) = 1;
		(*m_difference_filter[3])(2,2,1) = -1;		
		m_filter_weight[3] = 1.0 / sqrt(2.0);		
		(*m_difference_filter[4])(1,1,1) = 1;
		(*m_difference_filter[4])(0,2,1) = -1;
		m_filter_weight[4] = 1.0 / sqrt(2.0);		
		
		//
		(*m_difference_filter[5])(1,1,1) = 1;
		(*m_difference_filter[5])(2,1,2) = -1;	
		m_filter_weight[5] = 1.0 / sqrt(2.0);	
		(*m_difference_filter[6])(1,1,1) = 1;
		(*m_difference_filter[6])(0,1,2) = -1;	
		m_filter_weight[6] = 1.0 / sqrt(2.0);			
		//
		(*m_difference_filter[7])(1,1,1) = 1;
		(*m_difference_filter[7])(1,2,2) = -1;	
		m_filter_weight[7] = 1.0 / sqrt(2.0);	
		(*m_difference_filter[8])(1,1,1) = 1;
		(*m_difference_filter[8])(1,0,2) = -1;
		m_filter_weight[8] = 1.0 / sqrt(2.0);						
		
		//
		(*m_difference_filter[9])(1,1,1) = 1;
		(*m_difference_filter[9])(2,2,0) = -1;		
		m_filter_weight[9] = 1.0 / sqrt(3.0);
		(*m_difference_filter[10])(1,1,1) = 1;
		(*m_difference_filter[10])(0,2,0) = -1;		
		m_filter_weight[10] = 1.0 / sqrt(3.0);
		(*m_difference_filter[11])(1,1,1) = 1;
		(*m_difference_filter[11])(2,2,2) = -1;		
		m_filter_weight[11] = 1.0 / sqrt(3.0);
		(*m_difference_filter[12])(1,1,1) = 1;
		(*m_difference_filter[12])(0,2,2) = -1;		
		m_filter_weight[12] = 1.0 / sqrt(3.0);
		
	}

	// adjoint 
	m_difference_filter_abs.resize(m_difference_filter.size());
	m_difference_filter_adjoint.resize(m_difference_filter.size());
	m_difference_filter_adjoint_abs.resize(m_difference_filter.size());
	for (size_t n = 0; n < m_difference_filter.size(); n ++) {
		int fi = m_difference_filter[n]->getDimI();
		int fj = m_difference_filter[n]->getDimJ();
		int fk = m_difference_filter[n]->getDimK();		
		m_difference_filter_abs[n] = new Image<float>(fi,fj,fk);
		m_difference_filter_adjoint[n] = new Image<float>(fi,fj,fk);
		m_difference_filter_adjoint_abs[n] = new Image<float>(fi,fj,fk);
		for (int k = 0; k < fk; k ++) {
			for (int j = 0; j < fj; j ++) {
				for (int i = 0; i < fi; i ++) {
					(*m_difference_filter_adjoint[n])(i,j,k) = 
						(*m_difference_filter[n])(fi-i-1, fj-j-1, fk-k-1);
					(*m_difference_filter_abs[n])(i,j,k) = 
						fabs((*m_difference_filter[n])(i,j,k));
					(*m_difference_filter_adjoint_abs[n])(i,j,k) = 
						fabs((*m_difference_filter_adjoint[n])(i,j,k));
				}
			}
		}
	}

	// check
#if 0		
	for (size_t n = 0; n < m_difference_filter.size(); n ++) {
		SystemLog::write("#%d \n", n+1);
		for (int k = 0; k < 3; k ++) {
			for (int j = 0; j < 3; j ++) {
				for (int i = 0; i < 3; i ++) {
					SystemLog::write("%+.4f, ", (*(m_difference_filter_adjoint[n]))(i,j,k));
				}
				SystemLog::write("\n");
			}
			SystemLog::write("\n");
		}
	}	
#endif

}
*/

/*
void PairwiseMRFRegularizer::calculateNeighborhoodDifference(const Image<float>& x, 
										 					 const Image<float>& flt, 
										 					 Image<float>& diff)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	int fi = flt.getDimI();
	int fj = flt.getDimJ();
	int fk = flt.getDimK();
	int mid = fi / 2;
	
	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
				
				float s = 0.0; // sum
				
				for (int kk = 0; kk < fk; kk++) {
					int k0 = k + kk - mid;
					if (k0 >= 0 && k0 < nk) {
					
						for (int jj = 0; jj < fj; jj++) {
							int j0 = j + jj - mid;
							
							if (j0 >= 0 && j0 < nj) {
								for (int ii = 0; ii < fi; ii++) {
									int i0 = i + ii - mid;
									if (i0 >= 0 && i0 < ni) {
										s += x(i0,j0,k0) * flt(ii,jj,kk);
									}
								}
							}
							
						}
					
					}
				}
				
				diff(i,j,k) = s;				
			}
			
		}
	}
}	
*/

/*
void PairwiseMRFRegularizer::calculateGradientBasedOnOptimzationTransfer(const Image<float>& x,
        Image<float>& t1,
        Image<float>& t2)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	Image<float> all_ones(ni, nj, nk);
	all_ones.set(1.0);
	Image<float> temp1(ni, nj, nk);
	Image<float> temp2(ni, nj, nk);
	
	if (!m_is_isotropic) { // anisotropic
				
		Image<float> v1(ni, nj, nk);
		Image<float> v2(ni, nj, nk);
		
		t1.clear();
		t2.clear();
		for (size_t n = 0; n < m_difference_filter.size(); n ++) {
		
			// C * x
			temp1.clear();
			calculateNeighborhoodDifference(x, *m_difference_filter[n], temp1);
			
			// |C| * 1
			temp2.clear();
			calculateNeighborhoodDifference(all_ones, *m_difference_filter_abs[n], temp2);

			for (int j = 0; j < x.getSize(); j ++) {
				v1[j] = m_filter_weight[n] * dot_phi_over_t(temp1[j]) * temp2[j]; //diag(weight)*|C|*1
				v2[j] = m_filter_weight[n] * dot_phi_over_t(temp1[j]) * temp1[j]; //diag(weight)*C * x
			}

			temp2.clear();
			calculateNeighborhoodDifference(v1, *m_difference_filter_adjoint_abs[n], temp2);
			temp1.clear();
			calculateNeighborhoodDifference(v2, *m_difference_filter_adjoint[n], temp1);
			
			for (int j = 0; j < x.getSize(); j ++) {
				t1[j] += temp2[j];
				t2[j] += temp1[j];
			}
			
		}
			
	} else { // isotropic
		
		Image<float> weight(ni, nj, nk);
		Image<float> v1(ni, nj, nk);
		Image<float> v2(ni, nj, nk);
		
		weight.clear();
		for (size_t n = 0; n < m_difference_filter.size(); n ++) {
			temp1.clear();
			calculateNeighborhoodDifference(x, *m_difference_filter[n], temp1);
			for (int j = 0; j < x.getSize(); j ++) {
				weight[j] += temp1[j];
			}
		}
		for (int j = 0; j < x.getSize(); j ++) {
			weight[j] += dot_phi_over_t(weight[j]);
		}
	
		for (size_t n = 0; n < m_difference_filter.size(); n ++) {
			
			temp1.clear();
			calculateNeighborhoodDifference(x, *m_difference_filter[n], temp1);
			temp2.clear();
			calculateNeighborhoodDifference(all_ones, *m_difference_filter_abs[n], temp2);

			for (int j = 0; j < x.getSize(); j ++) {
				v1[j] = weight[j] * temp2[j]; //diag(weight)*|C|*1
				v2[j] = weight[j] * temp1[j]; //diag(weight)* C * x
			}
			
			temp1.clear();
			calculateNeighborhoodDifference(v1, *m_difference_filter_adjoint_abs[n], temp1);
			temp2.clear();
			calculateNeighborhoodDifference(v2, *m_difference_filter_adjoint[n], temp2);
			
			for (int j = 0; j < x.getSize(); j ++) {
				t1[j] += temp1[j];
				t2[j] += temp2[j];
			}

		}
	
	}
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
PatchBasedRegularizer::PatchBasedRegularizer() :
m_patch_size(-1)
{
}

PatchBasedRegularizer::~PatchBasedRegularizer()
{
	for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
		if (m_neighbor_patch_position_offset[n] != 0) {
			delete m_neighbor_patch_position_offset[n];
		}
	}
}

void PatchBasedRegularizer::initialize(int potential_function_type,
    		   		int neighborhood_size,
    		   		bool is_isotropic,
    		   		double* buildin_parms,
    		   		int number_of_buildin_parms_detected)
{
	SystemLog::write("initializing patch-based regularizer ...\n");

	Regularizer::initialize(potential_function_type, 
							neighborhood_size, 
							is_isotropic, 
							buildin_parms, 
							number_of_buildin_parms_detected);
							
	switch (neighborhood_size) {
		case 0:
			m_patch_size = 3; // meaning 3 x 3 x 3 (direct neighbor patches only)
			break;
		case 1:
			m_patch_size = 3; // 3 x 3 x 3 (more neighbors)
			break;
		default:
			m_patch_size = 3; // 3 x 3 x 3 (more neighbors)
			SystemLog::write("large neighborhood size chosen! but patch size not change, still 3 x 3 x 3!\n");
//			abort();
			break;
	}
	
	// position vectors
	SystemLog::write("creating offsets ...\n");
	createNeighborPatchPositionOffset();	

	SystemLog::write("[pot:%d, ns:%d, patch_size:%d, iso:%d], ok.\n",
		m_potential_function_type, m_neighborhood_size, m_patch_size, m_is_isotropic);
	
	if (is_isotropic) {
		SystemLog::write("note: isotropic form is enabled.\n");
		SystemLog::write("gradient of original penalty not implemented yet!\n");
	}
}

void PatchBasedRegularizer::createNeighborPatchPositionOffset()
{
	if (m_neighborhood_size == 0) {
	
		m_neighbor_patch_position_offset.resize(6);		
		for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
			m_neighbor_patch_position_offset[n] = new Image<float>(6);
			m_neighbor_patch_position_offset[n]->clear();
		}
		
		// first-order neighbor
		(*m_neighbor_patch_position_offset[0])(0) = 1; // i+1, j, k
		(*m_neighbor_patch_position_offset[1])(1) = 1; // i, j+1, k
		(*m_neighbor_patch_position_offset[2])(2) = 1; // i, j, k+1
		(*m_neighbor_patch_position_offset[3])(0) = -1; // i-1, j, k
		(*m_neighbor_patch_position_offset[4])(1) = -1; // i, j-1, k
		(*m_neighbor_patch_position_offset[5])(2) = -1; // i, j, k-1
		
	} else {
	
		if (m_neighborhood_size == 1) {
			
			m_neighbor_patch_position_offset.resize(26);		
			for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
				m_neighbor_patch_position_offset[n] = new Image<float>(3);
				m_neighbor_patch_position_offset[n]->clear();
			}
		
			// first-order neighbor
			(*m_neighbor_patch_position_offset[0])(0) = 1; // i+1, j, k
			(*m_neighbor_patch_position_offset[1])(1) = 1; // i, j+1, k
			(*m_neighbor_patch_position_offset[2])(2) = 1; // i, j, k+1
			// diagonals
			(*m_neighbor_patch_position_offset[3])(0) = 1;
			(*m_neighbor_patch_position_offset[3])(1) = 1;
			(*m_neighbor_patch_position_offset[4])(0) = 1;
			(*m_neighbor_patch_position_offset[4])(1) = -1;
			//
			(*m_neighbor_patch_position_offset[5])(1) = 1;
			(*m_neighbor_patch_position_offset[5])(2) = 1;
			(*m_neighbor_patch_position_offset[6])(1) = 1;
			(*m_neighbor_patch_position_offset[6])(2) = -1;
			//
			(*m_neighbor_patch_position_offset[7])(2) = 1;
			(*m_neighbor_patch_position_offset[7])(0) = 1;
			(*m_neighbor_patch_position_offset[8])(2) = 1;
			(*m_neighbor_patch_position_offset[8])(0) = -1;
		
			// more 
			(*m_neighbor_patch_position_offset[9])(0) = 1;
			(*m_neighbor_patch_position_offset[9])(1) = 1;
			(*m_neighbor_patch_position_offset[9])(2) = 1;
			//
			(*m_neighbor_patch_position_offset[10])(0) = 1;
			(*m_neighbor_patch_position_offset[10])(1) = 1;
			(*m_neighbor_patch_position_offset[10])(2) = -1;
			//
			(*m_neighbor_patch_position_offset[11])(0) = 1;
			(*m_neighbor_patch_position_offset[11])(1) = -1;
			(*m_neighbor_patch_position_offset[11])(2) = 1;
			//
			(*m_neighbor_patch_position_offset[12])(0) = 1;
			(*m_neighbor_patch_position_offset[12])(1) = -1;
			(*m_neighbor_patch_position_offset[12])(2) = -1;

			// another set of patches (added 04-21-2014)
			// first-order neighbor
			(*m_neighbor_patch_position_offset[0 + 13])(0) = -1; // i-1, j, k
			(*m_neighbor_patch_position_offset[1 + 13])(1) = -1; // i, j-1, k
			(*m_neighbor_patch_position_offset[2 + 13])(2) = -1; // i, j, k-1
			// diagonals
			(*m_neighbor_patch_position_offset[3 + 13])(0) = -1;
			(*m_neighbor_patch_position_offset[3 + 13])(1) = -1;
			(*m_neighbor_patch_position_offset[4 + 13])(0) = -1;
			(*m_neighbor_patch_position_offset[4 + 13])(1) = 1;
			//
			(*m_neighbor_patch_position_offset[5 + 13])(1) = -1;
			(*m_neighbor_patch_position_offset[5 + 13])(2) = -1;
			(*m_neighbor_patch_position_offset[6 + 13])(1) = -1;
			(*m_neighbor_patch_position_offset[6 + 13])(2) = 1;
			//
			(*m_neighbor_patch_position_offset[7 + 13])(2) = -1;
			(*m_neighbor_patch_position_offset[7 + 13])(0) = -1;
			(*m_neighbor_patch_position_offset[8 + 13])(2) = -1;
			(*m_neighbor_patch_position_offset[8 + 13])(0) = 1;
		
			// more 
			(*m_neighbor_patch_position_offset[9 + 13])(0) = -1;
			(*m_neighbor_patch_position_offset[9 + 13])(1) = -1;
			(*m_neighbor_patch_position_offset[9 + 13])(2) = -1;
			//
			(*m_neighbor_patch_position_offset[10 + 13])(0) = -1;
			(*m_neighbor_patch_position_offset[10 + 13])(1) = -1;
			(*m_neighbor_patch_position_offset[10 + 13])(2) = 1;
			//
			(*m_neighbor_patch_position_offset[11 + 13])(0) = -1;
			(*m_neighbor_patch_position_offset[11 + 13])(1) = 1;
			(*m_neighbor_patch_position_offset[11 + 13])(2) = -1;
			//
			(*m_neighbor_patch_position_offset[12 + 13])(0) = -1;
			(*m_neighbor_patch_position_offset[12 + 13])(1) = 1;
			(*m_neighbor_patch_position_offset[12 + 13])(2) = 1;
		
		} else {

			if (m_neighborhood_size >= 2) {

				// added by jnzhou@06-12-2014
				int hw = m_neighborhood_size;
				int ps = 2*m_neighborhood_size + 1;
				m_neighbor_patch_position_offset.resize(ps * ps * ps - 1); 
				for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
					m_neighbor_patch_position_offset[n] = new Image<float>(3);
					m_neighbor_patch_position_offset[n]->clear();
				}

				int c = 0;
				for (int k = -hw; k <= hw; k ++) {
					for (int j = -hw; j <= hw; j ++) {
						for (int i = -hw; i <= hw; i ++) {
							if ((i == 0) && (j == 0) && (k == 0)) {
								continue;
							}
							(*m_neighbor_patch_position_offset[c])(0) = i;
							(*m_neighbor_patch_position_offset[c])(1) = j;							
							(*m_neighbor_patch_position_offset[c])(2) = k;
							c ++;
						}	
					}					
				}

			} else {

				SystemLog::write("neighbor size not supported yet. "
							 "valid sizes are 0, 1 (3 x 3 x 3), 2 (5 x 5 x 5), ...\n");
				abort();
			}
		}
	
	}
	
}
									 					 
void PatchBasedRegularizer::calculateGradientBasedOnOptimizationTransfer(const Image<float>& x,
        																Image<float>& t1,
        																Image<float>& t2)
{
	// anisotropic case only
#if USE_OMP

	std::vector<Image<float>* > t1_list;
	std::vector<Image<float>* > t2_list;
	t1_list.resize(m_neighbor_patch_position_offset.size());
	t2_list.resize(m_neighbor_patch_position_offset.size());

	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	if (!m_is_isotropic) {

		size_t n = 0;
		#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
			t1_list[n] = new Image<float>(ni, nj, nk);
			t2_list[n] = new Image<float>(ni, nj, nk);
			calculateTerm1AndTerm2Anisotropic(x, *m_neighbor_patch_position_offset[n], 
											  *t1_list[n], *t2_list[n]);
		}
	
	} else {

		size_t n = 0;
		Image<float>* weight = calculateWeightIsotropic(x);

		#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
			t1_list[n] = new Image<float>(ni, nj, nk);
			t2_list[n] = new Image<float>(ni, nj, nk);
			calculateTerm1AndTerm2Isotropic(x, *m_neighbor_patch_position_offset[n], 
											*weight, *t1_list[n], *t2_list[n]);
		}

		delete weight;
	}

	for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
		int j = 0;
		#pragma omp parallel for private(j) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (j = 0; j < x.getSize(); j ++) {
			t1[j] += (*t1_list[n])[j];
			t2[j] += (*t2_list[n])[j];
		}
		delete t1_list[n];
		delete t2_list[n];
	}

#else
	
	if (!m_is_isotropic) {
		for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
			calculateTerm1AndTerm2Anisotropic(x, *m_neighbor_patch_position_offset[n], t1, t2);
		}
	} else {
		Image<float>* weight = calculateWeightIsotropic(x);
		for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
			calculateTerm1AndTerm2Isotropic(x, *m_neighbor_patch_position_offset[n], 
											*weight, t1, t2);
		}
		delete weight;
	}

#endif

}

float PatchBasedRegularizer::calculateWeightIsotropicSingleVoxel(const Image<float>& x, const int voxel_index)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	int nh = m_patch_size / 2;

	Image<float> patch_j = Image<float>(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_k = Image<float>(m_patch_size, m_patch_size, m_patch_size);	

	size_t i, j, k;
	x.ind2sub(size_t(voxel_index), i, j, k);

	// patch j
	for (int kk=0; kk<m_patch_size; kk++) {
		for (int jj=0; jj<m_patch_size; jj++) {
			for (int ii=0; ii<m_patch_size; ii++) {
						
				int i0 = i - nh + ii;
				int j0 = j - nh + jj;
				int k0 = k - nh + kk;
							
				patch_j(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
									 (j0<0) || (j0>=nj) || 
									 (k0<0) || (k0>=nk)) ? 
									 0.0 : x(i0, j0, k0);							
			}
		}
	}

	float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 : 
					(*m_spatial_variant_weight)(i, j, k);


	// patch k, neighborhood patch j
	float sum_all = 0.0;
	for (int nnn = 0; nnn < m_neighbor_patch_position_offset.size(); nnn ++) {

		Image<float>* neighbor_patch_position_offset =
			m_neighbor_patch_position_offset[nnn];	

		float s = 0.0;
		for (int kk=0; kk<m_patch_size; kk++) {
			for (int jj=0; jj<m_patch_size; jj++) {
				for (int ii=0; ii<m_patch_size; ii++) {
						
						int i0 = i - nh + ii + (*neighbor_patch_position_offset)(0);
						int j0 = j - nh + jj + (*neighbor_patch_position_offset)(1);
						int k0 = k - nh + kk + (*neighbor_patch_position_offset)(2);
							
						patch_k(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
											 (j0<0) || (j0>=nj) || 
											 (k0<0) || (k0>=nk)) ? 
											 0.0 : x(i0, j0, k0);
														
						// l2-norm
						float df = patch_j(ii,jj,kk) - patch_k(ii,jj,kk);							
						s += df * df;							
							
				}
			}
		}

		int iii = i + (*neighbor_patch_position_offset)(0);
		int jjj = j + (*neighbor_patch_position_offset)(1);
		int kkk = k + (*neighbor_patch_position_offset)(2);	
		/*			
		float kappa2 = ((m_spatial_variant_weight == 0) || (iii<0) || (iii>=ni) || 
						(jjj<0) || (jjj>=nj) || 
						(kkk<0) || (kkk>=nk)) ? 1.0 : 
						(*m_spatial_variant_weight)(iii, jjj, kkk);
						*/
/*
		float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
					( ((iii<0) || (iii>=ni) || 
							(jjj<0) || (jjj>=nj) ||  
						(kkk<0) || (kkk>=nk)) ? 0.0 : 
					(*m_spatial_variant_weight)(iii, jjj, kkk));
*/

		float kappa2 = ( ((iii<0) || (iii>=ni) || 
									(jjj<0) || (jjj>=nj) ||  
									(kkk<0) || (kkk>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(iii, jjj, kkk) );					


		// sum together;
		sum_all += s * kappa1 * kappa2;
	}
	return dot_phi_over_t(sqrt(sum_all));
}

Image<float>* PatchBasedRegularizer::calculateWeightIsotropic(const Image<float>& x)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	Image<float>* weight = new Image<float>(ni, nj, nk);

#if 0	
	int nh = m_patch_size / 2;

	Image<float> patch_j = Image<float>(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_k = Image<float>(m_patch_size, m_patch_size, m_patch_size);

	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// patch j
				for (int kk=0; kk<m_patch_size; kk++) {
					for (int jj=0; jj<m_patch_size; jj++) {
						for (int ii=0; ii<m_patch_size; ii++) {
						
							int i0 = i - nh + ii;
							int j0 = j - nh + jj;
							int k0 = k - nh + kk;
							
							patch_j(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : x(i0, j0, k0);							
						}
					}
				}

				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 : 
								(*m_spatial_variant_weight)(i, j, k);


				// patch k, neighborhood patch j
				float sum_all = 0.0;
				for (int nnn = 0; nnn < m_neighbor_patch_position_offset.size(); nnn ++) {

					Image<float>* neighbor_patch_position_offset =
						m_neighbor_patch_position_offset[nnn];	

					float s = 0.0;
					for (int kk=0; kk<m_patch_size; kk++) {
						for (int jj=0; jj<m_patch_size; jj++) {
							for (int ii=0; ii<m_patch_size; ii++) {
						
								int i0 = i - nh + ii + (*neighbor_patch_position_offset)(0);
								int j0 = j - nh + jj + (*neighbor_patch_position_offset)(1);
								int k0 = k - nh + kk + (*neighbor_patch_position_offset)(2);
							
								patch_k(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
													 (j0<0) || (j0>=nj) || 
													 (k0<0) || (k0>=nk)) ? 
													 0.0 : x(i0, j0, k0);
														
								// l2-norm
								float df = patch_j(ii,jj,kk) - patch_k(ii,jj,kk);							
								s += df * df;							
							
							}
						}
					}

					int iii = i + (*neighbor_patch_position_offset)(0);
					int jjj = j + (*neighbor_patch_position_offset)(1);
					int kkk = k + (*neighbor_patch_position_offset)(2);				
					float kappa2 = ((m_spatial_variant_weight == 0) || (iii<0) || (iii>=ni) || 
								(jjj<0) || (jjj>=nj) || 
								(kkk<0) || (kkk>=nk)) ? 1.0 : 
								(*m_spatial_variant_weight)(iii, jjj, kkk);

					// sum together;
					sum_all += s * kappa1 * kappa2;
				}

				(*weight)(i,j,k) = dot_phi_over_t(sqrt(sum_all));

			}
		}
	}
#else // fast version
	int n = 0;
	#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
	for (n = 0; n < x.getSize(); n ++) {
		(*weight)[n] = calculateWeightIsotropicSingleVoxel(x, n);
	}
#endif

	return weight;	
}

void PatchBasedRegularizer::calculateTerm1AndTerm2Isotropic(const Image<float>& x,
										 					const Image<float>& neighbor_patch_position_offset,
										 					const Image<float>& image_weight,
										 					Image<float>& term1,
										 					Image<float>& term2)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	int nh = m_patch_size / 2;
	Image<float> patch_j(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_k(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_diff(m_patch_size, m_patch_size, m_patch_size);
	
	//
	Image<float> patch_j1(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_k1(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_sum(m_patch_size, m_patch_size, m_patch_size);

	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// patch j
				for (int kk=0; kk<m_patch_size; kk++) {
					for (int jj=0; jj<m_patch_size; jj++) {
						for (int ii=0; ii<m_patch_size; ii++) {
						
							int i0 = i - nh + ii;
							int j0 = j - nh + jj;
							int k0 = k - nh + kk;
							
							patch_j(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : x(i0, j0, k0);
							
							patch_j1(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : 1.0;					 
							
						}
					}
				}
				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 : 
								(*m_spatial_variant_weight)(i, j, k);
			
				// patch k, neighborhood patch j
				float s = 0.0;
				for (int kk=0; kk<m_patch_size; kk++) {
					for (int jj=0; jj<m_patch_size; jj++) {
						for (int ii=0; ii<m_patch_size; ii++) {
						
							int i0 = i - nh + ii + neighbor_patch_position_offset(0);
							int j0 = j - nh + jj + neighbor_patch_position_offset(1);
							int k0 = k - nh + kk + neighbor_patch_position_offset(2);
							
							patch_k(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : x(i0, j0, k0);
							
							patch_k1(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : 1.0;					 
							
							// l2-norm
							float df = patch_j(ii,jj,kk) - patch_k(ii,jj,kk);							
							s += df * df;							
							patch_diff(ii,jj,kk) = df;
							
							float sm = patch_j1(ii,jj,kk) + patch_k1(ii,jj,kk);
							patch_sum(ii,jj,kk) = sm;
						}
					}
				}
				
				int iii = i + neighbor_patch_position_offset(0);
				int jjj = j + neighbor_patch_position_offset(1);
				int kkk = k + neighbor_patch_position_offset(2);
				/*				
				float kappa2 = ((m_spatial_variant_weight == 0) || (iii<0) || (iii>=ni) || 
								(jjj<0) || (jjj>=nj) || 
								(kkk<0) || (kkk>=nk)) ? 1.0 : 
								(*m_spatial_variant_weight)(iii, jjj, kkk);*/
/*
				float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
								( ((iii<0) || (iii>=ni) || 
									(jjj<0) || (jjj>=nj) ||  
									(kkk<0) || (kkk>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(iii, jjj, kkk));
*/

		float kappa2 = ( ((iii<0) || (iii>=ni) || 
									(jjj<0) || (jjj>=nj) ||  
									(kkk<0) || (kkk>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(iii, jjj, kkk) );		

				// compute weight
				// only difference between aniso and iso form
				// image_weight is precomputed by calculateWeightIsotropic()
				float weight = image_weight(i, j, k) * kappa1 * kappa2;
				
				// put patch back to image				
				for (int kk=0; kk<m_patch_size; kk++) {
					int k0 = k - nh + kk;
					if (k0 >= 0 && k0 < nk) {
						for (int jj=0; jj<m_patch_size; jj++) {
							int j0 = j - nh + jj;
							if (j0 >= 0 && j0 < nj) {
								for (int ii=0; ii<m_patch_size; ii++) {							
									int i0 = i - nh + ii;									
									if (i0 >= 0 && i0 < ni) {
										term2(i0, j0, k0) += patch_diff(ii, jj, kk) * weight;
										term1(i0, j0, k0) += patch_sum(ii, jj, kk) * weight;
										//
									}
								}
							}
						}
					}
				}
				
				// another patch, neighbor patch, so take `-'
				for (int kk=0; kk<m_patch_size; kk++) {
					int k0 = k - nh + kk + neighbor_patch_position_offset(2);
					if (k0 >= 0 && k0 < nk) {
						for (int jj=0; jj<m_patch_size; jj++) {
							int j0 = j - nh + jj + neighbor_patch_position_offset(1);
							if (j0 >= 0 && j0 < nj) {
								for (int ii=0; ii<m_patch_size; ii++) {
									int i0 = i - nh + ii + neighbor_patch_position_offset(0);						
									if (i0 >= 0 && i0 < ni) {							
										term2(i0, j0, k0) -= patch_diff(ii, jj, kk) * weight;
										term1(i0, j0, k0) += patch_sum(ii, jj, kk) * weight;
									}
									
								}
							}
						}
					}
				}
			}
		}
	}	
}

void PatchBasedRegularizer::calculateTerm1AndTerm2Anisotropic(const Image<float>& x,
										   		   			  Image<float>& neighbor_patch_position_offset,
										   		   			  Image<float>& term1,
										   		   			  Image<float>& term2)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	int nh = m_patch_size / 2;
	Image<float> patch_j(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_k(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_diff(m_patch_size, m_patch_size, m_patch_size);
	
	//
	Image<float> patch_j1(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_k1(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_sum(m_patch_size, m_patch_size, m_patch_size);

	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// patch j
				for (int kk=0; kk<m_patch_size; kk++) {
					for (int jj=0; jj<m_patch_size; jj++) {
						for (int ii=0; ii<m_patch_size; ii++) {
						
							int i0 = i - nh + ii;
							int j0 = j - nh + jj;
							int k0 = k - nh + kk;
							
							patch_j(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : x(i0, j0, k0);
							
							patch_j1(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : 1.0;					 
							
						}
					}
				}
				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 : 
								(*m_spatial_variant_weight)(i, j, k);
			
				// patch k, neighborhood patch j
				float s = 0.0;
				for (int kk=0; kk<m_patch_size; kk++) {
					for (int jj=0; jj<m_patch_size; jj++) {
						for (int ii=0; ii<m_patch_size; ii++) {
						
							int i0 = i - nh + ii + neighbor_patch_position_offset(0);
							int j0 = j - nh + jj + neighbor_patch_position_offset(1);
							int k0 = k - nh + kk + neighbor_patch_position_offset(2);
							
							patch_k(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : x(i0, j0, k0);
							
							patch_k1(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : 1.0;					 
							
							// l2-norm
							float df = patch_j(ii,jj,kk) - patch_k(ii,jj,kk);							
							s += df * df;							
							patch_diff(ii,jj,kk) = df;
							
							float sm = patch_j1(ii,jj,kk) + patch_k1(ii,jj,kk);
							patch_sum(ii,jj,kk) = sm;
						}
					}
				}
				
				int iii = i + neighbor_patch_position_offset(0);
				int jjj = j + neighbor_patch_position_offset(1);
				int kkk = k + neighbor_patch_position_offset(2);	

				/*			
				float kappa2 = ((m_spatial_variant_weight == 0) || (iii<0) || (iii>=ni) || 
								(jjj<0) || (jjj>=nj) || 
								(kkk<0) || (kkk>=nk)) ? 1.0 : 
								(*m_spatial_variant_weight)(iii, jjj, kkk);*/
/*
				float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
								( ((iii<0) || (iii>=ni) || 
									(jjj<0) || (jjj>=nj) ||  
									(kkk<0) || (kkk>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(iii, jjj, kkk));
*/

		float kappa2 = ( ((iii<0) || (iii>=ni) || 
									(jjj<0) || (jjj>=nj) ||  
									(kkk<0) || (kkk>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(iii, jjj, kkk) );		

				// for weight;
				float weight = dot_phi_over_t(sqrt(s)) * kappa1 * kappa2;
								
				for (int kk=0; kk<m_patch_size; kk++) {
					int k0 = k - nh + kk;
					if (k0 >= 0 && k0 < nk) {
						for (int jj=0; jj<m_patch_size; jj++) {
							int j0 = j - nh + jj;
							if (j0 >= 0 && j0 < nj) {
								for (int ii=0; ii<m_patch_size; ii++) {							
									int i0 = i - nh + ii;									
									if (i0 >= 0 && i0 < ni) {
										term2(i0, j0, k0) += patch_diff(ii, jj, kk) * weight;
										term1(i0, j0, k0) += patch_sum(ii, jj, kk) * weight;
										//
									}
								}
							}
						}
					}
				}
				
				for (int kk=0; kk<m_patch_size; kk++) {
					int k0 = k - nh + kk + neighbor_patch_position_offset(2);
					if (k0 >= 0 && k0 < nk) {
						for (int jj=0; jj<m_patch_size; jj++) {
							int j0 = j - nh + jj + neighbor_patch_position_offset(1);
							if (j0 >= 0 && j0 < nj) {
								for (int ii=0; ii<m_patch_size; ii++) {
									int i0 = i - nh + ii + neighbor_patch_position_offset(0);						
									if (i0 >= 0 && i0 < ni) {							
										term2(i0, j0, k0) -= patch_diff(ii, jj, kk) * weight;
										term1(i0, j0, k0) += patch_sum(ii, jj, kk) * weight;
									}
									
								}
							}
						}
					}
				}
			}
		}
	}
}




// // Mengxi's version
// void PatchBasedRegularizer::calculateTerm1AndTerm2AnisotropicCIP(const Image<float>& x,
// 										   		   			  Image<int>& neighbor_patch_position_offset,
// 										   		   			  Image<float>& term1,
// 										   		   			  Image<float>& term2)
// {
//     // <Guobao>
//     int gbi, gbj, gbk, fdi;
    
// 	int ni = x.getDimI();
// 	int nj = x.getDimJ();
// 	int nk = x.getDimK();
    

//     int ciplabel = 1;
//     if (neighbor_patch_position_offset.getDimJ()==1)
//     {
//         ciplabel = 0;
//     }
    
// //  <Guobao>   
// //  @guobao use term1 to transfer cip image into function
// //  and remember to clear term1. 
//     Image<float> xcip(ni,nj,nk);
//     for (int k = 0; k < nk; k ++) {
//         for (int j = 0; j < nj; j ++) {
//             for (int i = 0; i < ni; i ++) {
//                 xcip(i,j,k) = term1(i,j,k);
//                 term1(i,j,k) = 0;
//             }}}
    
//     Image<float> gb_patch_j(m_patch_size, m_patch_size, m_patch_size);
// 	Image<float> gb_patch_k(m_patch_size, m_patch_size, m_patch_size);
//     // </Guobao>
	
// 	int nh = m_patch_size / 2;
// 	Image<float> patch_j(m_patch_size, m_patch_size, m_patch_size);
// 	Image<float> patch_k(m_patch_size, m_patch_size, m_patch_size);
// 	Image<float> patch_diff(m_patch_size, m_patch_size, m_patch_size);
	
	
// 	Image<float> patch_j1(m_patch_size, m_patch_size, m_patch_size);
// 	Image<float> patch_k1(m_patch_size, m_patch_size, m_patch_size);
// 	Image<float> patch_sum(m_patch_size, m_patch_size, m_patch_size);

// 	for (int k = 0; k < nk; k ++) {
	
// 		for (int j = 0; j < nj; j ++) {
			
// 			for (int i = 0; i < ni; i ++) {
                
// //                 <Guobao>
//                 if (ciplabel == 0)
//                 {
//                     gbi = neighbor_patch_position_offset(0);
//                     gbj = neighbor_patch_position_offset(1);
//                     gbk = neighbor_patch_position_offset(2);
//                 }
//                 else
//                 {
//                     fdi = k*ni*nj + j*ni +i;
//                     gbi = neighbor_patch_position_offset(fdi,0);
//                     gbj = neighbor_patch_position_offset(fdi,1);
//                     gbk = neighbor_patch_position_offset(fdi,2);
//                 }
                
// //                 if (i==132 && j==139 && k==61)
// //                 {
// //                     printf("%d %d %d\n",gbi,gbj,gbk);
// //                 }
                
// //                 </Guobao>
			
// 				// patch j
// 				for (int kk=0; kk<m_patch_size; kk++) {
// 					for (int jj=0; jj<m_patch_size; jj++) {
// 						for (int ii=0; ii<m_patch_size; ii++) {
						
// 							int i0 = i - nh + ii;
// 							int j0 = j - nh + jj;
// 							int k0 = k - nh + kk;
							
// 							patch_j(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
// 												 (j0<0) || (j0>=nj) || 
// 												 (k0<0) || (k0>=nk)) ? 
// 												 0.0 : x(i0, j0, k0);
                                                 
//                             gb_patch_j(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
// 												 (j0<0) || (j0>=nj) || 
// 												 (k0<0) || (k0>=nk)) ? 
// 												 0.0 : xcip(i0, j0, k0);                     
                                                 
							
// 							patch_j1(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
// 												 (j0<0) || (j0>=nj) || 
// 												 (k0<0) || (k0>=nk)) ? 
// 												 0.0 : 1.0;					 
							
// 						}
// 					}
// 				}
// 				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 : 
// 								(*m_spatial_variant_weight)(i, j, k);
			
// 				// patch k, neighborhood patch j
// 				float s = 0.0;
//                 float gb_s = 0;
// 				for (int kk=0; kk<m_patch_size; kk++) {
// 					for (int jj=0; jj<m_patch_size; jj++) {
// 						for (int ii=0; ii<m_patch_size; ii++) {
						
// 							int i0 = i - nh + ii + gbi;
// 							int j0 = j - nh + jj + gbj;
// 							int k0 = k - nh + kk + gbk;
							
// 							patch_k(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
// 												 (j0<0) || (j0>=nj) || 
// 												 (k0<0) || (k0>=nk)) ? 
// 												 0.0 : x(i0, j0, k0);
                                                 
//                             gb_patch_k(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
// 												 (j0<0) || (j0>=nj) || 
// 												 (k0<0) || (k0>=nk)) ? 
// 												 0.0 : xcip(i0, j0, k0);                     
							
// 							patch_k1(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
// 												 (j0<0) || (j0>=nj) || 
// 												 (k0<0) || (k0>=nk)) ? 
// 												 0.0 : 1.0;					 
							
// 							// l2-norm
// 							float df = patch_j(ii,jj,kk) - patch_k(ii,jj,kk);							
// 							s += df * df;							
// 							patch_diff(ii,jj,kk) = df;
                           
// 							float sm = patch_j1(ii,jj,kk) + patch_k1(ii,jj,kk);
// 							patch_sum(ii,jj,kk) = sm;
                            
// //                             @ guobao
//                             gb_s += (gb_patch_j(ii,jj,kk) - gb_patch_k(ii,jj,kk)) * (gb_patch_j(ii,jj,kk) - gb_patch_k(ii,jj,kk));
                            
// 						}
// 					}
// 				}
				
// 				int iii = i + gbi;
// 				int jjj = j + gbj;
// 				int kkk = k + gbk;

// 		float kappa2 = ( ((iii<0) || (iii>=ni) || 
// 						 (jjj<0) || (jjj>=nj) ||  
// 						 (kkk<0) || (kkk>=nk)) ) ? 0.0 : 
// 							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
// 							  	(*m_spatial_variant_weight)(iii, jjj, kkk) );

// 				// for weight;
//                 float weight = dot_phi_over_t(sqrt(s)) * kappa1 * kappa2;
                
// //                 @ guobao
//                 if(ciplabel == 1)
//                 {
// //                     gb_s = sqrt(gb_s);
//                     weight = weight * exp( (-gb_s)/3/pow(10,-9) );
// //                     weight = weight*6/50;
//                 }
                                
								
// 				for (int kk=0; kk<m_patch_size; kk++) {
// 					int k0 = k - nh + kk;
// 					if (k0 >= 0 && k0 < nk) {
// 						for (int jj=0; jj<m_patch_size; jj++) {
// 							int j0 = j - nh + jj;
// 							if (j0 >= 0 && j0 < nj) {
// 								for (int ii=0; ii<m_patch_size; ii++) {							
// 									int i0 = i - nh + ii;									
// 									if (i0 >= 0 && i0 < ni) {
// 										term2(i0, j0, k0) += patch_diff(ii, jj, kk) * weight;
// 										term1(i0, j0, k0) += patch_sum(ii, jj, kk) * weight;
// 										//
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
				
// 				for (int kk=0; kk<m_patch_size; kk++) {
// 					int k0 = k - nh + kk + neighbor_patch_position_offset(2);
// 					if (k0 >= 0 && k0 < nk) {
// 						for (int jj=0; jj<m_patch_size; jj++) {
// 							int j0 = j - nh + jj + neighbor_patch_position_offset(1);
// 							if (j0 >= 0 && j0 < nj) {
// 								for (int ii=0; ii<m_patch_size; ii++) {
// 									int i0 = i - nh + ii + neighbor_patch_position_offset(0);						
// 									if (i0 >= 0 && i0 < ni) {							
// 										term2(i0, j0, k0) -= patch_diff(ii, jj, kk) * weight;
// 										term1(i0, j0, k0) += patch_sum(ii, jj, kk) * weight;
// 									}
									
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// }



void PatchBasedRegularizer::calculateGradient(const Image<float>& x,
        									  Image<float>& grad)
{
	// anisotropic case only
#if USE_OMP

	std::vector<Image<float>* > grad_list;
	grad_list.resize(m_neighbor_patch_position_offset.size());

	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	size_t n = 0;
	#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
	for (n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
		grad_list[n] = new Image<float>(ni, nj, nk);
		calculateGradientAnisotropic(x, *m_neighbor_patch_position_offset[n], *grad_list[n]);
	}
	
	for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
		int j = 0;
		#pragma omp parallel for private(j) num_threads(Regularizer::NUMBER_OF_THREADS)
		for (j = 0; j < x.getSize(); j ++) {
			grad[j] += (*grad_list[n])[j];
		}
		delete grad_list[n];
	}

#else
	
	for (size_t n = 0; n < m_neighbor_patch_position_offset.size(); n ++) {
		calculateGradientAnisotropic(x, *m_neighbor_patch_position_offset[n], grad);
	}	

#endif

}

void PatchBasedRegularizer::calculateGradientAnisotropic(const Image<float>& x,
										   		   		 Image<float>& neighbor_patch_position_offset,
										   		   		 Image<float>& grad)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	int nh = m_patch_size / 2;
	Image<float> patch_j(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_k(m_patch_size, m_patch_size, m_patch_size);
	Image<float> patch_diff(m_patch_size, m_patch_size, m_patch_size);
	
	for (int k = 0; k < nk; k ++) {
	
		for (int j = 0; j < nj; j ++) {
			
			for (int i = 0; i < ni; i ++) {
			
				// patch j
				for (int kk=0; kk<m_patch_size; kk++) {
					for (int jj=0; jj<m_patch_size; jj++) {
						for (int ii=0; ii<m_patch_size; ii++) {
						
							int i0 = i - nh + ii;
							int j0 = j - nh + jj;
							int k0 = k - nh + kk;
							
							patch_j(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : x(i0, j0, k0);
							
						}
					}
				}
				float kappa1 = (m_spatial_variant_weight == 0) ? 1.0 : 
								(*m_spatial_variant_weight)(i, j, k);
			
				// patch k, neighborhood patch j
				float s = 0.0;
				for (int kk=0; kk<m_patch_size; kk++) {
					for (int jj=0; jj<m_patch_size; jj++) {
						for (int ii=0; ii<m_patch_size; ii++) {
						
							int i0 = i - nh + ii + neighbor_patch_position_offset(0);
							int j0 = j - nh + jj + neighbor_patch_position_offset(1);
							int k0 = k - nh + kk + neighbor_patch_position_offset(2);
							
							patch_k(ii,jj,kk) = ((i0<0) || (i0>=ni) || 
												 (j0<0) || (j0>=nj) || 
												 (k0<0) || (k0>=nk)) ? 
												 0.0 : x(i0, j0, k0);
							
							// l2-norm
							float df = patch_j(ii,jj,kk) - patch_k(ii,jj,kk);							
							s += df * df;							
							patch_diff(ii,jj,kk) = df;
							
						}
					}
				}
				
				int iii = i + neighbor_patch_position_offset(0);
				int jjj = j + neighbor_patch_position_offset(1);
				int kkk = k + neighbor_patch_position_offset(2);
				/*

				float kappa2 = ((m_spatial_variant_weight == 0) || (iii<0) || (iii>=ni) || 
								(jjj<0) || (jjj>=nj) || 
								(kkk<0) || (kkk>=nk)) ? 1.0 : 
								(*m_spatial_variant_weight)(iii, jjj, kkk);
*/
/*								
				float kappa2 = (m_spatial_variant_weight == 0) ? (1.0) : 
								( ((iii<0) || (iii>=ni) || 
									(jjj<0) || (jjj>=nj) ||  
									(kkk<0) || (kkk>=nk)) ? 0.0 : 
								(*m_spatial_variant_weight)(iii, jjj, kkk));
*/

		float kappa2 = ( ((iii<0) || (iii>=ni) || 
									(jjj<0) || (jjj>=nj) ||  
									(kkk<0) || (kkk>=nk)) ) ? 0.0 : 
							   ( (m_spatial_variant_weight == 0) ? 1.0 : 
							  	(*m_spatial_variant_weight)(iii, jjj, kkk) );							
				// for weight;
				float weight = dot_phi_over_t(sqrt(s)) * kappa1 * kappa2;
								
				for (int kk=0; kk<m_patch_size; kk++) {
					int k0 = k - nh + kk;
					if (k0 >= 0 && k0 < nk) {
						for (int jj=0; jj<m_patch_size; jj++) {
							int j0 = j - nh + jj;
							if (j0 >= 0 && j0 < nj) {
								for (int ii=0; ii<m_patch_size; ii++) {							
									int i0 = i - nh + ii;									
									if (i0 >= 0 && i0 < ni) {
										grad(i0, j0, k0) += patch_diff(ii, jj, kk) * weight;
										//
									}
								}
							}
						}
					}
				}
				
				for (int kk=0; kk<m_patch_size; kk++) {
					int k0 = k - nh + kk + neighbor_patch_position_offset(2);
					if (k0 >= 0 && k0 < nk) {
						for (int jj=0; jj<m_patch_size; jj++) {
							int j0 = j - nh + jj + neighbor_patch_position_offset(1);
							if (j0 >= 0 && j0 < nj) {
								for (int ii=0; ii<m_patch_size; ii++) {
									int i0 = i - nh + ii + neighbor_patch_position_offset(0);						
									if (i0 >= 0 && i0 < ni) {							
										grad(i0, j0, k0) -= patch_diff(ii, jj, kk) * weight;
									}
									
								}
							}
						}
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
MeanBasedRegularizer::MeanBasedRegularizer() :
m_mean_image(0)
{
}

MeanBasedRegularizer::~MeanBasedRegularizer()
{
	if (m_mean_image != 0) {
		delete m_mean_image;
	}
}

void MeanBasedRegularizer::initialize(int potential_function_type,
    		   					 	  int neighborhood_size,
    		   					 	  bool is_isotropic,
    		   					 	  double* buildin_parms,
    		   					 	  int number_of_buildin_parms_detected)
{
	Regularizer::initialize(potential_function_type, 
						    neighborhood_size, 
						    is_isotropic, 
						    buildin_parms, 
						    number_of_buildin_parms_detected);
}

void MeanBasedRegularizer::setMeanImage(const Image<float>& mean_image)
{
	int ni = mean_image.getDimI();
	int nj = mean_image.getDimJ();
	int nk = mean_image.getDimK();
	
	if (m_mean_image == 0) {
		m_mean_image = new Image<float>(ni, nj, nk);
	} else {
		if (m_mean_image->getSize() != mean_image.getSize()) {
			SystemLog::write("unmatched image size, cannot create mean image.\n");
			abort();
		}
	}
	
	memcpy(m_mean_image->getPtr(), 
		   mean_image.getPtr(), 
		   sizeof(float)*mean_image.getSize());
}

void MeanBasedRegularizer::calculateGradientBasedOnOptimizationTransfer(const Image<float>& x,
            										 				   Image<float>& term1,
            										 				   Image<float>& term2)
{
	if (m_mean_image == 0) {
		SystemLog::write("error, mean image not assigned!\n");
		abort();
	}
	
	for (int j = 0; j < x.getSize(); j ++) {
		term1[j] = 1.0;
		term2[j] = x[j] - (*m_mean_image)[j];
	}
}            										 				    





// // @Guobao useless. Just for compiler.
// void MeanBasedRegularizer::gbcalg(const Image<float>& x,
//      Image<float>& t1,Image<float>& t2,const Image<float>& cip,Image<int>* fd)
// {
//     printf("\nSorry, currently the code dose not support mean based CIP.\n");
//     abort();
// }
// void PairwiseMRFRegularizer::gbcalg(const Image<float>& x,
//         Image<float>& t1,Image<float>& t2,const Image<float>& cip,Image<int>* fd)
// {
//     printf("\nSorry, currently the code dose not support pairwise CIP.\n");
//     abort();
// }



// void PatchBasedRegularizer::gbcalg(const Image<float>& x,
//         								 Image<float>& t1,
//         								 Image<float>& t2,
//                                    const Image<float>& cip,
//                                          Image<int>* fd)
// {
// //     ref line 1434
    
//     int ni = x.getDimI();
// 	int nj = x.getDimJ();
// 	int nk = x.getDimK();
    
// //     printf("%d %d %d\n",ni,nj,nk);
    
//     int nvoxel = fd->getDimI();      // number of voxel (300x300x96 for phantom)
//     int nneighbor = fd->getDimJ();   // number of neighbor (50 (out of 7x7x7) by default)
// //     Image<int> fdsub(nvoxel,3);
    
//     std::vector<Image<int>* > fdsub;
//     fdsub.resize(nneighbor);
        
//     std::vector<Image<float>* > t1_list;
// 	std::vector<Image<float>* > t2_list;
//     t1_list.resize(nneighbor);
// 	t2_list.resize(nneighbor);

// // 	ref line 1850
// // 	if (!m_is_isotropic) {
// // 
// 		size_t n = 0;
// 		#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
// 		for (n = 0; n < nneighbor; n ++) {
// 			t1_list[n] = new Image<float>(ni, nj, nk);
// 			t2_list[n] = new Image<float>(ni, nj, nk);
//             fdsub[n] = new Image<int>(nvoxel,3);
            
            
            
//             for (int k = 0; k < nk; k ++) {
//                 for (int j = 0; j < nj; j ++) {
//                     for (int i = 0; i < ni; i ++) {
//                         (*t1_list[n])(i,j,k) = cip(i,j,k);
//                         for(int ii=0;ii<3;ii++)
//                         {
//                             (*fdsub[n])(k*ni*nj+j*ni+i,ii) = (*fd)(k*ni*nj+j*ni+i,n,ii);
//                         }
                        
//                     }}}

// 			calculateTerm1AndTerm2AnisotropicCIP(x, *fdsub[n], *t1_list[n], *t2_list[n]);
//             delete fdsub[n];
// 		}
// // 	
// // 	} else {
// // 
// // // 		size_t n = 0;       ??????
// // 		Image<float>* weight = calculateWeightIsotropic(x);
// // 
// // 		#pragma omp parallel for private(n) num_threads(Regularizer::NUMBER_OF_THREADS)
// // 		for (int n = 0; n < nneighbor; n ++) {
// // 			t1_list[n] = new Image<float>(ni, nj, nk);
// // 			t2_list[n] = new Image<float>(ni, nj, nk);
// // 			calculateTerm1AndTerm2Isotropic(x, *m_neighbor_patch_position_offset[n], 
// // 											*weight, *t1_list[n], *t2_list[n]);
// // 		}
// // 
// // 		delete weight;
// // 	}
// // 
// 	for (size_t n = 0; n < nneighbor; n ++) {
// 		int j = 0;
// 		#pragma omp parallel for private(j) num_threads(Regularizer::NUMBER_OF_THREADS)
// 		for (j = 0; j < x.getSize(); j ++) {
// 			t1[j] += (*t1_list[n])[j];
// 			t2[j] += (*t2_list[n])[j];
// 		}
// 		delete t1_list[n];
// 		delete t2_list[n];
        
// 	}
    
// }




