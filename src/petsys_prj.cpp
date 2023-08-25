#include <petsys_prj.h>

// int Projector::NUMBER_OF_THREADS_FP = 32;
// int Projector::NUMBER_OF_THREADS_BP = 32;


int Projector::NUMBER_OF_THREADS_FP = 48;
int Projector::NUMBER_OF_THREADS_BP = 48;


Projector::XtalPair Projector::createCrystalPairs(const int num_of_detblocks_t,
        									      const int xtal_array_size_t,
        									      const int num_of_projs_per_angle,
										          const bool compr_enabled,
										          const bool half_block_rotation)
{
    int nx = num_of_detblocks_t * xtal_array_size_t;
    int na = nx / 2;

    if (num_of_detblocks_t % 2 != 0) {
#if DEBUG
        SystemLog::write("odd number of detector blocks. don't know how to create crystal pairs");
#endif
        abort();
    }

#if DEBUG
    SystemLog::write("total # of crystals (per ring): %d, # of angles: %d\n", nx, na);
#endif

    int nb = num_of_projs_per_angle;

    if (num_of_projs_per_angle % 2 == 0) {
        nb ++;
#if DEBUG
        SystemLog::write("even projection number! ");
        SystemLog::write("add one more to make it odd: %d\n", nb);
#endif
    }

#if DEBUG
    SystemLog::write("sinogram size = [%d, %d(a)], ", nb, na);

    if (half_block_rotation) {
        SystemLog::write("assume half block rotation because angle offset isn't equal to 0\n");
    }

#endif

    std::vector<int> half_0;
    std::vector<int> half_1;

    if (half_block_rotation) {
        int id0 = 0;
        int id1 = na;
        int nr = na / 2 * 2 - 2;

        for (int i = 0; i < nr; i ++) {
            //
            int id = id0 + nr / 2 - i - 1;

            if (id < 0) id += nx;

            half_0.push_back(id);
            //
            id = id1 - nr / 2 + i;
            half_1.push_back(id);
        }
    } else {
        // id's start from 0 (different from old version)
        int id0 = xtal_array_size_t / 2;
        // pair of center LOR at the first angle
        bool odd = xtal_array_size_t % 2 != 0;
        int id1 = odd ? (nx / 2 + xtal_array_size_t / 2) :
                  (nx / 2 + xtal_array_size_t / 2 - 1);
        int id, nr = odd ? (na / 2 * 2 - 1) : (na / 2 * 2);

        for (int i = 0; i < nr; i ++) {
            //
            id = odd ? (id0 + nr / 2 - i) : (id0 + nr / 2 - 1 - i);

            if (id < 0) id += nx;

            half_0.push_back(id);
            //
            id = odd ? (id1 - nr / 2 + i) : (id1 - nr / 2 + 1 + i);
            half_1.push_back(id);
#if 0
            SystemLog::write("%d | %d\n", half_0.back(), half_1.back());
#endif
        }
    }

    // all possible pairs in the first angle
    std::pair<int, int> xp;
    XtalPair xp_first_angle;
    std::size_t c = 0;

    while (c < half_0.size()) {
        xp.first = half_0[c];
        xp.second = half_1[c];
        xp_first_angle.push_back(xp);

        if ((c + 1) < half_1.size()) {
            xp.first = half_0[c];
            xp.second = half_1[c + 1];
            xp_first_angle.push_back(xp);
        }

        c ++;
    }

    // check if required projection number is greater than the maximum number
    if (xp_first_angle.size() < std::size_t(nb)) {
#if DEBUG
        SystemLog::write("exceed the maximum number of projections per angle: %d (> %d)\n",
              nb, int(xp_first_angle.size()));
#endif
        abort();
    }

    // for compression
    int na_c;
    int sym_num;

    if (num_of_detblocks_t % 4 == 0) { // 8-fold
//		int r = (na-2)/2;
//		na_c = (r%2 == 0) ? (r/2 + 1) : (r/2 + 2);
        // or simply at below
        na_c = na / 4 + 1;
#if DEBUG
        SystemLog::write("[%d (a compr)]\n", na_c);
#endif
        sym_num = 8;
    } else { // 4-fold
        na_c = (na - 2) / 2 + 2; // or simply na/2+1
        sym_num = 4;
    }

    // generate pairs for all angles
    // note: rotate CCW, so ++
    XtalPair xtal_pairs;

    if (compr_enabled) {
        xtal_pairs.resize(na_c * (nb / 2 + 1));
    } else {
        xtal_pairs.resize(na * nb);
    }

    int k = 0;

    for (int i = 0; i < na; i ++) {
#if 0
        SystemLog::write("angle #%d:\n", i + 1);
#endif

        for (int j = 0; j < nb; j ++) {
            xp = xp_first_angle[xp_first_angle.size() / 2 - nb / 2 + j];

            int p0 = (xp.first + i) % nx;
            int p1 = (xp.second + i) % nx;

            if (compr_enabled) {
                if (i < na_c && j <= (nb / 2)) {
                    xtal_pairs[k].first = p0;	// circular shift, so take modulo
                    xtal_pairs[k].second = p1;
                    k ++;
                }
            } else {
                xtal_pairs[k].first = p0;	// circular shift, so take modulo
                xtal_pairs[k].second = p1;
                k ++;
            }

            if (i < na_c && j <= (nb / 2)) {
#if 0 // output compressed pairs	
                output_compr << p0 << " " << p1 << std::endl;
#endif
            }
        }
    }

    return xtal_pairs;
}

void Projector::initializeIPSFModel(const SIZE& image_size, 
						 	 		const char* trans_psf_file, 
						 	 		const char* axial_psf_file)
{
	if (m_ipsf_model != 0) {
		SystemLog::write("error, IPSF model is already created.\n");
		abort();
	} else {
		m_ipsf_model = new IPSFModel(image_size);
	    if (!m_ipsf_model->initialize(trans_psf_file, axial_psf_file)) {
	    	SystemLog::write("psf intialization error ... %s, %s\n", 
    		trans_psf_file, axial_psf_file);
	    }
	}
}




// void Projector::initializeKernelModel(const SIZE& image_size, 
//                                     const char* kernel_matrix_file)
// {
//     if (m_ikernel_model != 0) {
//         SystemLog::write("error, IKernel model is already created.\n");
//         abort();
//     } else {
//         m_ikernel_model = new IKernelModel(image_size);
//         if (!m_ikernel_model->initialize(kernel_matrix_file)) {
//             SystemLog::write("kernel intialization error ... %s\n", 
//             kernel_matrix_file);
//         }
//     }
// }










