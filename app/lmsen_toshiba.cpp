#include <petsys_cscanner.h>

#include <petsys_log.h>

int main(int argc, char* argv[])
{
	SystemLog::open("calsen.log");
			
#if USE_TOF
	SystemLog::write("error, TOF projector not supported yet.\n");
	abort();
#endif

	CScanner scanner;
	
	const char* cfg_file = argv[1];
	scanner.initialize(cfg_file);
	
	const char* weight_file = argv[2]; // such as attenuation
	int number_of_radial_bins = scanner.getNumberOfRadialBins(); //atoi(argv[3]);
	if (number_of_radial_bins < 0) {
		SystemLog::write("must specify number of radial bins!\n");
		abort();
	}
	int number_of_angles = scanner.getSystemParameters().number_of_detector_modules.t * 
						   scanner.getSystemParameters().crystal_array_size.t / 2;
	int nring = scanner.getNumberOfCrystalRings();
	int number_of_planes = nring * nring;
	SystemLog::write("weight (sinogram) size: [%d(rad), %d(plane), %d(angle)]\n",
		number_of_radial_bins, number_of_planes, number_of_angles);
	
	int length = ListModeProjector::getFileLength(weight_file) / sizeof(float);
	if ((length / number_of_planes / number_of_angles) != number_of_radial_bins) {
		SystemLog::write("warning, the number of radial bins detected [%d] does "
						 "not match the expected value [%d].\n", 
						 length / number_of_planes / number_of_angles, number_of_radial_bins);
		//abort();
		SystemLog::write("either the file does not exist or the specified number of radial bins is incorrect!\n");
	}
	
	// read weight
#if 0
	Image<float> prj_weight(number_of_radial_bins, number_of_planes, number_of_angles);
	if (!prj_weight.read(weight_file)) {
		SystemLog::write("read weight failed. will use default weight.\n");
		prj_weight.set(1.0);
	}
#else // modified, so don't need permute outside
	SystemLog::write("assuming data is stored in terms of [#rad x #angle x #plane]\n");
	Image<float>* prj_tmp = new Image<float>(number_of_radial_bins, number_of_angles, number_of_planes);
	if (!prj_tmp->read(weight_file)) {
		SystemLog::write("read weight failed. will use default weight.\n");
	}
	SystemLog::write("permuting (also taking inverse) ...");
	Image<float> prj_weight(number_of_radial_bins, number_of_planes, number_of_angles);
	for (int k = 0; k < number_of_planes; k ++) {
		for (int j = 0; j < number_of_angles; j ++) {
			for (int i = 0; i < number_of_radial_bins; i ++) {
				prj_weight(i,k,j) = 1.0 / (*prj_tmp)(i,j,k);
			}
		}
	}
	SystemLog::write("done!\n");
	delete prj_tmp;
#endif

	// Toshiba's sinogram plane ordering
	std::vector<std::pair<int, int> > rp;
	for (int n = 0; n < nring; n ++) {
	
		SystemLog::write("ring difference = %d\n", n);

		if (n==0) {
			for (int m = 0; m < nring; m ++) {
				std::pair<int, int> p(m, m);
				rp.push_back(p);
				//
				SystemLog::write("[%d %d],", rp.back().first, rp.back().second);
				//
			}
		} else {
						
			for (int m = 0; m < nring - n; m ++) {
				std::pair<int, int> p(m, n + m);
				rp.push_back(p);
				//
				SystemLog::write("[%d %d],", rp.back().first, rp.back().second);
				//
			}
			for (int m = 0; m < nring - n; m ++) {
				std::pair<int, int> p(n + m, m);
				rp.push_back(p);
				//
				SystemLog::write("[%d %d],", rp.back().first, rp.back().second);
				//
			}
		}
		
		SystemLog::write("\n");
	
	}

	// create a listmode projector
	ListModeProjector* projector = new ListModeProjector();
	projector->initialize(scanner.getNumberOfCrystalRings(),
    					  scanner.getCrystalRingPitchSize(),
    					  scanner.getSystemParameters().number_of_detector_modules.t,
    					  scanner.getSystemParameters().crystal_array_size.t,
    					  scanner.getSystemParameters().detector_ring_diameter,
    					  scanner.getCrystalTransaxialPitchSize(),
    					  scanner.getSystemParameters().crystal_size.depth,
    					  scanner.getSystemParameters().image_size.i,
    					  scanner.getSystemParameters().image_size.k,
    					  scanner.getSystemParameters().voxel_size_i,
    					  scanner.getSystemParameters().voxel_size_k,
    					  -1, -1, 0);

	// create default xtal pairs (ids start from 0)
	Projector::XtalPair xp = projector->createCrystalPairs(
				scanner.getSystemParameters().number_of_detector_modules.t,
				scanner.getSystemParameters().crystal_array_size.t,
				number_of_radial_bins, 0, 0);

	//
	int seg = number_of_radial_bins * number_of_planes * 64;
	int nseg = number_of_angles / 64;
	std::vector<ListModeProjector::LMEVENT> lmdata(seg);
	Image<float> prjs(seg);
	Image<float> temp(scanner.getSystemParameters().image_size);
	Image<float> sensitivity(scanner.getSystemParameters().image_size);
	
	for (int m = 0; m < nseg; m ++) {

		SystemLog::write("processing seg #%d (%d in total [%d]) ...\n",
			m+1, nseg, seg);

		for (int n = 0; n < seg; n ++) {

			int gid = n + m * seg;
			int phi = gid / (number_of_radial_bins * number_of_planes);
			int r0 = gid % (number_of_radial_bins * number_of_planes);
			int nr = r0 / number_of_radial_bins;
			int rad = r0 % number_of_radial_bins;
			int nx = rad + phi * number_of_radial_bins;
		
			ListModeProjector::LMEVENT e;
			e.xtal_id_1 = xp[nx].first;
			e.xtal_id_2 = xp[nx].second;
			e.ring_id_1 = rp[nr].second; // switch
			e.ring_id_2 = rp[nr].first; 
			
			lmdata[n] = e;
		}
		
		memcpy(prjs.getPtr(), prj_weight.getPtr() + m * seg, sizeof(float) * seg);
		temp.clear();
		projector->doBackProj(temp, prjs, lmdata);
		
		for (int j = 0; j < sensitivity.getSize(); j ++) {
			sensitivity[j] += temp[j];
		}	
			
	}
	
	SystemLog::write("save to file: sensitivity.img.\n");
	sensitivity.write("sensitivity.img");
	
	SystemLog::close();

	return 0;
}
