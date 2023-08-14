#include <petsys_cscanner.h>

#include <petsys_log.h>

#include <cstdio>
#include <vector>
#include <fstream>
using namespace std;

typedef struct tagPRE_RECON_PARA {
    short sign_group;
    short symm_ind;
    short shift;
    short counts;
    int glbidx_p;
    float multi_fac;
    float add_fac; //<-not used!
} PRE_RECON_PARA;

typedef struct tagTOFBININFO {
    float dt;
    float factor;
} TOF_BIN_INFO;

int main(int argc, char* argv[])
{
	SystemLog::open("decdat.log");

	CScanner scanner;
	const char* cfg_file = argv[1];
	scanner.initialize(cfg_file);
	
	const char* prl_file = argv[2];

	int nring = scanner.getNumberOfCrystalRings();
	int num_of_rps = nring * nring;
	int num_of_phis = scanner.getSystemParameters().number_of_detector_modules.t * 
					  scanner.getSystemParameters().crystal_array_size.t / 2;
	int num_of_rads = scanner.getNumberOfRadialBins();
	if (num_of_rads < 0) {
		SystemLog::write("must specify number of radial bins!\n");
		abort();
	}
	
    float tw_spacing = scanner.getSystemParameters().tof_info.bin_size;
	SystemLog::write("%d,%d,%d,%d,%f\n", nring, num_of_rps, num_of_phis, num_of_rads, tw_spacing);

	// create default xtal pairs (ids start from 0)
	ListModeProjector* projector = new ListModeProjector();
	Projector::XtalPair xtal_pair = projector->createCrystalPairs(
				scanner.getSystemParameters().number_of_detector_modules.t,
				scanner.getSystemParameters().crystal_array_size.t, 
				num_of_rads, 0, 0);
    
//    for (Projector::XtalPair::iterator it = xtal_pair.begin(); it != xtal_pair.end(); it ++) {
//    	printf("%d %d\n", (*it).first, (*it).second);
//    }

    // make ring idx lookup tables
    Image<int> ring_ind(nring, nring);
    Image<std::pair<int, int> > ring_pair(nring*nring); 
    int offset = 0;    
    for (int i = 0; i < nring; i ++) {
        for (int  j = 0; j < (nring-i); j ++) {
            if (i == 0) { // direct plane
                ring_ind(j, j) = j;
                ring_pair[j].first = j;
                ring_pair[j].second = j;
            } else { // oblique planes
                ring_ind(j,j+i) = offset + j;
                ring_pair[offset+j].first = j;
                ring_pair[offset+j].second = j+i;
                
                ring_ind(j+i,j) = offset + nring-i+j;
                ring_pair[offset+nring-i+j].first = j+i;
                ring_pair[offset+nring-i+j].second = j;                
            }
        }
        offset += (i==0) ? nring : 2*(nring-i);
    } 
//    ring_ind.write("ro.map");
//    ring_pair.write("rp.map");
    
    // open lm file
    FILE* fp = fopen(prl_file, "rb");
    if (fp == 0) {
        SystemLog::write("open file failed!\n");
        abort();
    }

#if USE_NONTOF
	SystemLog::write("nontof not implemented yet!\n");
	abort();
#endif

	char str[512];
#if USE_NONTOF
//    FILE* mylm = fopen("my.nontof.lm", "wb");
#else    
//	sprintf(str, "%s.lm", prl_file);
    FILE* mylm = fopen("my.lm", "wb");
//    sprintf(str, "%s.add_fac", prl_file);
    FILE* mylm_add = fopen("my.add_fac", "wb");
//    sprintf(str, "%s.mul_fac", prl_file);
    FILE* mylm_mul = fopen("my.mul_fac", "wb");
    FILE* mylm_gid = fopen("my.gid", "wb");
    FILE* mylm_aid = fopen("my.aid", "wb");
#endif    
    if (mylm == 0 || mylm_add == 0 || mylm_mul == 0 || mylm_gid == 0 || mylm_aid == 0) {
    	SystemLog::write("open file failed!\n");
    	abort();
    }
    
    // PRE-RECON List-mode Head
    // check HeadTable
    int lors, events;
    int total_num_of_events = 0;
    for (int i = 0; i < num_of_phis; i ++) {
        fread(&lors, sizeof(int), 1, fp);
        fread(&events, sizeof(int), 1, fp);
//        if ((i+1) % 8 == 0) {
//            SystemLog::write("#%d: [%d (LORs), %d (events)]", i+1, lors, events); 
//        }
        total_num_of_events += events;
    }
	SystemLog::write("total: %d\n", total_num_of_events);

	float mul_fac, add_fac;
    
    // PRE-RECON BOXES
    PRE_RECON_PARA pbox;
    TOF_BIN_INFO tbin;
	ListModeProjector::LMEVENT e;
	
    int total_lm = 0;
    int valid_lm = 0;
    int cc = 0;
    
    while (total_lm < total_num_of_events) {//(!feof(fp)) {
    
        fread(&pbox, sizeof(PRE_RECON_PARA), 1, fp);
#if USE_NONTOF
        
//        tracer("%d\n", pbox.counts);
/*        
        total_lm += pbox.counts;
                
        if (pbox.counts > 0) {

            int phi = pbox.glbidx_p / (num_of_rads*num_of_rps);                
            int rm = pbox.glbidx_p % (num_of_rads*num_of_rps);
            int nr = rm / num_of_rads;
            int rad = rm % num_of_rads;                
            int nx = rad + phi * num_of_rads;

            for (int nn = 0; nn < pbox.counts; nn++) {
                e.xtal_id_1 = xtal_pair[nx].first - 1; // because of 1based
                e.xtal_id_2 = xtal_pair[nx].second - 1;                
                e.ring_id_2 = ring_pair[nr].first; ///!!!!
                e.ring_id_1 = ring_pair[nr].second;                
                e.tbin_id = 0; // not used
                e.mul_fac = pbox.multi_fac;
                e.add_fac = pbox.add_fac;
                ebuff[nn] = e;
                
                //
                ee[nn].xtal_id_1 = e.xtal_id_1;
                ee[nn].ring_id_1 = e.ring_id_1;
                ee[nn].xtal_id_2 = e.xtal_id_2;
                ee[nn].ring_id_2 = e.ring_id_2;
                ee[nn].tbin_id = 0;
            }
            
            fwrite(&ebuff[0], sizeof(MYLMEVENT)*pbox.counts, 1, myfp);
//            fwrite(&e, sizeof(MYLMEVENT), 1, myfp);

            valid_lm ++;
            
//            ee.xtal_id_1 = e.xtal_id_1;
//            ee.xtal_id_2 = e.xtal_id_2;
//            ee.ring_id_1 = e.ring_id_1;
//            ee.ring_id_2 = e.ring_id_2;
//            ee.tbin_id = 0; // not used
            
//            fwrite(&ee, sizeof(ListModeProjector::LMEVENT), pbox.counts, mylm);
            fwrite(&ee[0], sizeof(ListModeProjector::LMEVENT)*pbox.counts, 1, mylm);
        
        }
*/
        
#else                
        if (pbox.counts > 0) {
        
            if (cc < 8) {
                SystemLog::write("global_id = %d, counts = %d, multi_fac = %g\n", 
                    pbox.glbidx_p, pbox.counts, pbox.multi_fac);
                    cc ++;
            }       
        
            for (int c = 0; c < pbox.counts; c ++) {            
            
                fread(&tbin, sizeof(TOF_BIN_INFO), 1, fp);

                int phi = pbox.glbidx_p / (num_of_rads*num_of_rps);                
                int rm = pbox.glbidx_p % (num_of_rads*num_of_rps);
                int nr = rm / num_of_rads;
                int rad = rm % num_of_rads;                
                int nx = rad + phi * num_of_rads;
                
                //
                e.xtal_id_1 = xtal_pair[nx].first;// - 1; // because of 1based
                e.xtal_id_2 = xtal_pair[nx].second;// - 1;                
// old:
                e.ring_id_2 = ring_pair[nr].first; ///!!!!
                e.ring_id_1 = ring_pair[nr].second;                
// test: TEST is wrong!!!
//                e.ring_id_1 = ring_pair[nr].first; ///!!!!
//                e.ring_id_2 = ring_pair[nr].second;                

// old: -
                e.tbin_id = short(-tbin.dt / tw_spacing); // quantization
// test: +                
//                e.tbin_id = short(+tbin.dt / 5.0); // quantization


                mul_fac = pbox.multi_fac;
                add_fac = tbin.factor;
                                
                fwrite(&e, sizeof(ListModeProjector::LMEVENT), 1, mylm);
                fwrite(&mul_fac, sizeof(float), 1, mylm_mul);
                fwrite(&add_fac, sizeof(float), 1, mylm_add);
                fwrite(&pbox.glbidx_p, sizeof(int), 1, mylm_gid);   
                fwrite(&phi, sizeof(int), 1, mylm_aid);             
                valid_lm ++;
                
                if (cc < 8) {
                    SystemLog::write("\tdt=%g,\tadd_fac=%g\n", tbin.dt, tbin.factor);
                }

                if (((total_lm+1) % 1000000 == 0)) {
                    SystemLog::write("processing #%d [XTAL:[%d %d], RING:[%d %d]), %d, %.20f, %.20f] ...\n", 
                          total_lm+1, e.xtal_id_1, e.xtal_id_2, e.ring_id_1, e.ring_id_2,
                           e.tbin_id, mul_fac, add_fac);
                }
                
                total_lm ++;                
            }       
        }

#endif        
    }
    
    SystemLog::write("# of events written = %d\n", valid_lm);
        
    fclose(mylm);
    fclose(mylm_add);
    fclose(mylm_mul);
    fclose(mylm_gid);
    fclose(mylm_aid);

    return 0;
}
