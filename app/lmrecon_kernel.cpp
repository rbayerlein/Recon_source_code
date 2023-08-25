#include <petsys_cscanner.h>

#include <petsys_log.h>

int main(int argc, char* argv[])
{
	SystemLog::open("scanner.log");

	CScanner scanner;
	scanner.initialize(argv[1]);

	// paralle computing setup
#if USE_OMP
	if (argc > 2) {
		Projector::NUMBER_OF_THREADS_FP = atoi(argv[2]);		
	}
	if (argc > 3) {
		Projector::NUMBER_OF_THREADS_BP = atoi(argv[3]);		
	}
	if (argc > 4) {		
		Regularizer::NUMBER_OF_THREADS = atoi(argv[4]);
	}
	printf("initialized number of threads: %d (fp) %d (bp) %d (reg)\n", Projector::NUMBER_OF_THREADS_FP,
		Projector::NUMBER_OF_THREADS_BP, Regularizer::NUMBER_OF_THREADS);
#endif
	//

	// scanner.runRecon();


 //    const char* argv5 = argv[5];
	// char kmat_fname[512];
	// sprintf(kmat_fname, "%s", argv5);
	
	// PSpMatrix* kmat;	
	// std::size_t total_sz = 0;					
	// std::stringstream ostr;
	// ostr.str("");      // clear
	// ostr << kmat_fname;
	// kmat = new PSpMatrix();
	// kmat->read(ostr.str().c_str());	
	// total_sz += kmat->getNNZ()*sizeof(SSpMatrix::PACKAGE) +
	// 			(kmat->getRowNum()+1)*sizeof(int);

	// printf("OK! total memory space occupied: %lu (bytes)\n", total_sz);
	// printf("kmat->getNNZ(): %lu (bytes)\n", kmat->getNNZ());
	// printf("sizeof(SSpMatrix::PACKAGE): %lu (bytes)\n", sizeof(SSpMatrix::PACKAGE));
	// printf("(kmat->getRowNum()+1): %lu (bytes)\n", (kmat->getRowNum()+1));
	// printf("(kmat->getRowNum()+1)*sizeof(int): %lu (bytes)\n", (kmat->getRowNum()+1)*sizeof(int));



	scanner.runRecon2();
	// scanner.runRecon3();



	std::fprintf(stderr, "end of recon!\n");
	
	SystemLog::close();

	return 0;
}
