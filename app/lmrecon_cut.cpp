#include <petsys_cscanner_cut.h>

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
	SystemLog::write("initialized number of threads: %d (fp) %d (bp) %d (reg)\n", Projector::NUMBER_OF_THREADS_FP,
		Projector::NUMBER_OF_THREADS_BP, Regularizer::NUMBER_OF_THREADS);
#endif
	//

	// scanner.runRecon();
	
	scanner.runRecon2_cut();


	std::fprintf(stderr, "end of recon!\n");
	
	SystemLog::close();

	return 0;
}
