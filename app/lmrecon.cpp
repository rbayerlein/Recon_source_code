#include <petsys_cscanner.h>

#include <petsys_log.h>

int main(int argc, char* argv[])
{
	SystemLog::open("scanner.log");

	CScanner scanner;
	scanner.initialize(argv[1]);

	// paralle computing setup
	// pre initialize number of threads, rbayerlein, 2021-11-15
	int num_threads = 48;
#if USE_OMP
	Projector::NUMBER_OF_THREADS_FP = num_threads;
	if (argc > 2) {
		Projector::NUMBER_OF_THREADS_FP = atoi(argv[2]);		
	}
	Projector::NUMBER_OF_THREADS_BP = num_threads;
	if (argc > 3) {
		Projector::NUMBER_OF_THREADS_BP = atoi(argv[3]);		
	}
	Regularizer::NUMBER_OF_THREADS = num_threads;
	if (argc > 4) {		
		Regularizer::NUMBER_OF_THREADS = atoi(argv[4]);
	}
	SystemLog::write("initialized number of threads: %d (fp) %d (bp) %d (reg)\n", Projector::NUMBER_OF_THREADS_FP,
		Projector::NUMBER_OF_THREADS_BP, Regularizer::NUMBER_OF_THREADS);
#endif
	//

	// scanner.runRecon();
	
	scanner.runRecon2();


	std::fprintf(stderr, "end of recon! Finish\n");
	
	SystemLog::close();

	return 0;
}
