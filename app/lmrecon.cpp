#include <petsys_cscanner.h>

#include <petsys_log.h>
/**
 * @brief Main entry point for the PET system scanner.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return 0 on successful execution.
 */
int main(int argc, char* argv[])
{
	SystemLog::open("scanner.log"); ///< create a log file

	CScanner scanner;
	scanner.initialize(argv[1]); ///< create an initialize an instande of CScanner

	// paralle computing setup
	// pre initialize number of threads, rbayerlein, 2021-11-15
	int num_threads = 48;
#if USE_OMP
	// Code related to parallelization using OpenMP
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

	// Running an updated version of the reconstruction algorithm
	scanner.runRecon2(); ///< Run recon 2. Same as runRecon(), but there is option to invoke KEM reconstruction 


	std::fprintf(stderr, "end of recon! Finish\n");
	
	SystemLog::close(); ///< close log file

	return 0; ///< return 0. exit program.
}
