#include <petsys_cscanner.h>

#include <petsys_log.h>

int main(int argc, char* argv[])
{
	if(argc < 2){
		std::cout << "not enough input arguments. Usage:\n" << argv[0] 
		<< " [cfg-file](required) [NUMBER_OF_THREADS_FP](opt) [NUMBER_OF_THREADS_BP](opt) [NUMBER_OF_THREADS_REG](opt) " << std::endl;
		return 0;
	}
	std::cout << "config file name: " << argv[1] << std::endl;
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
	SystemLog::write("initialized number of threads: %d (fp) %d (bp) %d (reg)\n", 
		Projector::NUMBER_OF_THREADS_FP,
		Projector::NUMBER_OF_THREADS_BP, 
		Regularizer::NUMBER_OF_THREADS
	);
#endif
	//

	// scanner.runRecon();	
	// scanner.runRecon2();	
	scanner.runexplmfp();


	std::fprintf(stderr, "end of recon!\n");
	
	SystemLog::close();

	return 0;
}
