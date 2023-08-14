#include <petsys_log.h>
using namespace UCDPETSYS;

bool SystemLog::working = false;
bool SystemLog::screenshot = true;
char SystemLog::charbuff[1024] = {0};
FILE* SystemLog::log_file = 0;

SystemLog::SystemLog()
{
}

SystemLog::~SystemLog()
{
}

void SystemLog::open(const char* filename)
{
	if (!SystemLog::working) { 		
		log_file = fopen(filename, "wt");
		
		if (log_file != 0) {
			
			time_t t;
			time(&t);
			std::fprintf(log_file, "Log file is created. %s\n\n", ctime(&t));
			
			SystemLog::working = true;
		}
	}
}

void SystemLog::close()
{
	if (SystemLog::working) {
		time_t t;
		time(&t);
		std::fprintf(log_file, "\n\nLog file is closed. %s\n", ctime(&t));
		fclose(log_file);
		SystemLog::working = false;
	}
}

void SystemLog::write(const char* format, ...)
{
	if (SystemLog::working) {
		va_list args;
		va_start (args, format);
  		vsprintf(SystemLog::charbuff, format, args);
  		va_end (args);
		std::fprintf(log_file, "%s", charbuff);  		  		
		if (SystemLog::screenshot) {
			std::fprintf(stderr, "%s", charbuff);
		}
	}
}
