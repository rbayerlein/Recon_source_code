#ifndef PETSYS_LOG_H
#define	PETSYS_LOG_H

#include <cstdio>
#include <cstdarg>
#include <ctime>

namespace UCDPETSYS 
{

class SystemLog 
{

public:
	SystemLog();
	~SystemLog();

public:
	static bool screenshot;
	static void write(const char* format, ...);
	static void open(const char* log_file);
	static void close();

private:
	static char charbuff[1024];
	static bool working;
	static FILE* log_file;
};


} // end of UCDPETSYS

#endif // petsys_log.h
