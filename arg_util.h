#ifndef ARG_UTIL_H
#define ARG_UTIL_H

#include <stdlib.h>

int argtoi(const char *arg, int min, int max)
{
	int result = atoi(arg);
	if (result<min)
		result = min;
	
	if (result>max)
		result = max;

	return result;
}

float argtof(const char *arg, float min, float max)
{
	float result = (float)atof(arg);
	if (result<min)
		result = min;
	
	if (result>max)
		result = max;

	return result;
}

#pragma warning(disable:4800)
bool argtob(const char *arg, bool def)
{
	int i = atoi(arg);
	bool result = def;

	if (i==0 || i==1)
		result = (bool)i;

	return result;
}

#endif