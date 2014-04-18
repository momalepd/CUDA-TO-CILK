#include <stdio.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <list>

#include <time.h>
#include <fstream>
#include <string>
#include <iostream>

#include <unistd.h>
#include <cassert>
#include <inttypes.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>


#define MYINFINITY	100000000
#define SWAP(x, y)	{tmp = x; x = y; y = tmp;}

typedef unsigned foru;

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}
