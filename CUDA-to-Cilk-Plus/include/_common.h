#ifndef LSG_COMMON
#define LSG_COMMON

#include <stdio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <iostream>
#include <limits>
#include <string.h>

#include <unistd.h>
#include <cassert>
#include <inttypes.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <inttypes.h>

// For MAC and FreeBSD: by Rashid Kaleem.
#ifdef __APPLE__ 
#include <libkern/OSByteOrder.h>
#  define le64toh(x) OSSwapLittleToHostInt64(x)
#  define le32toh(x) OSSwapLittleToHostInt32(x)
#elif __FreeBSD__ 
#  include <sys/endian.h>
#elif __linux__ 
#  include <endian.h>
#  ifndef le64toh
#    if __BYTE_ORDER == __LITTLE_ENDIAN
#      define le64toh(x) (x)
#      define le32toh(x) (x)
#    else
#      define le64toh(x) __bswap_64 (x)
#    endif
#  endif
#endif

#ifndef LSGDEBUG
#define LSGDEBUG 0
#endif 

#define dprintf	if (debug) printf
unsigned const debug = LSGDEBUG;

typedef unsigned foru;
//typedef float foru;

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

#endif
