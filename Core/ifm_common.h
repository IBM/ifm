/* common defintions for ifm */

#ifndef _IFM_COMMON_H
#define _IFM_COMMON_H
#include <stdint.h>

/* This should change according to the target architecture */
#ifdef USE_FLOATS
 typedef float        real_t;
 #define ncReal       ncFloat
 #define real_sqrt    sqrtf
 #define real_fabs    fabsf
 #define real_pow     powf
 #define real_cbrt    cbrtf
 #define REAL_EPSILON FLT_EPSILON
#else
 typedef double       real_t;
 #define ncReal       ncDouble
 #define real_sqrt    sqrt
 #define real_fabs    fabs
 #define real_pow     pow
 #define real_cbrt    cbrt
 #define REAL_EPSILON DBL_EPSILON
#endif

/* just in case somebody wants to run this on a 32bit machine, need more testing though */
#if BITS == 64
#define MLONG long
#else
#define MLONG long long
#endif

/* macros for netCDF */
#define ERRCODE 2
#define NETCDF_ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

/* macro for testing float */
#define FLOAT_EQ(A, B) (fabs((A) - (B)) < 1e-5)
#define MIN2(A,B)  (A)<(B)?(A):(B)

#endif

/* end */
