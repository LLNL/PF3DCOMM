/* mytypes.h */
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#ifndef __MYTYPES_H__

#ifdef USE_DOUBLE
typedef double real;
#else
typedef float real;
#endif

typedef struct {
    real re;
    real im;
} complex;

/* stencil and downwind points */
#define NSTENCIL 6
#define NDOWNWIND 2

#define NCEQ   6                /* Need to fix this */
#define NMATMAX 4
#define MAXNEQ (NCEQ+NMATMAX-1)

/*******************************************************************************
 *                                                                mesh structure
 *******************************************************************************/
#if 0
typedef struct {
    real  *p;                    /* ion pressure */
    int    nxa, nya, nza;        /* allocated box dimensions */
} mesh;
#endif

#endif

#ifndef PI
#define PI        3.1415926535897932384626433832795029
#endif

#ifndef FALSE
#define FALSE     0
#define TRUE      1
#endif

#define min(a,b)    ( ( (a) < (b) ) ? (a) : (b) )
#define max(a,b)    ( ( (a) < (b) ) ? (b) : (a) )
#define sqr(a)      ( (a) * (a) )

#define __MYTYPES_H__
