/*
 * grid.h
 *
 */
#include <mpi.h>

#include "mytypes.h"

#ifndef FALSE
#define FALSE (0)
#endif

#ifndef TRUE
#define TRUE (1)
#endif

/************************************************************************
 * macros
 ************************************************************************/
#ifdef USE_DOUBLE
#define MPI_SINGLE MPI_DOUBLE
#else
#define MPI_SINGLE MPI_FLOAT
#endif

/************************************************************************
 * parallel structure for communications
 ************************************************************************/
typedef struct {
    int      nproc;              /* number of processors == P * Q * R */
    int      P, Q, R;            /* global size of 2-d grid */
    int      myP, myQ, myR;      /* local grid position */
    int      me;                 /* local rank number */
    int      nbr[6];             /* nbr nodes (0=W, 1=E, 2=S, 3=N, 4=D, 5=U) */
    int      znbr[2];            /* neighnors in z */
    int      ynbr[2];            /* neighbor in y */
    int      xnbr[2];            /* neighbor in x */
    int      diag;               /* diagonal in xycom */
    MPI_Comm gridcom;            /* communicator for grid */
    MPI_Comm xcom;               /* communicator for x-direction */
    MPI_Comm ycom;               /* communicator for y-direction */
    MPI_Comm zcom;               /* communicator for z-direction */
    MPI_Comm xycom;              /* communicator for xy-direction */
    MPI_Datatype MPI_STATE;      /* Data typemap for state vector */  
} grid3d;


/************************************************************************
 * manifest constants
 ************************************************************************/

/* cardinal directions */
#define WEST  (0)
#define EAST  (1)
#define SOUTH (2)
#define NORTH (3)
#define DOWN  (4)
#define UP    (5)

#define DIAG  (6)

/* neighbor indexes */
#define BEFORE (0)
#define AFTER  (1)

/* communicator ordering */
#ifdef ORIG_ORDER
#define MP_X 0
#define MP_Y 1
#define MP_Z 2
#else
#define MP_X 2
#define MP_Y 1
#define MP_Z 0
#endif

/* message tags */
#define LIGHT_TAG 68
#define SBS_TAG 69
#define SHIFT_TAG 70
#define INFO_TAG 76

/************************************************************************
 * global variables
 ************************************************************************/

#ifdef __parm_init__

/* parallel parameters */
int     master = 0;              /* master node rank */
grid3d  g3;

int     dfile=FALSE;             /* write output debug messages to files or stdout */

#else

/* parallel parameters */
extern int     master;           /* master node rank */
extern grid3d  g3;

extern int     dfile;            /* write output debug messages to files or stdout */

#endif


/************************************************************************
 * functions 
 ************************************************************************/

/* define a 3-D Cartesian topology atop the hardware */
void build_grid(void);

/* simulation is finished, exit the message passing. */
void simdone(void);

/* initialize parallel simulation */
void siminit(int *ac, char **av[]);

/* reduce scalars from each node into a scalar result on the master node */
void simreduce(double in, double *out);

/* synchronize entire grid */
void simsync(void);
