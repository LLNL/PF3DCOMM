/*
 * simMPI.c
 *
 * this file provides the basic communications routines for the pF3d/Z3
 * codes using MPI.
 *
 */

#include "grid.h"
#include "main.h"

int nodes, mp_p, mp_q, mp_r, mp_myp, mp_myq, mp_myr, mp_rank;

/*------------------------------------------------------------------------
 * build_grid
 *
 * define a 3-D Cartesian topology atop the hardware
 */

#define REORDER FALSE

void build_grid(void)
{
    int       dims[3] = {0, 0, 0};              /* global size of 3-d grid */
    int       coord[3] = {0, 0, 0};             /* local grid position */
    int       periods[3] = {TRUE, TRUE, TRUE};  /* periodic in each dimension */
    int       temp[3];                          /* for allocating x/y/z comms */

    MPI_Comm_rank(MPI_COMM_WORLD, &g3.me);
    MPI_Comm_size(MPI_COMM_WORLD, &g3.nproc);

#if defined(DEBUG) && (DEBUG>128)
    /* print this out for debugging purposes */
    (void) printf("bld[%d]: Process %d, size %d\n", g3.me, g3.me, g3.nproc);
#endif
    
    dims[MP_X] = mp_p;
    dims[MP_Y] = mp_q;
    dims[MP_Z] = mp_r;
    MPI_Dims_create(g3.nproc, 3, dims);
    g3.P = dims[MP_X];
    g3.Q = dims[MP_Y];
    g3.R = dims[MP_Z];
    
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, REORDER, &g3.gridcom);
    MPI_Cart_get(g3.gridcom, 3, dims, periods, coord);
    g3.myP = coord[MP_X];
    g3.myQ = coord[MP_Y];
    g3.myR = coord[MP_Z];
    
#if defined(DEBUG) && (DEBUG>128)
    if (g3.me == master)
        (void) printf("bld[0]: size of mesh is %d by %d by %d\n", g3.P, g3.Q, g3.R);
    (void) printf("bld[%d]: coord is (%d,%d,%d)\n", g3.me, g3.myP, g3.myQ, g3.myR);
#endif
    
    /* get six neighbors (0=W, 1=E, 2=S, 3=N, 4=D, 5=U) */
    /* guard against negative node numbers because IBM MPI doesn't handle that */
    coord[MP_X] = ( g3.myP - 1 + g3.P ) % g3.P;
    coord[MP_Y] = g3.myQ;
    coord[MP_Z] = g3.myR;
    MPI_Cart_rank(g3.gridcom, coord, &g3.nbr[0]);

    coord[MP_X] = ( g3.myP + 1 ) % g3.P;
    coord[MP_Y] = g3.myQ;
    coord[MP_Z] = g3.myR;
    MPI_Cart_rank(g3.gridcom, coord, &g3.nbr[1]);

    coord[MP_X] = g3.myP;
    coord[MP_Y] = ( g3.myQ - 1 + g3.Q ) % g3.Q;
    coord[MP_Z] = g3.myR;
    MPI_Cart_rank(g3.gridcom, coord, &g3.nbr[2]);

    coord[MP_X] = g3.myP;
    coord[MP_Y] = ( g3.myQ + 1 ) % g3.Q;
    coord[MP_Z] = g3.myR;
    MPI_Cart_rank(g3.gridcom, coord, &g3.nbr[3]);
    
    coord[MP_X] = g3.myP;
    coord[MP_Y] = g3.myQ;
    coord[MP_Z] = ( g3.myR - 1 + g3.R ) % g3.R;
    MPI_Cart_rank(g3.gridcom, coord, &g3.nbr[4]);
    
    coord[MP_X] = g3.myP;
    coord[MP_Y] = g3.myQ;
    coord[MP_Z] = ( g3.myR + 1 ) % g3.R;
    MPI_Cart_rank(g3.gridcom, coord, &g3.nbr[5]);
    
#if defined(DEBUG) && (DEBUG>128)
    (void) printf("bld[%d]: nbrs = %d/%d/%d/%d/%d/%d (EWNSUD)\n", g3.me,
                  g3.nbr[1], g3.nbr[0], g3.nbr[3], g3.nbr[2], g3.nbr[5], g3.nbr[4]);
#endif

    /* define z communicator */
    temp[MP_X] = FALSE; temp[MP_Y] = FALSE; temp[MP_Z] = TRUE;
    MPI_Cart_sub(g3.gridcom, temp, &g3.zcom);

    /* get z neighbors (DOWN and UP) */
    coord[0] = ( g3.myR - 1 + g3.R ) % g3.R;
    MPI_Cart_rank(g3.zcom, coord, &g3.znbr[BEFORE]);
    coord[0] = ( g3.myR + 1 ) % g3.R;
    MPI_Cart_rank(g3.zcom, coord, &g3.znbr[AFTER]);

    /* define y communicator */
    temp[MP_X] = FALSE; temp[MP_Y] = TRUE; temp[MP_Z] = FALSE;
    MPI_Cart_sub(g3.gridcom, temp, &g3.ycom);

    /* get y neighbors (SOUTH and NORTH) */
    coord[0] = ( g3.myQ - 1 + g3.Q ) % g3.Q;
    MPI_Cart_rank(g3.ycom, coord, &g3.ynbr[BEFORE]);
    coord[0] = ( g3.myQ + 1 ) % g3.Q;
    MPI_Cart_rank(g3.ycom, coord, &g3.ynbr[AFTER]);

    /* define x communicator */
    temp[MP_X] = TRUE; temp[MP_Y] = FALSE; temp[MP_Z] = FALSE;
    MPI_Cart_sub(g3.gridcom, temp, &g3.xcom);

    /* get y neighbors (WEST and EAST) */
    coord[0] = ( g3.myP - 1 + g3.P ) % g3.P;
    MPI_Cart_rank(g3.xcom, coord, &g3.xnbr[BEFORE]);
    coord[0] = ( g3.myP + 1 ) % g3.P;
    MPI_Cart_rank(g3.xcom, coord, &g3.xnbr[AFTER]);

    /* define xy-communicator for sub-global collective operations */
    temp[MP_X] = TRUE;temp[MP_Y] = TRUE; temp[MP_Z] = FALSE;
    MPI_Cart_sub(g3.gridcom, temp, &g3.xycom);

    /* get diagonal xy neighbor */
#ifdef ORIG_ORDER
    coord[0] = (g3.P -1 - g3.myP);
    coord[1] = (g3.Q -1 - g3.myQ);
#else
    coord[1] = (g3.P -1 - g3.myP);
    coord[0] = (g3.Q -1 - g3.myQ);
#endif
    MPI_Cart_rank(g3.xycom, coord, &g3.diag);

    /* define the MPI_STATE typemap for mfhydro communications */
    MPI_Type_contiguous(MAXNEQ, MPI_SINGLE, &g3.MPI_STATE);
    MPI_Type_commit(&g3.MPI_STATE);

    /* export some info */
    mp_p   = g3.P;   mp_q   = g3.Q;   mp_r   = g3.R;
    mp_myp = g3.myP; mp_myq = g3.myQ; mp_myr = g3.myR;
    mp_rank= g3.me;
}


/*------------------------------------------------------------------------
 * simdone
 *
 * simulation is finished, exit the message passing.
 */

void simdone(void)
{
    MPI_Finalize();
}


/*------------------------------------------------------------------------
 * siminit
 *
 * initialize parallel simulation 
 */

void siminit(int *ac, char **av[])
{
    /* initialize; get rank and size information */
    MPI_Init(ac,av);
    MPI_Comm_rank(MPI_COMM_WORLD, &g3.me);
    MPI_Comm_size(MPI_COMM_WORLD, &g3.nproc);
}


/*------------------------------------------------------------------------
 * simsync
 *
 * synchronize entire grid
 */
void simsync(void)
{
    double tim0, tim1;

    MPI_Barrier(g3.gridcom);
}
