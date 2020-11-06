/*
 * main.c
 *
 * driver for testing g3MPI
 *
 */
#define __parm_init__

#include <stdio.h>
#include "grid.h"
#include "main.h"
#include <sys/time.h>
#include <string.h>

#define NSAV 20

#define MAXTYPLEN 100
char RunType[MAXTYPLEN];

real *arry, *brry, *zrry;
long nxtot, nytot, nztot;
int  nxl, nyl, nzl;     /* size of one "tile" in the checkerboard decomp */
double fom;    /* bytes/sec/proc for most recently processed test */
FILE *file_csv;
char RUN_TYPE[100];
int do_csv= 0;

extern void print_one(int nsav, double *time_one, double *bytes_one,
                      char *label, int dir);
extern void print_full_fft(int nsav, double *time_one, double *bytes_one,
                           double *time_two, double *bytes_two,
                           char *label, int dir);
extern void getdims(int dims[5]);
extern int main(int argc, char **argv);
extern void print_time(char *msg);

void print_time(char *msg)
{
    struct timeval tv;
    int n;
    time_t      sec;     /* seconds */
    suseconds_t usec;    /* microseconds */
    double time;

    if(g3.me == master) {
      n= gettimeofday(&tv, 0);
      sec= tv.tv_sec;
      usec= tv.tv_usec;
      time= sec+1.0e-6*usec; /* currently not used */
      printf("%s %ld.%06ld\n", msg, sec, usec);
    }
}

int main(int argc, char **argv)
{
    int      nBuf=65536;            /* default buffer size */
    int      nTrip;                 /* numer of tripstimes to run each test */
    int      i, ip, j;              /* counter */
    int      procPerNode;
    double   startTime, stopTime, deltaT;  /* timing variables */
    double   sum, sumsq;            /* statistics variables */
    double   avgTime, errTime;      /* mean and standard error of time sample */
    double fom_2dfft, fom_fwd;
    int numproc, nmax, ival;
    long num_xyz;
    double time_one[NSAV], bytes_one[NSAV];
    double time_two[NSAV], bytes_two[NSAV];
    double *time_fulldata, *time_fulldata2;
    double time_lo[NSAV], time_hi[NSAV], time_sum[NSAV], time_sum2[NSAV];
    double time_avg[NSAV], time_std[NSAV];
    double timbase, time_begin, time_finish, time_total;

    int run_2Dfft= 1;
    int run_rowall= 0;
    int run_colall= 0;
    int run_colshift= 0;
    int run_rowshift= 0;
    int run_plane= 1;
    int run_hydro= 1;
#define MAXNAM 200
    char CSVNAM[MAXNAM];
  
    /* ------------------------------ Initialize MPI */
    siminit( &argc, &argv );
    
    time_fulldata= (double *) malloc(sizeof(double)*NSAV*g3.nproc);
    if(!time_fulldata) {
      puts("Not enough room to gather timing data from all processes\n");
      exit(-1);
    }
    time_fulldata2= (double *) malloc(sizeof(double)*NSAV*g3.nproc);
    if(!time_fulldata2) {
      puts("Not enough room to gather timing data from all processes\n");
      exit(-1);
    }

    
    /* ------------------------------ create cartesian grid and coordinators */
    /* default decomp dimensions make sure the code will run,
       but the FFT results are of no interest with a 1x1 decomp. */
    mp_p = 1;
    mp_q = 1;

    /* read the decomposition from the command line, if present */
    if(argc > 1) {
      mp_p = atoi(argv[1]);
    }
    if(argc > 2) {
      mp_q = atoi(argv[2]);
    }
    if(g3.nproc < mp_p*mp_q) {
      printf("Decomposition (%d, %d, NNN) inconsistent with %d processees",
             mp_p, mp_q, g3.nproc);
      exit(1);
    }
    /* compute mp_q to use all processors */
    mp_r = g3.nproc/(mp_p*mp_q);

    if(g3.nproc != mp_p*mp_q*mp_r) {
      printf("Decomposition (%d, %d, %d) inconsistent with %d processees",
             mp_p, mp_q, mp_r, g3.nproc);
      exit(1);
    }
    
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    build_grid();

    /* This is a variant of Bert Still's g3test4, a message passing test
       designed to evaluate the performance of the parallel 2D FFTs
       used in pf3d.
       It also includes a test of message passing performance
       designed to mimic the pf3d hydro package and the advection
       of light waves in the z-direction.

       The problem may be decomposed in all three directions.

       The x-, y-, and z-extent of a domain are hard-coded.
    */
    nxl= 1024;
    nyl= 768;
    nzl= 20;
    if(argc > 3) {
      nxl = atoi(argv[3]);
    }
    if(argc > 4) {
      nyl = atoi(argv[4]);
    }
    if(argc > 5) {
      nzl = atoi(argv[5]);
    }
    nxtot= nxl*mp_p;
    nytot= nyl*mp_q;
    nztot= mp_r;
    /* ------------------- set the number of trips */
    if(argc > 6) {
      nTrip = atoi(argv[6]);
    } else {
      nTrip = 20;
    }
    if(nTrip < NSAV) nTrip= NSAV;
    if(argc > 7) {
      nodes = atoi(argv[7]);
    } else {
      nodes= 4;
    }
    if(argc > 8) {
      strncpy(RunType, argv[8], MAXTYPLEN);
    } else {
      strncpy(RunType, "BaseRun", MAXTYPLEN);
    }
    procPerNode= mp_p*mp_q*mp_r/nodes;
    
    /* compute the number of zones in a full array */
    num_xyz= nxl*nyl*nzl;
    /* create storage buffers large enough for all tests.
       The alltoall test uses complex nxl-by-nyl-by-nzl arrays,
       so use that. */
    /* storage_init(2*num_xyz); */
    /* The FFT functions use real variables as of June 8, 2020,
       so allocate those. */
    storage_init(num_xyz);


    sum= sumsq= 0.0;

    /* ------------------------------ output parameters */
    if ( g3.me == master ) {
	printf("Running on %d CPUs, nTrip = %d\nDecomposition is mp_p=%d, mp_q=%d, mp_r=%d\n",
	       g3.nproc, nTrip, g3.P, g3.Q, g3.R);
        printf("\n");
	printf("------------------------------------------------------------------------------\n");
	printf("Test performance for the message passsing patterns used in pF3D\n");
	printf("------------------------------------------------------------------------------\n");
        printf("Each domain is %d by %d by %d zones\n\n", nxl, nyl, nzl);
    }
#ifdef BGQ
    /* determine the topology of the partition on a Blue Gene/Q system */
    {
      int dims[5];

      /* NOTE that this calls MPI functions so it must be run by all
         processes */
      getdims(dims);
      if ( g3.me == master ) {
        printf("The partition is %d by %d by %d by %d by %d nodes\n\n",
               dims[0]+1, dims[1]+1, dims[2]+1, dims[3]+1, dims[4]+1);
      }
    } 
#endif

    for (i=0; i < num_xyz; i++) {
	arry[i] = (real) rand() / (real) RAND_MAX;
	brry[i] = 0.0;
    }

    
    /* Only MPI rank zero should write a csv file */
    if(mp_rank == 0) do_csv= 1;
    else do_csv= 0;
    if(do_csv) {
      snprintf(CSVNAM, MAXNAM, "results_%d_nodes_%s_%d_per_node.csv", nodes, RunType, procPerNode);
      file_csv= fopen(CSVNAM, "w");
      if(file_csv) {
        fputs("Run type, nodes, proc_per_node, mp_p, mp_q, mp_r, nxl, nyl, nzl, 2D-fft time, std dev, 2D-fft rate GB/s, Advect time, std dev, Advect rate GB/s, Hydro time, std dev, Hydro rate GB/s\n", file_csv);
        fprintf(file_csv, "%s, %d, %d, %d, %d, %d, %d, %d, %d, ", RunType, nodes, procPerNode, mp_p, mp_q, mp_r, nxl, nyl, nzl);
      } else {
        puts("Failed to create results.csv");
        exit(-1);
      }
    }
    time_begin= MPI_Wtime();

  /* ------------------------------ start the shift test with Sendrecv calls */
  if( run_colshift) {
    print_time("START time is");
    MPI_Barrier(g3.gridcom);
    for (i=0; i<nTrip; i++) {
	startTime = MPI_Wtime();

        torows(arry, brry, nxl, nyl);

        stopTime = MPI_Wtime();
        if(i < NSAV) {
          time_one[i]= stopTime - startTime;
          bytes_one[i]= bytes_moved;
        }
    }
    MPI_Barrier(g3.gridcom);
    print_time("STOP time is");
    if ( g3.me == master ) {
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\tpass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
      print_one(NSAV, time_one, bytes_one, "Convert to rows", 1);
    }
  }

  /* ------------------------------ start the shift test with Sendrecv calls */
  if( run_rowshift) {
    print_time("START time is");
    MPI_Barrier(g3.gridcom);
    for (i=0; i<NSAV; i++) {
	startTime = MPI_Wtime();

        fromrows(arry, brry, nxl, nyl);

        stopTime = MPI_Wtime();
        if(i < NSAV) {
          time_one[i]= stopTime - startTime;
          bytes_one[i]= bytes_moved;
        }
    }
    MPI_Barrier(g3.gridcom);
    print_time("STOP time is");
    if ( g3.me == master ) {
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\tpass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
      print_one(NSAV, time_one, bytes_one, "Convert to rows", 1);
    }
  }

  /* ------------------------------ start the shift test with Sendrecv calls */
  if( run_rowshift) {
    if ( g3.me == master ) {
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\tpass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
    }

    /* ------------------------------ compute the timing statistics */
    avgTime = sum / (double) NSAV;
    errTime = sqrt( (NSAV * sumsq - sqr(sum)) / (sqr(NSAV) * (NSAV - 1)) );

    /* ------------------------------ output the results */
    if ( g3.me == master ) {
        print_one(NSAV, time_one, bytes_one, "Convert from rows", 1);
    }
  }

  /* ------------------------------ start the row shift test with Alltoall calls */
  if( run_rowall) {
    print_time("START time is");
    MPI_Barrier(g3.gridcom);
    for (i=0; i<nTrip; i++) {
	startTime = MPI_Wtime();

        torowsall(arry, brry, nxl, nyl);

        stopTime = MPI_Wtime();
        if(i < NSAV) {
          deltaT= stopTime - startTime;
	  if ( g3.me == master ) printf("torowsall time= %e\n", deltaT);
          time_one[i]= deltaT;
          bytes_one[i]= bytes_moved;
        }
    }
    MPI_Barrier(g3.gridcom);
    print_time("STOP time is");
    if ( g3.me == master ) {
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\ttorowsall pass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
      print_one(NSAV, time_one, bytes_one, "Convert to rows (Alltoall)", 1);
    }
  }

  /* ------------------------------ start the column shift test with Alltoall calls */
  if( run_colall) {
    print_time("START time is");
    MPI_Barrier(g3.gridcom);
    for (i=0; i<nTrip; i++) {
	startTime = MPI_Wtime();

        tocolsall(arry, brry, nxl, nyl);

        stopTime = MPI_Wtime();
        if(i < NSAV) {
          time_one[i]= stopTime - startTime;
          bytes_one[i]= bytes_moved;
        }
    }
    MPI_Barrier(g3.gridcom);
    print_time("STOP time is");
    if ( g3.me == master ) {
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\ttocolsall pass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
      print_one(NSAV, time_one, bytes_one, "Convert to columns (Alltoall)", 1);
    }
  }


  /* ------------- start the row+column shift tests with Alltoall calls */
  if( run_2Dfft) {
    for (i=0; i<NSAV; i++) {
        time_one[i]= 0.0;
        bytes_one[i]= 0;
        time_two[i]= 0.0;
        bytes_two[i]= 0;
    }
    print_time("START time is");
    MPI_Barrier(g3.gridcom);
    for (i=0; i<nTrip; i++) {
	startTime = MPI_Wtime();

        torowsall(arry, brry, nxl, nyl);
        if(i < NSAV) bytes_one[i] += bytes_moved;
        fromrowsall(arry, brry, nxl, nyl);

        stopTime = MPI_Wtime();
        if(i < NSAV) {
          time_one[i] = stopTime - startTime;
          bytes_one[i] += bytes_moved;
        }

	startTime = MPI_Wtime();

        tocolsall(arry, brry, nxl, nyl);
        if(i < NSAV) bytes_two[i] += bytes_moved;
        fromcolsall(arry, brry, nxl, nyl);

        stopTime = MPI_Wtime();
        if(i < NSAV) {
          time_two[i] = stopTime - startTime;
          bytes_two[i] += bytes_moved;
        }
    }
    MPI_Barrier(g3.gridcom);
    print_time("STOP time is");
    /* Gather the timing data from all processes back to rank zero.
       Use an Allgather, even though only rank zero needs the data,
       to keep the code uniform.
    */
    ival= MPI_Allgather(time_one, NSAV, MPI_DOUBLE, time_fulldata, NSAV, MPI_DOUBLE, g3.gridcom);
    ival= MPI_Allgather(time_two, NSAV, MPI_DOUBLE, time_fulldata2, NSAV, MPI_DOUBLE, g3.gridcom);
    if ( g3.me == master ) {
      /* Compute min, max, avg, and standard deviation of time over
         all MPI processes. */
      for(i=0; i<NSAV; i++) {
        time_lo[i]= 1.0e40;
        time_hi[i]= -1.0e40;
        time_sum[i]= 0.0;
        time_sum2[i]= 0.0;
        for(j= 0; j < g3.nproc; j++) {
          timbase= time_fulldata[i*g3.nproc+j];
          time_lo[i]= min(time_lo[i], timbase);
          time_hi[i]= max(time_hi[i], timbase);
          time_sum[i] += timbase;
          time_sum2[i] += timbase*timbase;
        }
        time_avg[i]= time_sum[i]/g3.nproc;
        time_std[i]= sqrt((time_sum2[i]-g3.nproc*time_avg[i]*time_avg[i]) / (g3.nproc-1));
      }
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\tX-fft pass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
      puts(" ");
      for(i=0; i<NSAV; i++) {
          printf("\tX-fft pass %2d avg %.3f sec min %.3f, max %.3f, std dev %.3f\n", i,
                 time_avg[i], time_lo[i], time_hi[i], time_std[i]);
      }
      puts(" ");
      for(i=0; i<NSAV; i++) {
        time_lo[i]= 1.0e40;
        time_hi[i]= -1.0e40;
        time_sum[i]= 0.0;
        time_sum2[i]= 0.0;
        for(j= 0; j < g3.nproc; j++) {
          timbase= time_fulldata2[i*g3.nproc+j];
          time_lo[i]= min(time_lo[i], timbase);
          time_hi[i]= max(time_hi[i], timbase);
          time_sum[i] += timbase;
          time_sum2[i] += timbase*timbase;
        }
        time_avg[i]= time_sum[i]/g3.nproc;
        time_std[i]= sqrt((time_sum2[i]-g3.nproc*time_avg[i]*time_avg[i]) / (g3.nproc-1));
      }
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\tY-fft pass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_two[i], 1.0e-9 * bytes_two[i] / time_two[i]);
      }
      puts(" ");
      for(i=0; i<NSAV; i++) {
          printf("\tY-fft pass %2d avg %.3f sec min %.3f, max %.3f, std dev %.3f\n", i,
                 time_avg[i], time_lo[i], time_hi[i], time_std[i]);
      }


      puts(" ");
      for(i=0; i<NSAV; i++) {
        time_lo[i]= 1.0e40;
        time_hi[i]= -1.0e40;
        time_sum[i]= 0.0;
        time_sum2[i]= 0.0;
        for(j= 0; j < g3.nproc; j++) {
          timbase= time_fulldata[i*g3.nproc+j]+time_fulldata2[i*g3.nproc+j];
          time_lo[i]= min(time_lo[i], timbase);
          time_hi[i]= max(time_hi[i], timbase);
          time_sum[i] += timbase;
          time_sum2[i] += timbase*timbase;
        }
        time_avg[i]= time_sum[i]/g3.nproc;
        time_std[i]= sqrt((time_sum2[i]-g3.nproc*time_avg[i]*time_avg[i]) / (g3.nproc-1));
      }
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\t2D-fft pass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 (time_one[i]+time_two[i]), 1.0e-9*(bytes_one[i]+bytes_two[i]) / (time_one[i]+time_two[i]) );
      }
      puts(" ");
      for(i=0; i<NSAV; i++) {
          printf("\t2D-fft pass %2d avg %.3f sec min %.3f, max %.3f, std dev %.3f\n", i,
                 time_avg[i], time_lo[i], time_hi[i], time_std[i]);
      }
      print_full_fft(NSAV, time_one, bytes_one, time_two, bytes_two, "2D FFT message passing (Alltoall)", 2);
#if 0
      print_full_fft(NSAV, time_fulldata, bytes_one, time_fulldata2, bytes_two, "2D FFT message passing (Alltoall)", 2);
#endif
      fom_2dfft= fom;
    }
  }


  
  /* ------------------------------ start the test of passing planes
     in the z-direction as in the light advection algorithm in pF3D */
  if( run_plane) {
    print_time("START time is");
    MPI_Barrier(g3.gridcom);
    for (i=0; i<nTrip; i++) {
	startTime = MPI_Wtime();

        syncplane(arry, nxl, nyl, NPLANE);
        
        stopTime = MPI_Wtime();
        if(i < NSAV) {
          time_one[i]= stopTime - startTime;
          bytes_one[i]= bytes_moved;
        }
    }
    MPI_Barrier(g3.gridcom);
    print_time("STOP time is");
    /* Gather the timing data from all processes back to rank zero.
       Use an Allgather, even though only rank zero needs the data,
       to keep the code uniform.
    */
    ival= MPI_Allgather(time_one, NSAV, MPI_DOUBLE, time_fulldata, NSAV, MPI_DOUBLE, g3.gridcom);
    if ( g3.me == master ) {
      /* Compute min, max, avg, and standard deviation of time over
         all MPI processes. */
      for(i=0; i<NSAV; i++) {
        time_lo[i]= 1.0e40;
        time_hi[i]= -1.0e40;
        time_sum[i]= 0.0;
        time_sum2[i]= 0.0;
        for(j= 0; j < g3.nproc; j++) {
          timbase= time_fulldata[i*g3.nproc+j];
          time_lo[i]= min(time_lo[i], timbase);
          time_hi[i]= max(time_hi[i], timbase);
          time_sum[i] += timbase;
          time_sum2[i] += timbase*timbase;
        }
        time_avg[i]= time_sum[i]/g3.nproc;
        time_std[i]= sqrt((time_sum2[i]-g3.nproc*time_avg[i]*time_avg[i]) / (g3.nproc-1));
      }
      puts("\n SyncForward");
      for(i=0; i<NSAV; i++) {
          /* bytes_moved for this pass. */
          printf("\tsyncforward pass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
      puts(" ");
      for(i=0; i<NSAV; i++) {
          printf("\tsyncforward test %2d avg %.3f sec min %.3f, max %.3f, std dev %.3f\n", i,
                 time_avg[i], time_lo[i], time_hi[i], time_std[i]);
      }
      print_one(NSAV, time_one, bytes_one, "Pass 3 xy-planes to neighbor", 1);
      fom_fwd= fom;   
    }
  }
    
    
    /* ------------------------------ start the test of passing planes
       as in the hydro package in pF3D */
  if( run_hydro) {
    print_time("START time is");
    MPI_Barrier(g3.gridcom);
    for (i=0; i<nTrip; i++) {
	startTime = MPI_Wtime();

        syncplane_hydro(arry, nxl, nyl, nzl);
        
        stopTime = MPI_Wtime();
        if(i < NSAV) {
          time_one[i]= stopTime - startTime;
          bytes_one[i]= bytes_moved;
        }
    }
    MPI_Barrier(g3.gridcom);
    print_time("STOP time is");
    /* Gather the timing data from all processes back to rank zero.
       Use an Allgather, even though only rank zero needs the data,
       to keep the code uniform.
    */
    ival= MPI_Allgather(time_one, NSAV, MPI_DOUBLE, time_fulldata, NSAV, MPI_DOUBLE, g3.gridcom);
    if ( g3.me == master ) {
      /* Compute min, max, avg, and standard deviation of time over
         all MPI processes. */
      for(i=0; i<NSAV; i++) {
        time_lo[i]= 1.0e40;
        time_hi[i]= -1.0e40;
        time_sum[i]= 0.0;
        time_sum2[i]= 0.0;
        for(j= 0; j < g3.nproc; j++) {
          timbase= time_fulldata[i*g3.nproc+j];
          time_lo[i]= min(time_lo[i], timbase);
          time_hi[i]= max(time_hi[i], timbase);
          time_sum[i] += timbase;
          time_sum2[i] += timbase*timbase;
        }
        time_avg[i]= time_sum[i]/g3.nproc;
        time_std[i]= sqrt((time_sum2[i]-g3.nproc*time_avg[i]*time_avg[i]) / (g3.nproc-1));
      }
      puts("\nHydro boundary exchange");
      for(i=0; i<NSAV; i++) {
          /* bytes_moved is for this pass. */
          printf("\tsyncplane_hydro pass %2d took %.3f sec (B/W %.3f GB/sec)\n", i,
                 time_one[i], 1.0e-9 * bytes_one[i] / time_one[i]);
      }
      puts(" ");
      for(i=0; i<NSAV; i++) {
          printf("\tsyncplane_hydro test %2d avg %.3f sec min %.3f, max %.3f, std dev %.3f\n", i,
                 time_avg[i], time_lo[i], time_hi[i], time_std[i]);
      }
      print_one(NSAV, time_one, bytes_one, "Pass hydro guard planes", 1);
    }
  }
    
    
    /* ------------------------------ release storage no longer needed */
    free(brry);
    brry= 0;
    free(arry);
    arry= 0;
    if ( g3.me == master ) {
        /* print summary information for key tests */
        printf("SUMMARY - decomp (%d,%d,%d) has 2D FFT rate %e GB/s and advection rate %e GB/s\n", mp_p, mp_q, mp_r, fom_2dfft, fom_fwd);
    }

    MPI_Barrier(g3.gridcom);
    time_finish= MPI_Wtime();
    time_total= time_finish-time_begin;
    if(!mp_rank)printf("Total elapsed time is %f\n", time_total);
    
    if(do_csv) {
      fclose(file_csv);
    }

    /* ------------------------------ Finalize MPI for clean exit */
    simdone();
}

void print_one(int nsav, double *time_one, double *bytes_one, char *label, int dir)
{
    double sum, sumsq, avgTime, errTime, avgRate, bytes_moved, fac, facx, facy;
    long i, ntripm;
  
    /* Note that bytes_one is the count for one pass. That matches well with
       using the average time. */
    if(nsav < 1) return;
    /* ------------------------------ compute the timing statistics */
    sum= sumsq= 0.0;
    for(i= 0; i < nsav; i++) {
        sum += time_one[i];
        sumsq += sqr(time_one[i]);
    }
    bytes_moved= bytes_one[0];
    avgTime = sum / (double) nsav;
    errTime = sqrt( (nsav * sumsq - sqr(sum)) / (sqr(nsav) * (nsav - 1.0)) );
    avgRate= 1.0e-9 * bytes_moved / avgTime;
    printf("\nAll tests\n\t%s in %.3f +/- %.3f sec (B/W %.3f GB/sec)\n", label,
    	   avgTime, errTime, avgRate);
    printf("\tMessage passing phases %.3f per sec, %ld total bytes moved\n\n",
           1.0 / avgTime,  (long) bytes_moved);
    if(nsav < 2) return;
    if(do_csv) {
      fprintf(file_csv, "%e, %e, %e, ", avgTime, errTime, avgRate);
    }

    /* --------- compute the timing statistics without the first pass */
    ntripm= nsav-1;
    sum= sumsq= 0.0;
    for(i= 1; i < nsav; i++) {
        sum += time_one[i];
        sumsq += sqr(time_one[i]);
    }
    avgTime = sum / (double) ntripm;
    errTime = sqrt( (ntripm * sumsq - sqr(sum)) / (sqr(ntripm) * (ntripm - 1.0)) );
    avgRate= 1.0e-9 * bytes_moved / avgTime;
    printf("\nAll but first test\n\t%s in %.3f +/- %.3f sec (B/W %.3f GB/sec)\n", label,
    	   avgTime, errTime, avgRate);
    printf("\tMessage passing phases %.3f per sec, %ld total bytes moved\n\n",
           1.0 / avgTime,  (long) bytes_moved);
    if(dir > 0) {
      /* passing in the X direction */
      if(mp_p > 1) fac= (mp_p-1.0)/mp_p;
      else fac= 1.0;
    } else if(dir < 0) {
      /* passing in the Y direction */
      if(mp_q > 1) fac= (mp_q-1.0)/mp_q;
      else fac= 1.0;
    } else {
      /* passing in both the x and y directions */
      if(mp_p > 1) facx= (mp_p-1.0)/mp_p;
      else facx= 1.0;
      if(mp_q > 1) facy= (mp_q-1.0)/mp_q;
      else facy= 1.0;
      /* Warning - this is not correct if FFT messages
         in x and y have different lengths */
      fac= (facx+facy)*0.5;
    }
    /* Print the average number of GB per second sent by one process. */
    fom= 1.0e-9*fac*bytes_moved/avgTime;
    printf("\tMessage rate per process is %f GB/s\n\n", fom);
}


void print_full_fft(int nsav, double *time_one, double *bytes_one, double *time_two, double *bytes_two, char *label, int dir)
{
  double sum, sumsq, avgTime, errTime, avgRate, bytes_sent, fac, facx, facy;
    double tim2;
    long i, ntripm;
  
    /* Note that bytes_one is the count for one pass. That matches well with
       using the average time. */
    if(nsav < 1) return;
    /* ------------------------------ compute the timing statistics */
    sum= sumsq= 0.0;
    bytes_sent= bytes_one[0]+bytes_two[0];
    for(i= 0; i < nsav; i++) {
        tim2 = time_one[i]+time_two[i];
        sum += tim2;
        sumsq += tim2*tim2;
    }
    printf("sum of times is %e\n", sum);
    avgTime = sum / (double) nsav;
    errTime = sqrt( (nsav * sumsq - sqr(sum)) / (sqr(nsav) * (nsav - 1.0)) );
    avgRate= 1.0e-9 * bytes_sent / avgTime;
    printf("\nAll tests\n\t%s in %.3f +/- %.3f sec (B/W %.3f GB/sec)\n", label,
    	   avgTime, errTime, avgRate);
    printf("\tMessage passing phases %.3f per sec, %ld total bytes moved\n\n",
           1.0 / avgTime,  (long) bytes_sent);
    if(nsav < 2) return;
    if(do_csv) {
      fprintf(file_csv, "%e, %e, %e, ", avgTime, errTime, avgRate);
    }

    /* --------- compute the timing statistics without the first pass */
    ntripm= nsav-1;
    sum= sumsq= 0.0;
    for(i= 1; i < nsav; i++) {
        tim2 = time_one[i]+time_two[i];
        sum += tim2;
        sumsq += tim2*tim2;
    }
    avgTime = sum / (double) ntripm;
    errTime = sqrt( (ntripm * sumsq - sqr(sum)) / (sqr(ntripm) * (ntripm - 1.0)) );
    avgRate= 1.0e-9 * bytes_sent / avgTime;
    printf("\nAll but first test\n\t%s in %.3f +/- %.3f sec (B/W %.3f GB/sec)\n", label,
    	   avgTime, errTime, avgRate);
    printf("\tMessage passing phases %.3f per sec, %ld total bytes moved\n\n",
           1.0 / avgTime,  (long) bytes_sent);
    if(dir > 0) {
      /* passing in the X direction */
      if(mp_p > 1) fac= (mp_p-1.0)/mp_p;
      else fac= 1.0;
    } else if(dir < 0) {
      /* passing in the Y direction */
      if(mp_q > 1) fac= (mp_q-1.0)/mp_q;
      else fac= 1.0;
    } else {
      /* passing in both the x and y directions */
      if(mp_p > 1) facx= (mp_p-1.0)/mp_p;
      else facx= 1.0;
      if(mp_q > 1) facy= (mp_q-1.0)/mp_q;
      else facy= 1.0;
      /* Warning - this is not correct if FFT messages
         in x and y have different lengths */
      fac= (facx+facy)*0.5;
    }
    /* Print the average number of GB per second sent by one process. */
    fom= 1.0e-9*fac*bytes_sent/avgTime;
    printf("\tMessage rate per process is %f GB/s\n\n", fom);
}

